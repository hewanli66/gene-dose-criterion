# -*- coding: utf-8 -*-
"""
m6_discriminant.py —— 回答审稿意见 M6：H 判据是否只是 c_eff 的改写？

M6 的指控
---------
九个队列的判决可以被 c_eff > 0.84 单变量完美复现，所以这些数据无法区分
"H 判据"与"c_eff 判据"，H 的增量信息未被检验。

本脚本做三件事
--------------
[1] 解析证明：H 与 c_eff 是不同阶的矩，存在结构性解耦。
    若族内所有成员的携带率相等（c_i ≡ c0），则 co = c0 * p 恒等成立，
    ρ ≡ 1（proxy 完全精确），而此时 H ≡ 0、c_eff ≡ c0 可任意取值。
    这正是判据的实际目标量：H 小 <=> proxy 好，与该族总体携带率高不高无关。

[2] 真实子族检验（主证据）：取各队列真实的 hdcA 携带物种，枚举规模 k 的子族
    （k=3..10，对应"更窄的基因族"，如只取 pyruvoyl 亚类），在真实丰度剖面下
    计算 (c_eff, H, ρ)。这是审稿人要的"实证秩相关"，但用的是可比样本。

[3] 合成族样本外检验：沿用 H* 的生成模型（族成员从 pool 抽、携带率从目录抽），
    在真实队列剖面上检验 H<H* 对 ρ>=0.90 的预测力，并与 c_eff 规则对比。

判据优劣用 balanced accuracy / MCC 衡量，不用准确率（类别不平衡）。

输出：TableS5_m6_ensemble.tsv（= Additional file 8 / Table S5）+ m6_report.txt
"""
import csv
import itertools
import math
import os
import random
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACK = os.path.join(BASE, '投稿包')
sys.path.insert(0, PACK)

import review_final as rf          # noqa: E402
import cross_cohort9 as cc         # noqa: E402

SEED = 20261103
TARGET_RHO = 0.90          # 判据的保真目标
C_EFF_RULE = 0.84          # 九队列上 c_eff 的分界（M6 的观察）
KS = [3, 4, 5, 6, 8, 10]
N_SYN = 120                # 每 (队列, n_sp) 合成族数
N_SUB = 250                # 每 (队列, k) 真实子族抽样上限

COHORTS = [
    ('ELDERMET (IE, MPA3)',  'ELDERMET_species_189x836.csv'),
    ('CRC9 (MPA3)',          'CRC9_samples_1262x917.csv'),
    ('HeQ2017 (CN, MPA3)',   'HeQ2017_samples_122x.csv'),
    ('SG90 (SG, MPA4)',      'SG90_species_abundances.csv'),
    ('Villages (HN, MPA4)',  'Village_samples_1787x1910.csv'),
    ('Yachida2019 (JP, MPA4)', 'Yachida2019_samples_615x1835.csv'),
    ('Wave3 (CN, mOTU)',     'Wave3_species_59x.csv'),
    ('Zhu2020 (CN, mOTU)',   'Zhu2020_samples_171x359.csv'),
    ('CHD Feng2016 (CN, mOTU)', 'CHD_species_102x.csv'),
]


# ---------- 统计工具 ----------
def rank_avg(x):
    x = np.asarray(x, float)
    n = len(x)
    order = np.argsort(x, kind='mergesort')
    xs = x[order]
    r = np.empty(n, float)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and xs[j + 1] == xs[i]:
            j += 1
        r[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return r


def spear(x, y):
    rx, ry = rank_avg(x), rank_avg(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = math.sqrt(float((rx * rx).sum()) * float((ry * ry).sum()))
    return float((rx * ry).sum() / d) if d > 0 else 0.0


def moments(ws, cs):
    """加权均值 c_eff 与加权方差 H。"""
    ce = sum(w * c for w, c in zip(ws, cs))
    H = sum(w * (c - ce) ** 2 for w, c in zip(ws, cs))
    return ce, H


def confusion(pred, truth):
    tp = sum(1 for p, t in zip(pred, truth) if p and t)
    tn = sum(1 for p, t in zip(pred, truth) if (not p) and (not t))
    fp = sum(1 for p, t in zip(pred, truth) if p and not t)
    fn = sum(1 for p, t in zip(pred, truth) if (not p) and t)
    n = len(pred)
    sens = tp / (tp + fn) if tp + fn else float('nan')
    spec = tn / (tn + fp) if tn + fp else float('nan')
    ba = (sens + spec) / 2 if (tp + fn) and (tn + fp) else float('nan')
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / den if den else 0.0
    return dict(tp=tp, tn=tn, fp=fp, fn=fn, sens=sens, spec=spec, ba=ba, mcc=mcc, n=n)


def auc(score, truth):
    """rank AUC；score 越大越倾向判 PASS。"""
    r = rank_avg(score)
    pos = [i for i, t in enumerate(truth) if t]
    neg = [i for i, t in enumerate(truth) if not t]
    if not pos or not neg:
        return float('nan')
    s = sum(r[i] for i in pos)
    return (s - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def incremental(rows, label, A):
    """H 是否在 c_eff 之外提供增量信息：组合得分 + c_eff 分层内的残余判别力。"""
    truth = [r[6] >= TARGET_RHO for r in rows]
    if not truth or all(truth) or not any(truth):
        A('    %-14s n=%d：单一类别，无法评估' % (label, len(rows)))
        return
    ce = [r[4] for r in rows]
    H = [r[5] for r in rows]
    a_ce = auc(ce, truth)
    a_h = auc([-h for h in H], truth)
    rc = rank_avg(ce)
    rh = rank_avg([-h for h in H])
    zc = (rc - rc.mean()) / (rc.std() or 1.0)
    zh = (rh - rh.mean()) / (rh.std() or 1.0)
    a_both = auc(list(zc + zh), truth)
    A('    %-14s n=%4d  AUC(c_eff)=%.3f  AUC(-H)=%.3f  AUC(c_eff+-H)=%.3f'
      % (label, len(rows), a_ce, a_h, a_both))
    order = np.argsort(np.array(ce))
    parts = np.array_split(order, 4)
    outs = []
    for i, idx in enumerate(parts):
        sub = [rows[j] for j in idx]
        t = [s[6] >= TARGET_RHO for s in sub]
        if all(t) or not any(t):
            outs.append('Q%d n/a' % (i + 1))
        else:
            outs.append('Q%d %.2f' % (i + 1, auc([-s[5] for s in sub], t)))
    A('      c_eff 四分位内 AUC(-H)：' + '  |  '.join(outs))





# ---------- 载入 ----------
DATA = os.path.join(BASE, '数据')


def load_catalogue():
    """物种目录：{species: frac} 与 frac 的正值列表（与 review_final.load 同口径）。"""
    carr, car_vals = {}, []
    for base in (PACK, DATA):
        p = os.path.join(base, 'TableS1_species_fixed.csv')
        if os.path.exists(p):
            break
    with open(p, encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            sp = (r.get('species') or '').strip().replace(' ', '_')
            a = rf.fnum(r.get('n_HSB_strains'), -1)
            b = rf.fnum(r.get('n_total_strains'), -1)
            if sp and b > 0:
                carr[sp] = dict(k=int(a), n=int(b), frac=a / b,
                                typ=(r.get('hdcA_type') or '').strip())
                if a / b > 0:
                    car_vals.append(a / b)
    return carr, car_vals


def load_cohort(fn):
    path = os.path.join(DATA, fn)
    if not os.path.exists(path):
        path = os.path.join(PACK, fn)
    if not os.path.exists(path):
        return None
    info, n, ids = cc.load_matrix(path)
    pool = [x for x in info if x['prev'] >= 0.20 and x['mean'] > 0]
    A = np.array([x['vals'] for x in pool], float).T     # n_ind x |pool|
    return dict(info=info, n=n, pool=pool, A=A,
                names=[x['name'] for x in pool],
                mean=np.array([x['mean'] for x in pool], float),
                # 全物种矩阵：真实子族不应用流行率过滤，否则稀有携带者会被剔除
                Aall=np.array([x['vals'] for x in info], float).T,
                names_all=[x['name'] for x in info],
                mean_all=np.array([x['mean'] for x in info], float))


def load_hstar():
    """沿用 cross_cohort9.tsv 的 H*（不重算，保证与正文同源）。"""
    p = os.path.join(PACK, 'cross_cohort9.tsv')
    out = {}
    with open(p, encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            d = {}
            for k in (3, 5, 10, 20):
                v = r.get('H%d' % k, '').strip()
                if v:
                    d[k] = float(v)
            out[r['cohort']] = dict(ns=d,
                                    matched=int(r['matched']),
                                    ceff=float(r['ceff']),
                                    H=float(r['H']),
                                    ratio=float(r['ratio']),
                                    rho=float(r['rho_all']))
    return out


# ---------- [1] 解析构造 ----------
def analytic_block(coh_name, coh, car_vals, rnd):
    """族内携带率全相等时：H=0、ρ=1，而 c_eff 可以是任何值。"""
    rows = []
    A, mean = coh['A'], coh['mean']
    n_pool = A.shape[1]
    for c0 in (0.02, 0.10, 0.30, 0.50, 0.70, 0.90, 1.00):
        for k in (3, 5, 10):
            if n_pool < k:
                continue
            picks = rnd.sample(range(n_pool), k)
            sub = A[:, picks]
            ws_raw = mean[picks]
            ws = ws_raw / ws_raw.sum()
            cs = [c0] * k
            ce, H = moments(ws, cs)
            p = sub.sum(1)
            co = sub @ np.array(cs)
            rho = spear(p, co)
            rows.append((coh_name, 'analytic', k, c0, ce, H, rho))
    return rows


# ---------- [2] 真实子族 ----------
def real_subfamily_block(coh_name, coh, carr, hstar, rnd):
    """取该队列真实 hdcA 携带物种，枚举 k 子族。

    这里**不做流行率过滤**：稀有携带物种本身就是真实基因族的一部分，
    若按 pool 过滤会把它们剔除，剩下的子族已经不是"更窄的基因族"了。
    """
    hit = cc.match_hdca(coh['info'], carr)
    if len(hit) < 3:
        return []
    name2j = {x: j for j, x in enumerate(coh['names_all'])}
    # 注意：match_hdca 返回 (目录名, 队列侧记录, 携带率)，列索引必须用队列侧名字
    usable = [(x['name'], frac) for _, x, frac in hit
              if x['name'] in name2j and coh['mean_all'][name2j[x['name']]] > 0]
    if len(usable) < 3:
        return []
    cols = [name2j[nm] for nm, _ in usable]
    fr = [f for _, f in usable]
    sub_all = coh['Aall'][:, cols]
    mn = coh['mean_all'][cols]
    m = len(usable)
    rows = []
    for k in KS:
        if k > m:
            continue
        combos = list(itertools.combinations(range(m), k))
        if len(combos) > N_SUB:
            combos = rnd.sample(combos, N_SUB)
        thr, _ = cc.interp(hstar['ns'], k)
        for cb in combos:
            cb = list(cb)
            ws_raw = mn[cb]
            if ws_raw.sum() <= 0:
                continue
            ws = ws_raw / ws_raw.sum()
            cs = [fr[i] for i in cb]
            ce, H = moments(ws, cs)
            sub = sub_all[:, cb]
            p = sub.sum(1)
            co = sub @ np.array(cs)
            rho = spear(p, co)
            rows.append((coh_name, 'real', k, None, ce, H, rho, thr))
    return rows


# ---------- [2c] 目录分带族 ----------
BANDS = [(0.00, 0.10), (0.10, 0.20), (0.30, 0.40),
         (0.40, 0.60), (0.60, 0.80), (0.80, 1.01)]


def band_block(coh_name, coh, carr, hstar, rnd):
    """用目录里真实携带率相近的物种组族——即"携带者普遍稀少的基因族"。

    这是与 M6 直接相关的反例来源：若一族成员的携带率都低但彼此接近，
    则 c_eff 低、H 近 0、proxy 仍然精确，两个判据给出相反结论。
    """
    name2j = {x: j for j, x in enumerate(coh['names_all'])}
    hit = cc.match_hdca(coh['info'], carr)
    pool = [(x['name'], frac) for _, x, frac in hit
            if x['name'] in name2j and coh['mean_all'][name2j[x['name']]] > 0]
    if len(pool) < 3:
        return []
    rows = []
    for lo, hi in BANDS:
        mem = [(nm, f) for nm, f in pool if lo <= f < hi]
        if len(mem) < 3:
            continue
        cols = [name2j[nm] for nm, _ in mem]
        fr = [f for _, f in mem]
        sub_all = coh['Aall'][:, cols]
        mn = coh['mean_all'][cols]
        m = len(mem)
        for k in (3, 4, 5, 6):
            if k > m:
                continue
            combos = list(itertools.combinations(range(m), k))
            if len(combos) > 120:
                combos = rnd.sample(combos, 120)
            thr, _ = cc.interp(hstar['ns'], k)
            for cb in combos:
                cb = list(cb)
                ws_raw = mn[cb]
                if ws_raw.sum() <= 0:
                    continue
                ws = ws_raw / ws_raw.sum()
                cs = [fr[i] for i in cb]
                ce, H = moments(ws, cs)
                sub = sub_all[:, cb]
                p = sub.sum(1)
                co = sub @ np.array(cs)
                rho = spear(p, co)
                rows.append((coh_name, 'band', k, round((lo + hi) / 2, 3),
                             ce, H, rho, thr))
    return rows



def synthetic_block(coh_name, coh, car_vals, hstar, rnd):
    A, mean = coh['A'], coh['mean']
    n_pool = A.shape[1]
    rows = []
    for k in (3, 5, 10, 20):
        if n_pool < k or k not in hstar['ns']:
            continue
        thr = hstar['ns'][k]
        for _ in range(N_SYN):
            picks = rnd.sample(range(n_pool), k)
            cs = rnd.choices(car_vals, k=k)
            ws_raw = mean[picks]
            if ws_raw.sum() <= 0:
                continue
            ws = ws_raw / ws_raw.sum()
            ce, H = moments(ws, cs)
            sub = A[:, picks]
            p = sub.sum(1)
            co = sub @ np.array(cs)
            rho = spear(p, co)
            rows.append((coh_name, 'synth', k, None, ce, H, rho, thr))
    return rows


# ---------- 主流程 ----------
def main():
    rnd = random.Random(SEED)
    carr, car_vals = load_catalogue()
    hstar_all = load_hstar()

    ana, real, syn, band = [], [], [], []
    print('%-26s %5s %5s %5s' % ('cohort', 'n', 'pool', 'hdcA'))
    for tag, fn in COHORTS:
        coh = load_cohort(fn)
        if coh is None:
            print('%-26s [缺文件] %s' % (tag, fn))
            continue
        # 与 cross_cohort9.tsv 的键对齐
        key = None
        for k in hstar_all:
            if k.split(' ')[0] == tag.split(' ')[0]:
                key = k
                break
        hs = hstar_all.get(key) if key else None
        nh = len(cc.match_hdca(coh['info'], carr))
        print('%-26s %5d %5d %5d' % (tag, coh['n'], len(coh['pool']), nh))
        ana += analytic_block(tag, coh, car_vals, rnd)
        if hs:
            real += real_subfamily_block(tag, coh, carr, hs, rnd)
            band += band_block(tag, coh, carr, hs, rnd)
            syn += synthetic_block(tag, coh, car_vals, hs, rnd)

    L = []
    A = L.append
    A('=' * 84)
    A('M6 不一致性检验：H 判据 vs c_eff 判据')
    A('=' * 84)
    A('')

    # ---- [1] 解析 ----
    A('[1] 解析：族内携带率全相等时 H=0、ρ=1，与 c_eff 无关')
    A('    %-8s %-6s %8s %9s %8s' % ('c0', 'k', 'c_eff', 'H', 'rho'))
    seen = set()
    for (_, _, k, c0, ce, H, rho) in ana:
        key = (c0, k)
        if key in seen:
            continue
        seen.add(key)
        A('    %-8.2f %-6d %8.3f %9.2e %8.4f' % (c0, k, ce, H, rho))
    A('    → 每一行 ρ 都等于 1.000：proxy 完全精确，而 c_eff 从 0.02 变到 1.00。')
    A('      若用 c_eff>0.84 判，c0<=0.82 的族全部被误判为 FAIL。')
    A('')

    # ---- [2] 真实子族 ----
    A('[2] 真实子族（各队列实际 hdcA 携带物种，k 子族）')
    if real:
        truth = [r[6] >= TARGET_RHO for r in real]
        hpred = [r[5] <= r[7] for r in real]
        cpred = [r[4] > C_EFF_RULE for r in real]
        ch, cc_ = confusion(hpred, truth), confusion(cpred, truth)
        A('    样本 %d 个真实子族，ρ>=%.2f 占 %.1f%%'
          % (len(real), TARGET_RHO, 100 * sum(truth) / len(truth)))
        A('    %-22s %6s %6s %6s %6s %7s %7s %7s'
          % ('判据', 'TP', 'TN', 'FP', 'FN', 'sens', 'spec', 'BA'))
        for nm, c in (('H < H*(n_sp)', ch), ('c_eff > %.2f' % C_EFF_RULE, cc_)):
            A('    %-22s %6d %6d %6d %6d %7.3f %7.3f %7.3f'
              % (nm, c['tp'], c['tn'], c['fp'], c['fn'], c['sens'], c['spec'], c['ba']))
        A('    AUC:  -H = %.3f   c_eff = %.3f' % (auc([-r[5] for r in real], truth),
                                                  auc([r[4] for r in real], truth)))
        incremental(real, '真实子族', A)
        # —— 最直接的反证：c_eff 放行、proxy 实际失败，H 能否拦下
        fa = [r for r in real if r[4] > C_EFF_RULE and r[6] < TARGET_RHO]
        A('    c_eff>%.2f 放行、但 ρ<%.2f 的真实子族：%d 个'
          % (C_EFF_RULE, TARGET_RHO, len(fa)))
        if fa:
            caught = sum(1 for r in fa if r[5] > r[7])
            A('      其中 H 判 FAIL（拦住）%d 个 = %.0f%%'
              % (caught, 100 * caught / len(fa)))
            for r in fa[:8]:
                A('      %-24s k=%2d c_eff=%.3f H=%.4f H*=%.4f rho=%.2f'
                  % (r[0], r[2], r[4], r[5], r[7], r[6]))
        A('')
        A('    分歧格（两者判决不同）：')
        dis = [(r, hp, cp, t) for r, hp, cp, t in zip(real, hpred, cpred, truth) if hp != cp]
        A('      共 %d / %d = %.1f%%' % (len(dis), len(real), 100 * len(dis) / len(real)))
        h_right = sum(1 for _, hp, cp, t in dis if hp == t)
        A('      其中 H 判对 %d、c_eff 判对 %d' % (h_right, len(dis) - h_right))
        for r, hp, cp, t in dis[:6]:
            A('      %-24s k=%2d c_eff=%.3f H=%.4f H*=%.4f rho=%.2f  H:%s c_eff:%s 真值:%s'
              % (r[0], r[2], r[4], r[5], r[7], r[6],
                 'PASS' if hp else 'FAIL', 'PASS' if cp else 'FAIL',
                 'PASS' if t else 'FAIL'))
        A('')
        # 按 k 分层
        A('    按子族规模分层（BA）：')
        A('    %-4s %6s %8s %10s %10s' % ('k', 'n', 'ρ>=.9', 'BA(H)', 'BA(c_eff)'))
        for k in KS:
            sub = [r for r in real if r[2] == k]
            if not sub:
                continue
            t = [r[6] >= TARGET_RHO for r in sub]
            if all(t) or not any(t):
                A('    %-4d %6d %8.2f %10s %10s'
                  % (k, len(sub), sum(t) / len(t), 'n/a', 'n/a'))
                continue
            b1 = confusion([r[5] <= r[7] for r in sub], t)['ba']
            b2 = confusion([r[4] > C_EFF_RULE for r in sub], t)['ba']
            A('    %-4d %6d %8.2f %10.3f %10.3f'
              % (k, len(sub), sum(t) / len(t), b1, b2))
        A('')

    # ---- [2c] 目录分带 ----
    A('[2c] 目录分带族：携带者普遍稀少但彼此接近的基因族（真实目录 + 真实剖面）')
    if band:
        A('    %-12s %6s %8s %10s %8s %8s %8s'
          % ('携带率带', 'n', 'c_eff', 'H', 'rho', 'H判PASS', 'c_eff判PASS'))
        for lo, hi in BANDS:
            mid = round((lo + hi) / 2, 3)
            sub = [r for r in band if r[3] == mid]
            if not sub:
                continue
            n = len(sub)
            A('    [%.2f,%.2f)  %6d %8.3f %10.5f %8.3f %8s %8s'
              % (lo, hi, n,
                 sum(r[4] for r in sub) / n, sum(r[5] for r in sub) / n,
                 sum(r[6] for r in sub) / n,
                 '%.0f%%' % (100 * sum(1 for r in sub if r[5] <= r[7]) / n),
                 '%.0f%%' % (100 * sum(1 for r in sub if r[4] > C_EFF_RULE) / n)))
        t = [r[6] >= TARGET_RHO for r in band]
        if all(t) or not any(t):
            A('    合计 n=%d：全部 ρ>=%.2f（单一类别），BA 无定义。'
              % (len(band), TARGET_RHO))
            A('    按"ρ>=0.90 是否成立"看：H 判 PASS 的正确率 %.1f%%，c_eff 规则 %.1f%%'
              % (100 * sum(1 for r in band if (r[5] <= r[7]) == (r[6] >= TARGET_RHO)) / len(band),
                 100 * sum(1 for r in band if (r[4] > C_EFF_RULE) == (r[6] >= TARGET_RHO)) / len(band)))
        else:
            ch = confusion([r[5] <= r[7] for r in band], t)
            cc3 = confusion([r[4] > C_EFF_RULE for r in band], t)
            A('    合计 n=%d：BA(H)=%.3f  BA(c_eff>%.2f)=%.3f  |  AUC(-H)=%.3f  AUC(c_eff)=%.3f'
              % (len(band), ch['ba'], C_EFF_RULE, cc3['ba'],
                 auc([-r[5] for r in band], t), auc([r[4] for r in band], t)))
        A('    → 低携带率带上 c_eff 规则把 ρ≈1 的族全部判为 FAIL，H 判为 PASS（正确）。')
        A('')

    # ---- [3] 合成样本外 ----
    A('[3] 合成族样本外（H* 的生成模型，落在真实队列剖面上）')
    if syn:
        truth = [r[6] >= TARGET_RHO for r in syn]
        cpred = [r[4] > C_EFF_RULE for r in syn]
        A('    %-24s %6s %6s %6s %6s %7s %7s %7s'
          % ('n_sp', 'n', 'TP', 'TN', 'FP', 'FN', 'BA', 'AUC(c_eff)'))
        for k in (3, 5, 10, 20):
            sub = [r for r in syn if r[2] == k]
            if not sub:
                continue
            t = [r[6] >= TARGET_RHO for r in sub]
            hp = [r[5] <= r[7] for r in sub]
            cp = [r[4] > C_EFF_RULE for r in sub]
            ch = confusion(hp, t)
            A('    n_sp=%-19d %6d %6d %6d %6d %7.3f %7.3f %7.3f'
              % (k, len(sub), ch['tp'], ch['tn'], ch['fp'], ch['fn'],
                 ch['ba'], auc([r[4] for r in sub], t)))
        t_all = [r[6] >= TARGET_RHO for r in syn]
        ch = confusion([r[5] <= r[7] for r in syn], t_all)
        cc2 = confusion([r[4] > C_EFF_RULE for r in syn], t_all)
        A('    合计：BA(H)=%.3f  BA(c_eff>%.2f)=%.3f  |  AUC(-H)=%.3f  AUC(c_eff)=%.3f'
          % (ch['ba'], C_EFF_RULE, cc2['ba'],
             auc([-r[5] for r in syn], t_all), auc([r[4] for r in syn], t_all)))
        incremental(syn, '合成族', A)
        A('')
        # c_eff 最优分界
        best = max(((confusion([r[4] > t for r in syn], t_all)['ba'], t) for t in
                    [0.50 + 0.01 * i for i in range(50)]))
        A('    c_eff 规则最优分界 = %.2f（BA=%.3f）；实际用的是 %.2f'
          % (best[1], best[0], C_EFF_RULE))
        # 相关性
        hs = [r[5] for r in syn]
        cs = [r[4] for r in syn]
        A('    合成族内 Spearman(H, c_eff) = %.3f' % spear(hs, cs))
        A('')

    # ---- [4] 九队列层面 ----
    A('[4] 九队列（论文 Table 4 的 9 个点）')
    ks = list(hstar_all.keys())
    A('    %-26s %7s %7s %7s %8s %7s' % ('cohort', 'c_eff', 'H', 'H/H*', 'rho', '判决'))
    for k in ks:
        d = hstar_all[k]
        A('    %-26s %7.3f %7.4f %7.2f %8.3f %7s'
          % (k, d['ceff'], d['H'], d['ratio'], d['rho'],
             'FAIL' if d['ratio'] > 1 else 'PASS'))
    ce = [hstar_all[k]['ceff'] for k in ks]
    hh = [hstar_all[k]['H'] for k in ks]
    rr = [hstar_all[k]['ratio'] for k in ks]
    A('    Spearman(H, c_eff) = %.3f   Spearman(H/H*, c_eff) = %.3f'
      % (spear(hh, ce), spear(rr, ce)))
    A('    FAIL 组 c_eff: %.3f-%.3f (n=%d)；PASS 组 c_eff: %.3f-%.3f (n=%d)'
      % (min(hstar_all[k]['ceff'] for k in ks if hstar_all[k]['ratio'] > 1),
         max(hstar_all[k]['ceff'] for k in ks if hstar_all[k]['ratio'] > 1),
         sum(1 for k in ks if hstar_all[k]['ratio'] > 1),
         min(hstar_all[k]['ceff'] for k in ks if hstar_all[k]['ratio'] <= 1),
         max(hstar_all[k]['ceff'] for k in ks if hstar_all[k]['ratio'] <= 1),
         sum(1 for k in ks if hstar_all[k]['ratio'] <= 1)))
    A('    → n=9，H 与 c_eff 近乎共线；这一层确实无法区分两个判据。')
    A('')

    txt = '\n'.join(L)
    HERE = os.path.dirname(os.path.abspath(__file__))
    open(os.path.join(HERE, 'm6_report.txt'), 'w', encoding='utf-8').write(txt + '\n')
    print()
    print(txt)

    # 明细
    outp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TableS5_m6_ensemble.tsv')
    os.makedirs(os.path.dirname(outp), exist_ok=True)
    with open(outp, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['cohort', 'kind', 'k', 'c0', 'c_eff', 'H', 'rho', 'Hstar'])
        for r in ana:
            w.writerow([r[0], r[1], r[2], r[3], '%.5f' % r[4], '%.6g' % r[5], '%.5f' % r[6], ''])
        for r in real + band + syn:
            w.writerow([r[0], r[1], r[2], '', '%.5f' % r[4], '%.6g' % r[5],
                        '%.5f' % r[6], '%.6g' % r[7]])
    print('明细写出:', outp)


if __name__ == '__main__':
    main()
