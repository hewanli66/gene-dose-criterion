# -*- coding: utf-8 -*-
"""
cross_cohort3.py —— 第三个独立队列（中国长寿带，Wave 3）上的 H* 重算 + 样本量匹配检验。

队列（三个，互相独立）：
  A. ELDERMET   爱尔兰，189 例，836 物种（MetaPhlAn，eLife 50240 supp2）
  B. SG90       新加坡华人，516 例 / ≥77 岁 221 例（MetaPhlAn，Nat Commun 15:7751 supp）
  C. Wave3      中国长寿带，59 例（65–103 岁；mOTU 目录）—— 本文件新增

三个队列用了三套不同的物种定量流程（MetaPhlAn3 / MetaPhlAn4 / mOTU），
所以 H* 若仍落在同一量级、且 hdcA 的判定一致，就不能用"流程特异的假象"解释。

**本脚本额外做一件事**：样本量匹配。
H* 的条件概率判据是 P(ρ ≥ 0.90 | H ≈ h) ≥ 0.95，而 ρ 是"族内跨个体"的秩相关。
个体数越多，ρ 的抽样噪声越小，同一个 H 下越容易满足判据 —— 也就是 H* 本身会随
n_ind 上升。三个队列 n = 189 / 516 / 59，直接比数值会把"样本量"读成"群落结构"。
因此把三个队列都无放回抽到公共最小 n（59），重复若干次，报告 H* 的分布。

运行环境：Python 3.8+，仅标准库。
输入：ELDERMET_species_189x836.csv、SG90_species_abundances.csv、SG90_metadata.csv、
      Wave3_species_59x.csv、TableS1_species_fixed.csv
输出：cross_cohort3.txt、cross_cohort3.tsv
"""
import csv
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import review_final as rf

BASE = os.path.dirname(os.path.abspath(__file__))
NS_LIST = [3, 5, 10, 20]
N_DRAW = 12          # 匹配样本量时的重复抽取次数
SEEDS = (20261029, 7, 13, 29)


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return 0.0


def load_matrix(path, keep_ids=None):
    """样本 × 物种 的 CSV（首列 SampleID）。返回 (info, n_ind, ids)。"""
    with open(path, encoding='utf-8-sig') as f:
        rd = csv.reader(f)
        hdr = next(rd)
        cols = hdr[1:]
        rows = [r for r in rd if r and r[0].strip() and
                (keep_ids is None or r[0].strip() in keep_ids)]
    n = len(rows)
    info = []
    for j, c in enumerate(cols):
        v = [fnum(r[j + 1]) for r in rows]
        info.append(dict(name=c, mean=sum(v) / n,
                         prev=sum(1 for x in v if x > 0) / n, vals=v))
    return info, n, [r[0].strip() for r in rows]


def subset(info, idx):
    """按个体下标取子集，重建 info。"""
    out = []
    n = len(idx)
    for x in info:
        v = [x['vals'][i] for i in idx]
        out.append(dict(name=x['name'], mean=sum(v) / n,
                        prev=sum(1 for y in v if y > 0) / n, vals=v))
    return out


def summarise(tag, info, n_ind, car_vals, quiet=False):
    pool = [x for x in info if x['prev'] >= 0.20 and x['mean'] > 0]
    out = dict(tag=tag, n_ind=n_ind, pool=len(pool), ns={}, cum={})
    if not quiet:
        print('%-30s n_ind=%4d  species=%4d  pool(prev>=20%%)=%3d'
              % (tag, n_ind, len(info), len(pool)))
    for n_sp in NS_LIST:
        if len(pool) < n_sp:
            continue
        fam = rf.families(car_vals, pool, n_sp, n_ind)
        out['ns'][n_sp] = rf.h_cond(fam)
        out['cum'][n_sp] = rf.h_cum(fam)
        if not quiet:
            print('   n_sp=%2d  families=%4d  H*_cond=%.5f  H*_cum=%.5f'
                  % (n_sp, len(fam), out['ns'][n_sp], out['cum'][n_sp]))
    return out


def match_hdca(info, carr):
    by_key = {x['name'].lower().replace('_', ' '): x for x in info}
    hit = []
    for sp, d in carr.items():
        k = sp.lower().replace('_', ' ')
        if k in by_key:
            hit.append((sp, by_key[k], d['frac']))
    return hit


def interp(ns, n_sp):
    """在 log H* vs log n_sp 上线性插值；返回 (阈值, 是否查表值)。"""
    import math
    ks = sorted(ns)
    if n_sp in ns:
        return ns[n_sp], True
    if n_sp <= ks[0]:
        return ns[ks[0]], False      # 外推下限：取端点
    if n_sp >= ks[-1]:
        return ns[ks[-1]], False     # 外推上限：取端点
    lo = max(k for k in ks if k <= n_sp)
    hi = min(k for k in ks if k >= n_sp)
    t = (math.log(n_sp) - math.log(lo)) / (math.log(hi) - math.log(lo))
    return math.exp(math.log(ns[lo]) + t * (math.log(ns[hi]) - math.log(ns[lo]))), False


def hdca_state(tag, info, carr, res):
    hit = match_hdca(info, carr)
    if not hit:
        print('%-30s hdcA: 无匹配' % tag)
        return None
    tot = sum(x['mean'] for _, x, _ in hit) or 1.0
    ws = [x['mean'] / tot for _, x, _ in hit]
    cs = [c for _, _, c in hit]
    c_eff = sum(w * c for w, c in zip(ws, cs))
    H = sum(w * (c - c_eff) ** 2 for w, c in zip(ws, cs))
    thr, exact = interp(res['ns'], len(hit))
    ratio = H / thr
    nz = sum(1 for _, x, _ in hit if x['mean'] > 0)
    print('%-30s hdcA: 命中=%2d 非零=%-2d c_eff=%.3f H=%.4f H*(%d%s)=%.5f H/H*=%.2f'
          % (tag, len(hit), nz, c_eff, H, len(hit),
             '' if exact else ',插值', thr, ratio))
    return dict(tag=tag, matched=len(hit), nz=nz, c_eff=c_eff, H=H,
                thr=thr, ratio=ratio)


def main():
    carr, car_vals, _, _, _, _ = rf.load()
    print('carriage 观测值 n=%d\n' % len(car_vals))

    # --- SG90 元数据 ---
    met = {}
    with open(os.path.join(BASE, 'SG90_metadata.csv'), encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            met[r['Subject ID']] = fnum(r['Age (in years)'])

    e_info, e_n, e_ids = load_matrix(os.path.join(BASE, 'ELDERMET_species_189x836.csv'))
    s_info, s_n, s_ids = load_matrix(os.path.join(BASE, 'SG90_species_abundances.csv'))
    o_ids = {k for k, a in met.items() if a >= 77}
    o_info, o_n, _ = load_matrix(os.path.join(BASE, 'SG90_species_abundances.csv'), o_ids)
    w_info, w_n, w_ids = load_matrix(os.path.join(BASE, 'Wave3_species_59x.csv'))

    cohorts = [('ELDERMET (IE)', e_info, e_n),
               ('SG90 all (SG)', s_info, s_n),
               ('SG90 77+ (SG)', o_info, o_n),
               ('Wave3 65+ (CN)', w_info, w_n)]

    res = {}
    print('=== 各队列全样本 H* ===')
    for tag, info, n in cohorts:
        res[tag] = summarise(tag, info, n, car_vals)

    print('\n=== hdcA 家族状态 ===')
    states = []
    for tag, info, n in cohorts:
        st = hdca_state(tag, info, carr, res[tag])
        if st:
            states.append(st)

    print('\n=== 同一 n_sp 下的跨队列比值（相对 ELDERMET） ===')
    base = res['ELDERMET (IE)']
    for tag, _, _ in cohorts[1:]:
        for n_sp in NS_LIST:
            if n_sp in res[tag]['ns'] and n_sp in base['ns']:
                print('  n_sp=%2d  ELDERMET=%.5f  %-14s=%.5f  ratio=%.2f'
                      % (n_sp, base['ns'][n_sp], tag, res[tag]['ns'][n_sp],
                         res[tag]['ns'][n_sp] / base['ns'][n_sp]))

    # --- 样本量匹配 ---
    print('\n=== 样本量匹配（三个队列都抽到 n=59，%d 次）===' % N_DRAW)
    nm = min(e_n, s_n, w_n)
    match_rows = []
    for tag, info, n in (('ELDERMET (IE)', e_info, e_n),
                         ('SG90 all (SG)', s_info, s_n),
                         ('SG90 77+ (SG)', o_info, o_n),
                         ('Wave3 65+ (CN)', w_info, w_n)):
        if n == nm:
            sub = info
            draws = [summarise(tag + ' [全样本]', sub, nm, car_vals, quiet=True)]
        else:
            rnd = random.Random(SEEDS[0])
            draws = []
            for d in range(N_DRAW):
                idx = rnd.sample(range(n), nm)
                draws.append(summarise('%s [抽#%d]' % (tag, d + 1),
                                       subset(info, idx), nm, car_vals, quiet=True))
        for n_sp in NS_LIST:
            vals = [x['ns'][n_sp] for x in draws if n_sp in x['ns']]
            if not vals:
                continue
            vals.sort()
            med = vals[len(vals) // 2]
            print('  %-16s n_sp=%2d  H* 中位=%.5f  [%.5f, %.5f]  %d 次'
                  % (tag, n_sp, med, vals[0], vals[-1], len(vals)))
            match_rows.append((tag, n_sp, med, vals[0], vals[-1], len(vals)))

    # --- 种子稳健性（Wave3） ---
    print('\n=== Wave3 种子稳健性 ===')
    wpool = [x for x in w_info if x['prev'] >= 0.20 and x['mean'] > 0]
    for s in SEEDS:
        fam5 = rf.families(car_vals, wpool, 5, w_n, seed=s)
        fam20 = rf.families(car_vals, wpool, 20, w_n, seed=s)
        print('  seed=%-10d  H*(5)=%.5f  H*(20)=%.5f'
              % (s, rf.h_cond(fam5), rf.h_cond(fam20)))

    # --- 输出 ---
    with open(os.path.join(BASE, 'cross_cohort3.tsv'), 'w', encoding='utf-8') as f:
        f.write('# cohort\tn_ind\tpool\tn_sp\tH_star_cond\tH_star_cum\n')
        for tag, _, _ in cohorts:
            r = res[tag]
            for n_sp in NS_LIST:
                if n_sp in r['ns']:
                    f.write('%s\t%d\t%d\t%d\t%.5f\t%.5f\n'
                            % (tag, r['n_ind'], r['pool'], n_sp,
                               r['ns'][n_sp], r['cum'][n_sp]))
        f.write('# matched-n (n=%d)\n# cohort\tn_sp\tmedian\tmin\tmax\tdraws\n' % nm)
        for row in match_rows:
            f.write('%s\t%d\t%.5f\t%.5f\t%.5f\t%d\n' % row)
        f.write('# hdcA\n# cohort\tmatched\tnonzero\tc_eff\tH\tHstar\tH/Hstar\n')
        for st in states:
            f.write('%s\t%d\t%d\t%.3f\t%.4f\t%.5f\t%.2f\n'
                    % (st['tag'], st['matched'], st['nz'], st['c_eff'],
                       st['H'], st['thr'], st['ratio']))
    print('\nwritten cross_cohort3.tsv')


if __name__ == '__main__':
    main()
