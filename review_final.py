# -*- coding: utf-8 -*-
"""Final numbers: every quantity cited in the text, computed in one place.

本脚本在修复 H* 的定义错误之后，重算正文需要引用的每一个数：

  §1  H* 表（条件概率口径，主）与累积口径（对照）
  §2  M11 的敏感性网格（保真门槛 × 保留率）
  §3  M3 逐样本检出
  §4  M4 携带率不确定性的传播
  §5  M5 的 P 值区间
  §6  M6 的基因长度归一化
  §7  Table 4 的最终比值

定义说明（为什么改口径）
------------------------
旧实现把 H* 取成 "H 的 30% 分位数"：搜索从 q=30 起步，而该处达标率恒为 1，
循环立刻返回，因此那张表与 "≥95% 达标" 几乎无关，却产生了随 n_sp 单调增、
跨两个数量级的假象。

正文的用途是判定**单个**基因族能否用 proxy，因此正确的量是
    P(ρ ≥ 0.90 | H ≈ h) ≥ 0.95
的最大 h（条件口径），而不是"用 h 筛掉一批族后整体达标"（累积口径）。
两者在 n_sp 依赖上差别很大：条件口径 2.7 倍，累积口径 1.1 倍。
本脚本主报条件口径，同时给出累积口径以便读者核对。
"""
import csv
import math
import os
import random

BASE = os.path.dirname(os.path.abspath(__file__))


def data_path(fn):
    for d in (BASE, os.path.join(os.path.dirname(BASE), '数据')):
        p = os.path.join(d, fn)
        if os.path.exists(p):
            return p
    return os.path.join(BASE, fn)
N_TRIAL = 1200
N_BOOT = 4000
SEED = 20261028
WIN = 0.02
TARGET = 0.90
KEEP = 0.95

GENOME_MB = {
    'Eggerthella_lenta': 3.2, 'Streptococcus_thermophilus': 1.9,
    'Lactobacillus_fermentum': 2.1, 'Bacteroides_nordii': 5.5,
    'Clostridium_perfringens': 3.256683,   # CP000246, 3,256,683 bp
}
GENOME_DEFAULT = 3.5
N_READS = 6_300_000
GENE_LEN = 1000
ABUND_UNIT = 100.0
EPS = 0.71
LAM_THR = 5.0


def fnum(s, d=0.0):
    try:
        return float(str(s).strip())
    except (TypeError, ValueError):
        return d


def rankdata(x):
    n = len(x)
    idx = sorted(range(n), key=lambda i: x[i])
    r = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and x[idx[j + 1]] == x[idx[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[idx[k]] = avg
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    n = len(rx)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else 0.0


def load():
    carr, car_vals = {}, []
    with open(data_path('TableS1_species_fixed.csv'), encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            sp = (r.get('species') or '').strip().replace(' ', '_')
            a, b = fnum(r.get('n_HSB_strains'), -1), fnum(r.get('n_total_strains'), -1)
            if sp and b > 0:
                carr[sp] = dict(k=int(a), n=int(b), frac=a / b,
                                typ=(r.get('hdcA_type') or '').strip())
                if a / b > 0:
                    car_vals.append(a / b)
    with open(data_path('ELDERMET_species_189x836.csv'), encoding='utf-8-sig') as f:
        rd = csv.reader(f)
        hdr = next(rd)
        cols = hdr[1:]
        rows = [r for r in rd]
    mat, info = {}, []
    n_ind = len(rows)
    for j, c in enumerate(cols):
        v = [fnum(r[j + 1]) for r in rows]
        mat[c] = v
        info.append(dict(mean=sum(v) / n_ind,
                         prev=sum(1 for x in v if x > 0) / n_ind, vals=v))
    return carr, car_vals, mat, info, cols, n_ind


def families(car_vals, info, n_sp, n_ind, seed=SEED):
    cand = [x for x in info if x['prev'] >= 0.20 and x['mean'] > 0]
    rnd = random.Random(seed)
    fam = []
    for _ in range(N_TRIAL):
        picks = rnd.sample(cand, n_sp)
        cs = rnd.choices(car_vals, k=n_sp)
        px, co = [], []
        for i in range(n_ind):
            p = s = 0.0
            for k, sp in enumerate(picks):
                a = sp['vals'][i]
                p += a
                s += a * cs[k]
            px.append(p)
            co.append(s)
        if max(co) <= 0 or max(px) <= 0:
            continue
        rho = spearman(px, co)
        tot = sum(sp['mean'] for sp in picks) or 1.0
        ws = [sp['mean'] / tot for sp in picks]
        c_eff = sum(w * c for w, c in zip(ws, cs))
        H = sum(w * (c - c_eff) ** 2 for w, c in zip(ws, cs))
        fam.append((H, rho))
    fam.sort()
    return fam


def h_cond(fam, target=TARGET, keep=KEEP, width=WIN, min_n=40):
    best = fam[0][0]
    hs = [h for h, _ in fam]
    i = 0
    for j in range(len(fam)):
        while hs[j] - hs[i] > width:
            i += 1
        if j - i + 1 < min_n:
            continue
        ok = sum(1 for k in range(i, j + 1) if fam[k][1] >= target) / (j - i + 1)
        if ok >= keep:
            best = hs[j]
        else:
            break
    return best


def h_cum(fam, target=TARGET, keep=KEEP):
    best = fam[0][0]
    for i in range(2, len(fam) + 1):
        sub = fam[:i]
        if sum(1 for _, r in sub if r >= target) / len(sub) >= keep:
            best = fam[i - 1][0]
        else:
            break
    return best


def q(xs, p):
    s = sorted(xs)
    i = min(len(s) - 1, max(0, int(round(p * (len(s) - 1)))))
    return s[i]


def pois_le(k, lam):
    return sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(k + 1))


if __name__ == '__main__':
    carr, car_vals, mat, info, cols, n_ind = load()
    L = []
    A = L.append
    A('=' * 76)
    A('Final numbers')
    A('=' * 76)
    A('')

    # ---- §1 H*
    fams = {}
    hstar = {}
    A('§1  H* 表（条件口径为主，累积口径对照；每档 %d 族）' % N_TRIAL)
    A('    %-6s %10s %10s %10s' % ('n_sp', 'H*(条件)', 'H*(累积)', '条件比值'))
    for n_sp in (3, 5, 10, 20):
        fam = families(car_vals, info, n_sp, n_ind)
        fams[n_sp] = fam
        hc, hu = h_cond(fam), h_cum(fam)
        hstar[n_sp] = hc
        A('    %-6d %10.5f %10.5f' % (n_sp, hc, hu))
    A('    条件口径范围 %.5f–%.5f，跨 %.1f 倍'
      % (min(hstar.values()), max(hstar.values()),
         max(hstar.values()) / min(hstar.values())))
    A('')

    # ---- §2 M11 网格（条件口径）
    A('§2  M11 敏感性网格（条件口径，n_sp=5，同一批家族）')
    A('    %-8s %10s %10s %10s' % ('保真门槛', '保留90%', '保留95%', '保留99%'))
    grid = {}
    for target in (0.85, 0.90, 0.95):
        row = []
        for keep in (0.90, 0.95, 0.99):
            h = h_cond(fams[5], target=target, keep=keep)
            grid[(target, keep)] = h
            row.append(h)
        A('    %-8.2f %10.5f %10.5f %10.5f' % (target, row[0], row[1], row[2]))
    A('    网格范围 %.5f–%.5f；hdcA(H=0.152) 比值 %.1f–%.1f 倍'
      % (min(grid.values()), max(grid.values()),
         0.152 / max(grid.values()), 0.152 / min(grid.values())))
    A('')
    with open(os.path.join(BASE, 'review_M11_grid.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['fidelity_threshold', 'keep_rate', 'H_star', 'hdcA_ratio'])
        for (t, k), h in sorted(grid.items()):
            w.writerow(['%.2f' % t, '%.2f' % k, '%.5f' % h, '%.1f' % (0.152 / h)])

    # ---- §3 M3
    members = []
    for s, v in carr.items():
        if s in mat and v['typ'].startswith('pyruvoyl'):
            a = sum(mat[s]) / n_ind
            if a > 0:
                members.append((s, a, v['frac'], GENOME_MB.get(s, GENOME_DEFAULT)))
    members.sort(key=lambda z: -z[1])

    # 家族 λ_raw（队列均值口径），供 §5 的 ε 敏感性用，不写死
    TOT_RAW = sum(N_READS * (a / ABUND_UNIT) * c * GENE_LEN / (g * 1e6)
                  for (s, a, c, g) in members)

    lam_det_i = []
    for i in range(n_ind):
        lam = 0.0
        for (s, a_m, c, g) in members:
            lam += N_READS * (mat[s][i] / ABUND_UNIT) * c * GENE_LEN / (g * 1e6)
        lam_det_i.append(lam * EPS)
    p_i = [1 - math.exp(-x) for x in lam_det_i]
    p_mean = sum(p_i) / len(p_i)
    A('§3  M3 逐样本检出（hdcA 家族 %d 个物种，189 人）' % len(members))
    A('    平均 λ_det = %.3f（与队列均值口径一致，因线性）' % (sum(lam_det_i) / len(lam_det_i)))
    A('    个体间 λ_det 分位：10%%=%.3f  50%%=%.3f  90%%=%.3f'
      % (q(lam_det_i, .1), q(lam_det_i, .5), q(lam_det_i, .9)))
    A('    P(≥1 read) 均值 = %.3f；189 人预期阳性 = %.0f 人（%.1f%%）'
      % (p_mean, sum(p_i), 100 * sum(p_i) / n_ind))
    A('    完全测不到（λ_det<0.01）的个体占 %.1f%%'
      % (100 * sum(1 for x in lam_det_i if x < 0.01) / n_ind))
    A('    任取 5 样本：预期阳性 %.2f 个；至少 1 个阳性概率 %.2f'
      % (5 * p_mean, 1 - (1 - p_mean) ** 5))
    A('')

    # ---- §4 M4
    rnd = random.Random(SEED)
    tot_ab = sum(m[1] for m in members)
    ws = [m[1] / tot_ab for m in members]
    Hs, ratios, lams = [], [], []
    both = fid = det = 0
    for _ in range(N_BOOT):
        cs = [rnd.betavariate(carr[m[0]]['k'] + .5, carr[m[0]]['n'] - carr[m[0]]['k'] + .5)
              for m in members]
        ce = sum(w * c for w, c in zip(ws, cs))
        H = sum(w * (c - ce) ** 2 for w, c in zip(ws, cs))
        lam = sum(N_READS * (m[1] / ABUND_UNIT) * c * GENE_LEN / (m[3] * 1e6)
                  for m, c in zip(members, cs)) * EPS
        Hs.append(H); ratios.append(H / hstar[5]); lams.append(lam)
        f_fail = H / hstar[5] > 1
        d_fail = lam < LAM_THR
        if f_fail and d_fail:
            both += 1
        elif f_fail:
            fid += 1
        elif d_fail:
            det += 1
    A('§4  M4 携带率不确定性（Beta 后验，%d 次）' % N_BOOT)
    A('    H 中位 %.3f（90%% 区间 %.3f–%.3f）'
      % (q(Hs, .5), q(Hs, .05), q(Hs, .95)))
    A('    H/H*(n=5) 中位 %.2f（区间 %.2f–%.2f）'
      % (q(ratios, .5), q(ratios, .05), q(ratios, .95)))
    A('    λ_det 中位 %.3f（区间 %.3f–%.3f）'
      % (q(lams, .5), q(lams, .05), q(lams, .95)))
    A('    P(判定 both fail) = %.4f' % (both / N_BOOT))
    A('')

    # ---- §5 M5
    A('§5  M5 P(X≤1) 随 ε 变化（5 样本期望 = %.3f·ε·5）' % TOT_RAW)
    A('    %-8s %12s %12s' % ('ε', '期望', 'P(X≤1)'))
    for eps, tag in ((0.48, '下限'), (0.54, '逐样本下限'), (0.71, '点估计'), (0.94, '逐样本上限')):
        lam5 = TOT_RAW * eps * 5
        A('    %-8.2f %12.2f %12.4f  %s' % (eps, lam5, pois_le(1, lam5), tag))
    A('')

    # ---- §6 M6
    A('§6  M6 基因长度归一化')
    try:
        import json
        ref = json.load(open(data_path('rpob_ref.json'), encoding='utf-8'))
        aa = [v['length'] for v in ref.values()]
        kb = 3 * (sum(aa) / len(aa)) / 1000
        A('    rpoB 均值 %.0f aa → %.2f kb；hdcA 1.00 kb；长度因子 %.2f' % (sum(aa) / len(aa), kb, kb))
        A('    原始 19,691:1 → 长度归一化 %.0f:1' % (19691 / kb))
    except Exception as e:
        A('    （失败：%s）' % e)
    A('')

    # ---- §7 Table 4
    A('§7  Table 4 的 H/H*（条件口径）')
    for label, n_sp, H in (('hdcA all 20', 20, 0.152), ('hdcA contributing 5', 5, 0.152),
                           ('hdcA carriage≥40%', 3, 0.081)):
        A('    %-22s n=%2d  H=%.3f  H*=%.5f  ratio=%.2f'
          % (label, n_sp, H, hstar[n_sp], H / hstar[n_sp]))
    A('')

    txt = '\n'.join(L)
    open(os.path.join(BASE, 'review_response.txt'), 'w', encoding='utf-8').write(txt + '\n')
    print(txt)
