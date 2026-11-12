# -*- coding: utf-8 -*-
"""用 HMP2 的 HUMAnN3 实测输出，测「物种代理」的误差有多大。

数据来源（2026-11-11 实测可下载）
----------------------------------
物种表  HMP2/MGX/2018-05-04/taxonomic_profiles_3.tsv.gz   MetaPhlAn 3，相对丰度 %
功能表  HMP2/MGX/func_profile3_WGS/<SID>_humann3.tar.bz2  内 *genefamilies_cpm.tsv
        含分层行 `UniRef90_XXXX|g__Genus.s__Species` = per-species 贡献（CPM）

为什么必须先归一化（2026-11-11 踩坑）
-------------------------------------
HUMAnN 的 contrib 是 CPM、物种表是百分数，c_s = contrib/p_s 的量纲是 CPM/%，
不是正文那种 [0,1] 的株水平携带率。直接拿它和 §5.4 的 H*（0.024–0.094）比，
H/H* 会虚高到 10³–10⁶，毫无意义。故对每个家族在族内归一：
    c_s      = contrib_{f,s} / p_s
    r_s      = c_s / max_{s}(c_s)          ∈ (0,1]，最大者记为"典型携带"
    w_s      = p_s / Σ p_s                 （携带物种内）
    r_eff    = Σ w_s r_s
    H        = Σ w_s (r_s − r_eff)²        ← 与 H* 同量纲，可比
代理误差（尺度无关，故不必归一化）：
    err_s    = |r_eff − r_s| / r_s         （即 |c_eff−c_s|/c_s，p_s 约掉）
    err_f    = Σ_s 贡献_s · err_s / Σ_s 贡献_s     按贡献加权
另给 CV = sd(r)/mean(r)，作为不受 n 影响的离散度（min–max 会随物种数变宽）。

输出
----
out/humann3_by_sample.tsv   逐样本逐家族
out/humann3_families.tsv    家族级（跨样本取中位）
out/humann3_summary.txt     汇总
"""
import argparse
import csv
import glob
import gzip
import math
import os
import statistics
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
SPTAB = os.path.join(BASE, 'hmp2_h3', 'taxonomic_profiles_3.tsv.gz')
WORK = os.path.join(BASE, 'work')
OUT = os.path.join(BASE, 'out')

HSTAR = {3: 0.0240, 5: 0.0378, 10: 0.0732, 20: 0.0936}


def hstar_at(n):
    ks = sorted(HSTAR)
    if n <= ks[0]:
        return HSTAR[ks[0]]
    if n >= ks[-1]:
        return HSTAR[ks[-1]]
    for a, b in zip(ks, ks[1:]):
        if a <= n <= b:
            t = (math.log(n) - math.log(a)) / (math.log(b) - math.log(a))
            return math.exp(math.log(HSTAR[a]) * (1 - t) + math.log(HSTAR[b]) * t)
    return HSTAR[ks[-1]]


def read_species(sptab):
    out = {}
    with gzip.open(sptab, 'rt', encoding='utf-8', errors='replace') as f:
        r = csv.reader(f, delimiter='\t')
        hdr = next(r)
        samples = [h.replace('_profile', '') for h in hdr[1:]]
        for row in r:
            feat = row[0]
            if '|s__' not in feat:
                continue
            sp = feat.split('|s__')[-1].strip()
            if not sp or sp.endswith('_unclassified'):
                continue
            for j, s in enumerate(samples):
                if j + 1 >= len(row):
                    continue
                try:
                    v = float(row[j + 1])
                except ValueError:
                    continue
                if v > 0:
                    out.setdefault(s, {})[sp] = v
    return out


def read_genefamilies(path):
    fam = {}
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 2 or '|' not in p[0]:
                continue
            base, tax = p[0].split('|', 1)
            if not tax.startswith('g__') or '.s__' not in tax:
                continue
            sp = tax.split('.s__')[-1]
            if not sp or sp == 'unclassified':
                continue
            try:
                v = float(p[1])
            except ValueError:
                continue
            if v > 0:
                fam.setdefault(base, {})[sp] = v
    return fam


def family_stats(contrib, ptab, min_p, min_sp, min_cpm, min_contrib=0.5):
    """单个样本、单个家族 -> (n_sp, r_eff, H, err, cv) 或 None

    三处必要的稳定性处理（2026-11-11 实测踩坑）：
      · c_s = contrib/p 在 p 极小时不稳定，会抬高 max(c_s) 使归一化失真
        → 只保留 contrib ≥ min_contrib CPM 的物种参与；
      · r_eff 极小时 /r_eff 形式的误差会爆炸（见过 6581）
        → 误差改用绝对差，有界于 [0,1]；
      · 携带物种太少或家族总量太小 → 直接跳过。
    """
    carriers = {sp: v for sp, v in contrib.items()
                if ptab.get(sp, 0) >= min_p and v >= min_contrib}
    if len(carriers) < min_sp:
        return None
    tot = sum(carriers.values())
    if tot < min_cpm:
        return None
    tot_p = sum(ptab[sp] for sp in carriers)
    c = {sp: v / ptab[sp] for sp, v in carriers.items()}
    cmax = max(c.values())
    if cmax <= 0:
        return None
    r = {sp: c[sp] / cmax for sp in c}
    w = {sp: ptab[sp] / tot_p for sp in carriers}
    reff = sum(w[s] * r[s] for s in r)
    H = sum(w[s] * (r[s] - reff) ** 2 for s in r)
    # 绝对差（以"典型携带强度"=1 为单位），按贡献加权
    err = sum(carriers[s] * abs(reff - r[s]) for s in r) / tot
    mean_r = reff
    sd = math.sqrt(H)
    cv = sd / mean_r if mean_r else float('nan')
    return len(r), reff, H, err, cv
    mean_r = reff
    sd = math.sqrt(H)
    cv = sd / mean_r if mean_r else float('nan')
    return len(r), reff, H, err, cv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-p', type=float, default=0.01)
    ap.add_argument('--min-sp', type=int, default=3)
    ap.add_argument('--min-cpm', type=float, default=1.0)
    ap.add_argument('--exclude-shallow', type=int, default=30000,
                    help='分层行少于该数的样本视为浅测序，跳过')
    a = ap.parse_args()

    species = read_species(SPTAB)
    files = sorted(glob.glob(os.path.join(WORK, '*_genefamilies_cpm.tsv')))
    print('物种表 %d 样本；功能表 %d 样本' % (len(species), len(files)), flush=True)

    per_sample = []
    skipped = []
    for path in files:
        sid = os.path.basename(path).replace('_genefamilies_cpm.tsv', '')
        ptab = species.get(sid)
        if not ptab:
            skipped.append((sid, 'no_species_table'))
            continue
        fam = read_genefamilies(path)
        nstrat = sum(len(d) for d in fam.values())
        if nstrat < a.exclude_shallow:
            skipped.append((sid, 'shallow(%d)' % nstrat))
            continue
        n = 0
        for f, d in fam.items():
            st = family_stats(d, ptab, a.min_p, a.min_sp, a.min_cpm)
            if not st:
                continue
            nsp, reff, H, err, cv = st
            hs = hstar_at(nsp)
            per_sample.append([sid, f, nsp, '%.4f' % reff, '%.5f' % H,
                               '%.5f' % hs, '%.3f' % (H / hs), '%.4f' % err,
                               '%.3f' % cv])
            n += 1
        print('  %-14s 家族 %-6d 合格 %d' % (sid, len(fam), n), flush=True)
    if not per_sample:
        sys.exit('无合格数据')

    os.makedirs(OUT, exist_ok=True)
    p1 = os.path.join(OUT, 'humann3_by_sample.tsv')
    with open(p1, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['sample', 'family', 'n_sp', 'r_eff', 'H', 'H_star', 'H_over_Hstar',
                    'proxy_err', 'cv_r'])
        w.writerows(per_sample)

    # 家族级：跨样本取中位
    byfam = {}
    for row in per_sample:
        byfam.setdefault(row[1], []).append(row)
    famrows = []
    for f, rows in byfam.items():
        if len(rows) < 1:
            continue
        med = lambda i: statistics.median(float(r[i]) for r in rows)
        famrows.append([f, len(rows), int(med(2)), '%.4f' % med(3), '%.5f' % med(4),
                        '%.5f' % med(5), '%.3f' % med(6), '%.4f' % med(7),
                        '%.3f' % med(8)])
    famrows.sort(key=lambda r: -float(r[4]))
    p2 = os.path.join(OUT, 'humann3_families.tsv')
    with open(p2, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['family', 'n_samples', 'n_sp', 'r_eff', 'H', 'H_star',
                    'H_over_Hstar', 'proxy_err', 'cv_r'])
        w.writerows(famrows)

    def pct(v, q):
        s = sorted(v)
        return s[min(len(s) - 1, int(q * len(s)))]

    H = [float(r[4]) for r in famrows]
    E = [float(r[7]) for r in famrows]
    R = [float(r[6]) for r in famrows]
    C = [float(r[8]) for r in famrows]
    lines = []
    lines.append('样本数 %d（跳过 %d）' % (len(files) - len(skipped), len(skipped)))
    for s, why in skipped:
        lines.append('  skip %-14s %s' % (s, why))
    lines.append('家族数 %d，逐样本记录 %d' % (len(famrows), len(per_sample)))
    lines.append('H       中位 %.5f  p90 %.5f  最大 %.5f' %
                 (statistics.median(H), pct(H, .9), max(H)))
    lines.append('proxy_err 中位 %.3f  p90 %.3f  最大 %.3f' %
                 (statistics.median(E), pct(E, .9), max(E)))
    lines.append('H/H*    中位 %.3f  p90 %.3f  最大 %.3f' %
                 (statistics.median(R), pct(R, .9), max(R)))
    lines.append('CV(r)   中位 %.3f  p90 %.3f' % (statistics.median(C), pct(C, .9)))
    lines.append('H/H*>1 的家族占比 %.1f%%' % (100 * sum(1 for x in R if x > 1) / len(R)))
    lines.append('proxy_err>0.5 的家族占比 %.1f%%' %
                 (100 * sum(1 for x in E if x > 0.5) / len(E)))

    # 秩相关
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        rk = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                rk[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return rk

    def spear(x, y):
        rx, ry = rank(x), rank(y)
        mx, my = statistics.mean(rx), statistics.mean(ry)
        num = sum((p - mx) * (q - my) for p, q in zip(rx, ry))
        den = math.sqrt(sum((p - mx) ** 2 for p in rx) * sum((q - my) ** 2 for q in ry))
        return num / den if den else float('nan')

    lines.append('Spearman(H, proxy_err) = %.3f' % spear(H, E))
    lines.append('Spearman(H/H*, proxy_err) = %.3f' % spear(R, E))
    # 判据的实际用途：H/H*>1 能否挑出 err>0.5 的家族
    tp = sum(1 for r, e in zip(R, E) if r > 1 and e > 0.5)
    fp = sum(1 for r, e in zip(R, E) if r > 1 and e <= 0.5)
    fn = sum(1 for r, e in zip(R, E) if r <= 1 and e > 0.5)
    tn = sum(1 for r, e in zip(R, E) if r <= 1 and e <= 0.5)
    lines.append('H/H*>1 vs err>0.5：TP=%d FP=%d FN=%d TN=%d' % (tp, fp, fn, tn))
    if tp + fp:
        lines.append('  精确率 %.3f' % (tp / (tp + fp)))
    if tp + fn:
        lines.append('  召回率 %.3f' % (tp / (tp + fn)))
    txt = '\n'.join(lines)
    with open(os.path.join(OUT, 'humann3_summary.txt'), 'w', encoding='utf-8') as f:
        f.write(txt + '\n')
    print('\n' + txt)
    print('\n写出\n  %s\n  %s' % (p1, p2))


if __name__ == '__main__':
    main()
