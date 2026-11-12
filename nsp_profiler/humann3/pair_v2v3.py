# -*- coding: utf-8 -*-
"""同一样本的 HUMAnN2 / HUMAnN3 配对对比：升级是否让族内异质性变小。

背景
----
HMP2 同时提供同一样本的两种处理：
  2018-05-04/func_profiles/<SID>_humann2.tar.bz2   -> <SID>_genefamilies.tsv（RPK）
  func_profile3_WGS/<SID>_humann3.tar.bz2          -> <SID>_genefamilies_cpm.tsv（CPM）
两者都带 per-species 分层行，所以能在同一个体上比较 H。

口径说明（与 analyze_humann3.py 有一处不同，是为了配对可比）
--------------------------------------------------------------
绝对阈值（如 0.5 CPM）不能跨版本用：v2 是 RPK、v3 是 CPM，量纲不同，同一个
阈值会把两版的携带者集合筛得不一样，差异全来自阈值而非流程。本脚本改用
**相对阈值**：一个物种的贡献须 ≥ 该家族总贡献的 --rel（默认 0.1%），两版同用，
与单位无关。H 本身经族内 max 归一，也单位无关，故可直接比。

用法：python pair_v2v3.py [--rel 0.001] [--min-p 0.1] [--min-sp 3]
"""
import argparse
import glob
import math
import os
import statistics
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
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


def read_fam(path):
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


def stats(contrib, ptab, rel, min_p, min_sp):
    tot_all = sum(contrib.values())
    if tot_all <= 0:
        return None
    carriers = {sp: v for sp, v in contrib.items()
                if ptab.get(sp, 0) >= min_p and v / tot_all >= rel}
    if len(carriers) < min_sp:
        return None
    tot = sum(carriers.values())
    tp = sum(ptab[sp] for sp in carriers)
    c = {sp: v / ptab[sp] for sp, v in carriers.items()}
    cmax = max(c.values())
    if cmax <= 0:
        return None
    r = {s: c[s] / cmax for s in c}
    w = {s: ptab[s] / tp for s in carriers}
    reff = sum(w[s] * r[s] for s in r)
    H = sum(w[s] * (r[s] - reff) ** 2 for s in r)
    return len(r), reff, H


def load_species():
    import gzip
    import csv
    p = os.path.join(BASE, 'hmp2_h3', 'taxonomic_profiles_3.tsv.gz')
    out = {}
    with gzip.open(p, 'rt', encoding='utf-8', errors='replace') as f:
        r = csv.reader(f, delimiter='\t')
        hdr = next(r)
        names = [h.replace('_profile', '') for h in hdr[1:]]
        rows = [x for x in r if '|s__' in x[0]]
    for i, nm in enumerate(names):
        d = {}
        for row in rows:
            sp = row[0].split('|s__')[-1]
            try:
                v = float(row[i + 1])
            except (ValueError, IndexError):
                continue
            if v > 0:
                d[sp] = v
        out[nm] = d
    return out


def spearman(a, b):
    n = len(a)
    if n < 3:
        return float('nan')
    ra = {v: i for i, v in enumerate(sorted(a))}
    rb = {v: i for i, v in enumerate(sorted(b))}
    A = [ra[v] for v in a]
    B = [rb[v] for v in b]
    ma, mb = sum(A) / n, sum(B) / n
    num = sum((A[i] - ma) * (B[i] - mb) for i in range(n))
    den = math.sqrt(sum((A[i] - ma) ** 2 for i in range(n)) *
                    sum((B[i] - mb) ** 2 for i in range(n)))
    return num / den if den else float('nan')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rel', type=float, default=0.001)
    ap.add_argument('--min-p', type=float, default=0.1)
    ap.add_argument('--min-sp', type=int, default=3)
    a = ap.parse_args()

    sp = load_species()
    v3files = {os.path.basename(p).replace('_genefamilies_cpm.tsv', ''): p
               for p in glob.glob(os.path.join(WORK, '*_genefamilies_cpm.tsv'))}
    v2files = {os.path.basename(p).replace('_genefamilies.tsv', ''): p
               for p in glob.glob(os.path.join(WORK, '*_genefamilies.tsv'))
               if not p.endswith('_cpm.tsv')}
    # v3 的样本名带 _P / _TR 后缀，v2 不带（实测 HSMA33LH_TR vs HSMA33LH），
    # 直接取交集会漏配，故按去掉后缀的 base 名配。
    import re as _re
    def _base(s):
        return _re.sub(r'_(P|TR)$', '', s)
    v3bybase = {_base(k): k for k in v3files}
    pairs = [(v3bybase[b], b) for b in sorted(set(v3bybase) & set(v2files))]
    print('v3 %d，v2 %d，可配对 %d' % (len(v3files), len(v2files), len(pairs)), flush=True)
    if not pairs:
        raise SystemExit('无配对样本')

    rows = []
    for sid3, sid2 in pairs:
        sid = sid3
        ptab = sp.get(sid3) or sp.get(sid2)
        if not ptab:
            print('  %-14s 无物种表，跳过' % sid, flush=True)
            continue
        f3, f2 = read_fam(v3files[sid3]), read_fam(v2files[sid2])
        common = 0
        H3s, H2s, ratios = [], [], []
        for f in set(f3) & set(f2):
            s3 = stats(f3[f], ptab, a.rel, a.min_p, a.min_sp)
            s2 = stats(f2[f], ptab, a.rel, a.min_p, a.min_sp)
            if not s3 or not s2:
                continue
            common += 1
            H3s.append(s3[2])
            H2s.append(s2[2])
            r3, r2 = s3[2] / hstar_at(s3[0]), s2[2] / hstar_at(s2[0])
            ratios.append((r3, r2))
        if not common:
            continue
        dh = [b - a2 for a2, b in zip(H2s, H3s)]
        flips = sum(1 for r3, r2 in ratios if (r3 > 1) != (r2 > 1))
        rows.append([sid, len(f3), len(f2), common,
                     '%.4f' % statistics.median(H2s),
                     '%.4f' % statistics.median(H3s),
                     '%.4f' % statistics.median(dh),
                     '%.3f' % spearman(H2s, H3s),
                     flips, '%.1f' % (100 * flips / common)])
        print('  %-14s 共有 %-6d H2中位 %.4f  H3中位 %.4f  ρ=%.3f  判决翻转 %d(%.1f%%)'
              % (rows[-1][0], common, rows[-1][4] and float(rows[-1][4]),
                 float(rows[-1][5]), float(rows[-1][7]), flips,
                 100 * flips / common), flush=True)

    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, 'paired_v2v3.tsv')
    with open(p, 'w', encoding='utf-8', newline='') as f:
        f.write('sample\tn_fam_v3\tn_fam_v2\tn_common\tH_v2_median\tH_v3_median'
                '\tH_diff_median\tspearman_H\tn_verdict_flip\tpct_flip\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')
    print('\n写入 %s' % p)


if __name__ == '__main__':
    main()
