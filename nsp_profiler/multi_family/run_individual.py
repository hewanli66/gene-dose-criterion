# -*- coding: utf-8 -*-
"""个体层面的判决稳定性：H 用队列均值权重算，但失真发生在个体身上。

要回答的问题
------------
§3.3 已披露：「H is computed with cohort-mean weights while the distortion it
predicts is individual-level」。这条限制此前只有模拟能碰。HMP2 的 1,638 个
个体各有一份物种丰度表，于是可以对每个个体单独算 n_sp、c_eff、H、H*，看
同一个家族的判决是否随个体而变——这是实测，不是模拟。

与队列均值口径的差别
--------------------
队列均值：一个家族一个判决，权重是 1,638 人的平均丰度。
个体层面：一个家族 1,638 个判决，每人用自己的丰度；n_sp 也可能变（某人缺某物种）。

输出
----
out/individual_hmp2.tsv      家族级：判决为 pass / fail 的个体比例、H 的分位
out/individual_summary.txt   汇总

用法：python run_individual.py [--min-ab 0.001]
"""
import argparse
import csv
import gzip
import math
import os
import statistics
from collections import Counter

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))
HMP2 = os.path.join(ROOT, 'experiments-p2-p3', 'humann3', 'hmp2_h3',
                    'taxonomic_profiles_3.tsv.gz')
PANEL = os.path.join(BASE, 'out', 'carriage_panel.tsv')
OUT = os.path.join(BASE, 'out')

HSTAR = {3: 0.0240, 5: 0.0378, 10: 0.0732, 20: 0.0936}
READS = 6.3e6
EPS = 0.71


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


def load_hmp2(keep):
    """返回 sample -> {species: abundance%}，仅保留 keep 中的物种。"""
    per = {}
    with gzip.open(HMP2, 'rt', encoding='utf-8', errors='replace') as f:
        r = csv.reader(f, delimiter='\t')
        hdr = next(r)
        samples = hdr[1:]
        for row in r:
            if '|s__' not in row[0]:
                continue
            nm = row[0].split('|s__')[-1].replace('_', ' ')
            if nm not in keep:
                continue
            for i, v in enumerate(row[1:]):
                try:
                    x = float(v)
                except ValueError:
                    continue
                if x > 0:
                    per.setdefault(samples[i], {})[nm] = x
    return per


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-ab', type=float, default=0.001)
    a = ap.parse_args()

    byfam = {}
    allsp = set()
    with open(PANEL, encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            byfam.setdefault(r['pfam'], []).append(
                (r['species'], float(r['c_s']),
                 float(r['median_hit_len']) * 3, float(r['genome_bp'])))
            allsp.add(r['species'])
    per = load_hmp2(allsp)
    print('家族 %d，物种 %d，个体 %d' % (len(byfam), len(allsp), len(per)), flush=True)

    rows = []
    for pf, carriers in sorted(byfam.items()):
        npass = nfail = ndet = 0
        Hs = []
        for smp, ab in per.items():
            cs = []
            for sp, c, L, G in carriers:
                p = ab.get(sp, 0.0)
                if p > a.min_ab:
                    cs.append((p, c, L, G))
            if not cs:
                continue
            wsum = sum(x[0] for x in cs)
            ceff = sum(x[0] / wsum * x[1] for x in cs)
            H = sum(x[0] / wsum * (x[1] - ceff) ** 2 for x in cs)
            hs = hstar_at(len(cs))
            lam = sum(READS * x[0] / 100.0 * x[1] * x[2] / x[3] for x in cs if x[3])
            Hs.append(H)
            if H <= hs and lam * EPS >= 5:
                npass += 1
            elif H > hs and lam * EPS >= 5:
                nfail += 1
            else:
                ndet += 1
        tot = npass + nfail + ndet
        if not tot:
            continue
        rows.append([pf, tot, npass, nfail, ndet,
                     '%.3f' % (npass / tot),
                     '%.4f' % statistics.median(Hs),
                     '%.4f' % sorted(Hs)[int(0.95 * len(Hs))]])

    p = os.path.join(OUT, 'individual_hmp2.tsv')
    with open(p, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['pfam', 'n_ind', 'n_pass', 'n_fail_fid', 'n_fail_det',
                    'frac_pass', 'H_median', 'H_p95'])
        w.writerows(rows)

    # 汇总：判决一致性的分布
    cats = Counter()
    for r in rows:
        fp = float(r[5])
        if fp >= 0.99:
            cats['≥99% 个体一致'] += 1
        elif fp >= 0.95:
            cats['95–99%'] += 1
        elif fp >= 0.5:
            cats['50–95%'] += 1
        else:
            cats['<50%（多数个体给出相反判决）'] += 1
    # 队列均值判决 vs 个体多数判决
    fam = {}
    with open(os.path.join(OUT, 'families_panel.tsv'), encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            fam[r['pfam']] = r
    agree = 0
    n = 0
    for r in rows:
        e = fam.get(r[0])
        if not e:
            continue
        n += 1
        v_mean = float(e['H_over_Hstar']) <= 1 and float(e['lambda_det']) >= 5
        v_maj = float(r[5]) >= 0.5
        agree += (v_mean == v_maj)
    lines = ['个体层面判决稳定性（HMP2，%d 个体 × %d 家族）' % (len(per), len(rows)),
             '',
             '%-32s %6s' % ('判决一致度', '家族数')]
    for k in ['≥99% 个体一致', '95–99%', '50–95%', '<50%（多数个体给出相反判决）']:
        lines.append('%-32s %6d' % (k, cats[k]))
    lines += ['',
              '队列均值判决 与 个体多数判决 一致：%d/%d（%.1f%%）'
              % (agree, n, 100 * agree / max(1, n))]
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(OUT, 'individual_summary.txt'), 'w',
              encoding='utf-8') as f:
        f.write(txt + '\n')
    print('\n写入 %s' % p)


if __name__ == '__main__':
    main()
