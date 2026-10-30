# -*- coding: utf-8 -*-
"""补充：条件口径 H* 在三个队列子集内的稳定性。

2026-10-29 改：N_TRIAL / MIN_N 可由环境变量覆盖。
原值 N_TRIAL=900、min_n=30 与主表（cross_cohort9.py，N_TRIAL=1200、min_n=40）
不一致，导致正文那句"三子集跨度 1.4 倍"来自另一套集成。默认已改为主表口径，
输出同时落盘到 supplement_hstar.txt，便于复核。
用法：python supplement_hstar.py [输出文件名]
"""
import csv
import math
import os
import random
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
N_TRIAL = int(os.environ.get('N_TRIAL', '1200'))
MIN_N = int(os.environ.get('MIN_N', '40'))
WIN = 0.02
SEED0 = 20261028
SEEDS = [SEED0, SEED0 + 101, SEED0 + 2027, SEED0 + 40009,
         SEED0 + 7919, SEED0 + 104729, SEED0 + 1299709, SEED0 + 15485863]


def fnum(s):
    try:
        return float(str(s).strip())
    except (TypeError, ValueError):
        return 0.0


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


def h_cond(fam, target=0.90, keep=0.95, width=WIN, min_n=None):
    if min_n is None:
        min_n = MIN_N
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


def families(car_vals, info, n_sp, n_ind, seed=20261028):
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
        fam.append((sum(w * (c - c_eff) ** 2 for w, c in zip(ws, cs)), rho))
    fam.sort()
    return fam


def main():
    car_vals = []
    with open(os.path.join(BASE, 'TableS1_species_fixed.csv'), encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            a, b = fnum(r.get('n_HSB_strains')), fnum(r.get('n_total_strains'))
            if b > 0 and a / b > 0:
                car_vals.append(a / b)
    with open(os.path.join(BASE, 'ELDERMET_species_189x836.csv'), encoding='utf-8-sig') as f:
        rd = csv.reader(f)
        hdr = next(rd)
        cols = hdr[1:]
        rows = [r for r in rd]
    strat = {}
    with open(os.path.join(BASE, 'ELDERMET_clinical_189.csv'), encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            strat[r.get('SampleID', '')] = (r.get('Stratification') or '').strip()

    sids = [r[0] for r in rows]
    groups = {
        'all': list(range(len(rows))),
        'community': [i for i, s in enumerate(sids) if strat.get(s) == '1'],
        'care': [i for i, s in enumerate(sids) if strat.get(s) in ('3', '4')],
    }
    lines = []
    lines.append('=== 条件口径 H* 的三子集稳定性（八种子中位数）===')
    lines.append('N_TRIAL=%d  min_n=%d  width=%.2f' % (N_TRIAL, MIN_N, WIN))
    lines.append('seeds=%s' % (SEEDS,))
    lines.append('')
    lines.append('--- 单个种子（逐格原始值）---')
    out = {}
    for name, idx in groups.items():
        sub = [rows[i] for i in idx]
        info = []
        for j in range(len(cols)):
            v = [fnum(r[j + 1]) for r in sub]
            info.append(dict(mean=sum(v) / len(v),
                             prev=sum(1 for x in v if x > 0) / len(v), vals=v))
        vals = []
        for ns in (5, 10, 20):
            vv = [h_cond(families(car_vals, info, ns, len(sub), seed=s)) for s in SEEDS]
            vals.append(vv)
            lines.append('%-10s n=%3d n_sp=%2d  %s' %
                         (name, len(sub), ns,
                          '  '.join('%.5f' % x for x in vv)))
        out[name] = (len(sub), vals)
    lines.append('')
    lines.append('--- 每格八种子中位数 [极差] ---')
    lines.append('%-12s %5s %22s %22s %22s' % ('subset', 'n', 'H*(5)', 'H*(10)', 'H*(20)'))
    med = {}
    for name, (nsub, vals) in out.items():
        cells = []
        m = []
        for vv in vals:
            s = sorted(vv)
            md = (s[3] + s[4]) / 2.0
            m.append(md)
            cells.append('%.5f[%.0f%%]' % (md, 100.0 * (max(vv) - min(vv)) / md))
        med[name] = m
        lines.append('%-12s %5d %22s %22s %22s' % (name, nsub, cells[0], cells[1], cells[2]))
    lines.append('')
    lines.append('--- 跨子集跨度（按中位数算）---')
    for k, lab in ((0, 'H*(5)'), (1, 'H*(10)'), (2, 'H*(20)')):
        v = [m[k] for m in med.values()]
        lines.append('%s 跨子集：%.5f–%.5f，最大/最小 = %.2f 倍'
                     % (lab, min(v), max(v), max(v) / min(v)))
    txt = '\n'.join(lines)
    print(txt)
    out_name = sys.argv[1] if len(sys.argv) > 1 else 'supplement_hstar.txt'
    with open(os.path.join(BASE, out_name), 'w', encoding='utf-8') as f:
        f.write(txt + '\n')


if __name__ == '__main__':
    main()
