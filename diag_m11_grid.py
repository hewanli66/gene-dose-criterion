# -*- coding: utf-8 -*-
"""
diag_m11_grid.py —— M11 敏感性网格（保真门槛 × 保留率）按八种子重算。

背景：v3 的网格用单一种子生成的家族（review_final.py 里的 fams[5]），
中位数版全文主值是八种子中位数，网格若仍用单种子，会与
"every tabulated value is the median of the eight" 的约定冲突（第三轮审稿指出）。

本脚本对每个种子各生成一批 n_sp = 5 的家族，在 9 个 (target, keep) 组合下
各求一次 h_cond，报每格的中位数与全距，以及 hdcA（H = 0.1522）相对阈值的倍数。

用法：
    python diag_m11_grid.py

输出：
    diag_m11_grid.txt
    review_M11_grid_median.tsv

运行环境：Python 3.8+，仅标准库；耗时约 1 分钟。
"""
import csv
import os
import statistics
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import review_final as rf
from cross_cohort9_seed import SEEDS, load_all, pool_of

BASE = os.path.dirname(os.path.abspath(__file__))
TARGETS = (0.85, 0.90, 0.95)
KEEPS = (0.90, 0.95, 0.99)
N_SP = 5
COHORT = 'ELDERMET (IE, MPA3)'


def main():
    t0 = time.time()
    carr, car_vals, _, _, _, _ = rf.load()
    info, n = load_all()[COHORT]
    pool = pool_of(info)
    n_ind = len(info[0]['vals'])

    # hdcA 的 H：与 Table 4 / cross_cohort9_seed_t5.tsv 同源
    H = None
    with open(os.path.join(BASE, 'cross_cohort9_seed_t5.tsv'), encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['cohort'] == COHORT:
                H = float(r['H'])
    assert H is not None

    # 每个种子一批家族，9 个组合共用（与 review_final.py 的做法一致）
    per_seed = []
    for s in SEEDS:
        fams = rf.families(car_vals, pool, N_SP, n_ind, seed=s)
        per_seed.append({(t, k): rf.h_cond(fams, target=t, keep=k)
                         for t in TARGETS for k in KEEPS})

    lines = []
    lines.append('=== M11 敏感性网格（八种子，%s，n_sp = %d）===' % (COHORT, N_SP))
    lines.append('N_TRIAL=%d  min_n=%d  n_ind=%d  hdcA H=%.4f'
                 % (rf.N_TRIAL, 40, n_ind, H))
    lines.append('')
    lines.append('%-10s %-12s %-12s %-12s' % ('保真门槛', '保留 90%', '保留 95%', '保留 99%'))

    med = {}
    for t in TARGETS:
        cells = []
        for k in KEEPS:
            vals = [d[(t, k)] for d in per_seed]
            m = statistics.median(vals)
            med[(t, k)] = m
            cells.append('%.4f[%.4f-%.4f]' % (m, min(vals), max(vals)))
        lines.append('%-10.2f %-12s %-12s %-12s' % (t, *cells))

    lo = min(med.values())
    hi = max(med.values())
    lines.append('')
    lines.append('网格中位数的范围：%.4f–%.4f（%.1f 倍）' % (lo, hi, hi / lo))
    lines.append('hdcA 相对网格阈值的倍数：%.1f – %.1f（全部 %d 格均 > 1）'
                 % (H / hi, H / lo, len(med)))

    # 每一格是否在每个种子下都 > 1
    worst = min(H / d[(t, k)] for d in per_seed for t in TARGETS for k in KEEPS)
    lines.append('最不利的单格单种子值：%.1f 倍（%s）'
                 % (worst, '仍 > 1，判定不随约定翻转' if worst > 1 else '**< 1，需改写**'))

    lines.append('')
    lines.append('耗时 %.0f 秒' % (time.time() - t0))
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(BASE, 'diag_m11_grid.txt'), 'w', encoding='utf-8') as f:
        f.write(txt + '\n')
    with open(os.path.join(BASE, 'review_M11_grid_median.tsv'), 'w',
              encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['fidelity_threshold', 'keep_rate', 'H_star_median',
                    'H_star_min', 'H_star_max', 'hdcA_ratio'])
        for t in TARGETS:
            for k in KEEPS:
                vals = [d[(t, k)] for d in per_seed]
                w.writerow(['%.2f' % t, '%.2f' % k, '%.5f' % med[(t, k)],
                            '%.5f' % min(vals), '%.5f' % max(vals),
                            '%.1f' % (H / med[(t, k)])])


if __name__ == '__main__':
    main()
