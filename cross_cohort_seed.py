# -*- coding: utf-8 -*-
"""
cross_cohort_seed.py —— 跨队列 H* 的种子稳健性检验。

背景：cross_cohort.py 在 SG90 八旬子集上得到非单调的规模依赖
（H*(20) = 0.088 < H*(10) = 0.109），与 ELDERMET 的单调上升不一致。
必须先排除“模拟噪声”这一解释，才能把它写成结论。

做法：每个队列 × 每个 n_sp 用 4 个不同随机种子各跑一次，
报告 H* 的中位数与全距；若全距跨过规模间差值，则该规模依赖不成立。

运行环境：Python 3.8+ 仅标准库。约 2 分钟。
输出：数据/cross_cohort_seed.tsv + 屏幕表
"""
import csv
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import review_final as rf
from cross_cohort import load_matrix, pool_of, fnum

BASE = os.path.dirname(os.path.abspath(__file__))
SEEDS = [rf.SEED, rf.SEED + 101, rf.SEED + 2027, rf.SEED + 40009]
NS_LIST = [5, 10, 20]
AGE_CUT = 77


def hstar_over_seeds(pool, n_sp, n_ind):
    vals = []
    for s in SEEDS:
        fam = rf.families(rf.load()[1], pool, n_sp, n_ind, seed=s)
        vals.append(rf.h_cond(fam))
    return vals


def main():
    car_vals = rf.load()[1]
    met = {}
    with open(os.path.join(BASE, 'SG90_metadata.csv'), encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            met[r['Subject ID']] = fnum(r['Age (in years)'])
    old_ids = {k for k, a in met.items() if a >= AGE_CUT}

    e_info, e_n, _ = load_matrix(os.path.join(BASE, 'ELDERMET_species_189x836.csv'))
    s_info, s_n, _ = load_matrix(os.path.join(BASE, 'SG90_species_abundances.csv'))
    o_info, o_n, _ = load_matrix(os.path.join(BASE, 'SG90_species_abundances.csv'), old_ids)

    cohorts = (('ELDERMET', e_info, e_n),
               ('SG90_all', s_info, s_n),
               ('SG90_77plus', o_info, o_n))

    rows = []
    print('%-12s %5s  %-42s %8s %8s' % ('cohort', 'n_sp', 'H* over 4 seeds', 'median', 'range%'))
    for tag, info, n in cohorts:
        pool = pool_of(info)
        prev_med = None
        for n_sp in NS_LIST:
            vals = hstar_over_seeds(pool, n_sp, n)
            med = statistics.median(vals)
            lo, hi = min(vals), max(vals)
            rng = (hi - lo) / med * 100
            flag = ''
            if prev_med is not None and (hi - lo) > abs(med - prev_med):
                flag = '  <-- 种子波动大于规模间差值'
            print('%-12s %5d  %-42s %8.5f %7.1f%%%s'
                  % (tag, n_sp, ' '.join('%.4f' % v for v in vals), med, rng, flag))
            rows.append((tag, n, len(pool), n_sp, vals, med, rng))
            prev_med = med

    with open(os.path.join(BASE, 'cross_cohort_seed.tsv'), 'w', encoding='utf-8') as f:
        f.write('cohort\tn_ind\tpool\tn_sp\t' + '\t'.join('H_seed%d' % i for i in range(len(SEEDS)))
                + '\tmedian\trange_pct\n')
        for tag, n, npool, n_sp, vals, med, rng in rows:
            f.write('%s\t%d\t%d\t%d\t%s\t%.5f\t%.1f\n'
                    % (tag, n, npool, n_sp, '\t'.join('%.5f' % v for v in vals), med, rng))
    print('\nwritten cross_cohort_seed.tsv')


if __name__ == '__main__':
    main()
