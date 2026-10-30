# -*- coding: utf-8 -*-
"""
diag_window_width.py —— n_sp = 3 端点"弱识别"的窗口宽度敏感性，按八种子重算。

背景：v3（单种子）报的是 0.01/0.02/0.03 三个窗口给 0.0118/0.0343/0.0318。
中位数版主值为八种子中位数（n_sp=3 处 0.0240），若仍引单种子的 0.034 与
0.012，同一句里会出现两套口径。本脚本对每个种子生成同一批家族，再在三个
窗口宽度下各求一次 h_cond，报八种子的中位数与全距。

用法：
    python diag_window_width.py

输出：
    diag_window_width.txt

运行环境：Python 3.8+，仅标准库；耗时约 1 分钟（ELDERMET 189 人，N_TRIAL=1200）。
"""
import os
import statistics
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import review_final as rf
from cross_cohort9_seed import SEEDS, load_all, pool_of

BASE = os.path.dirname(os.path.abspath(__file__))
WIDTHS = [0.01, 0.02, 0.03]
SIZES = [3, 5, 10, 20]
COHORT = 'ELDERMET (IE, MPA3)'


def main():
    t0 = time.time()
    carr, car_vals, _, _, _, _ = rf.load()
    info, n = load_all()[COHORT]
    pool = pool_of(info)
    n_ind = len(info[0]['vals'])

    lines = []
    lines.append('=== 窗口宽度敏感性（八种子，%s）===' % COHORT)
    lines.append('N_TRIAL=%d  min_n=%d  n_ind=%d  pool=%d'
                 % (rf.N_TRIAL, 40, n_ind, len(pool)))
    lines.append('')

    res = {}
    for n_sp in SIZES:
        # 每个种子生成一次家族，三个窗口共用同一批家族（与主表口径一致）
        fams = [rf.families(car_vals, pool, n_sp, n_ind, seed=s) for s in SEEDS]
        lines.append('n_sp = %d' % n_sp)
        for w in WIDTHS:
            vals = [rf.h_cond(f, width=w) for f in fams]
            med = statistics.median(vals)
            res[(n_sp, w)] = med
            lines.append('  width=%.2f  median=%.4f  min=%.4f  max=%.4f  raw=%s'
                         % (w, med, min(vals), max(vals),
                            ' '.join('%.4f' % v for v in vals)))
        lines.append('')

    for n_sp in SIZES:
        lo = res[(n_sp, 0.01)]
        mid = res[(n_sp, 0.02)]
        hi = res[(n_sp, 0.03)]
        lines.append('n_sp = %d：width 0.01 → %.4f，0.02 → %.4f，0.03 → %.4f'
                     % (n_sp, lo, mid, hi))
        lines.append('       0.01 相对 0.02 变化 %.0f%%，0.03 相对 0.02 变化 %.0f%%'
                     % (abs(lo - mid) / mid * 100, abs(hi - mid) / mid * 100))
        lines.append('       三窗口极差 / 中位 = %.2f 倍'
                     % (max(lo, mid, hi) / min(lo, mid, hi)))

    lines.append('')
    lines.append('耗时 %.0f 秒' % (time.time() - t0))
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(BASE, 'diag_window_width.txt'), 'w', encoding='utf-8') as f:
        f.write(txt + '\n')


if __name__ == '__main__':
    main()
