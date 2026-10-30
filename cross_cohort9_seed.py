# -*- coding: utf-8 -*-
"""
cross_cohort9_seed.py —— 九队列 + SG90 八旬子集，H* 的「多种子中位数」版。

用途：给"报中位数"那版稿件提供主值。与 cross_cohort9.py 的唯一区别是
每个 (队列 × n_sp) 不再只跑一个种子，而是跑 8 个种子取中位数
（前 4 个与旧表相同，便于逐格对照；后 4 个为新增，因为 n_sp=3 处四种子
中位数本身不稳）。除此之外，load_matrix / hdcA 匹配 / 对数插值全部沿用
cross_cohort9.py，保证两版稿件除"种子口径"外没有其他差异。

用法：
    python cross_cohort9_seed.py t5     # Table 5：10 行 × 四档 × 八种子（约 20 分钟）
    python cross_cohort9_seed.py t6     # Table 6：n=59 对齐，12 抽样 × 八种子（约 25 分钟）

输出：
    cross_cohort9_seed_t5.tsv
    cross_cohort9_seed_t6.tsv

运行环境：Python 3.8+，仅标准库。
"""
import os
import statistics
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import review_final as rf
from cross_cohort9 import load_matrix, match_hdca, interp, subset, NS_LIST, N_DRAW, fnum

BASE = os.path.dirname(os.path.abspath(__file__))

# 前四个与 cross_cohort_seed.py 完全相同，便于与旧表逐格对照；
# 后四个为本次新增（冒烟测试显示 n_sp=3 的四种子中位数本身不稳，中间两值差 8%，
# 故扩到 8 个种子让中位数站得住）。全部为质数间隔，彼此不构成等比。
SEEDS = [rf.SEED, rf.SEED + 101, rf.SEED + 2027, rf.SEED + 40009,
         rf.SEED + 7919, rf.SEED + 104729, rf.SEED + 1299709, rf.SEED + 15485863]

SPEC = [
    ('ELDERMET (IE, MPA3)',           'ELDERMET_species_189x836.csv',  None),
    ('CRC9 pool (MPA3)',              'CRC9_samples_1262x917.csv',     None),
    ('HeQ2017 (CN, MPA3)',            'HeQ2017_samples_122x.csv',      None),
    ('SG90 (SG, MPA4)',               'SG90_species_abundances.csv',   None),
    ('SG90 77+ (SG, MPA4)',           'SG90_species_abundances.csv',   'age77'),
    ('Villages (HN, MPA4)',           'Village_samples_1787x1910.csv', None),
    ('Yachida2019 (JP, MPA4)',        'Yachida2019_samples_615x1835.csv', None),
    ('Wave3 65+ (CN, mOTU)',          'Wave3_species_59x.csv',         None),
    ('Zhu2020 (CN, mOTU)',            'Zhu2020_samples_171x359.csv',   None),
    ('CHD Feng2016 (CN, mOTU)',       'CHD_species_102x.csv',          None),
]

AGE_CUT = 77


def load_all():
    """返回 {tag: (info, n)}，含 SG90 的 77 岁以上子集（Table 5 第 5 行）。"""
    import csv
    met = {}
    with open(os.path.join(BASE, 'SG90_metadata.csv'), encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            met[r['Subject ID']] = fnum(r['Age (in years)'])
    old = {k for k, a in met.items() if a >= AGE_CUT}

    out = {}
    for tag, fn, keep in SPEC:
        path = os.path.join(BASE, fn)
        if not os.path.exists(path):
            print('%-28s [缺文件] %s' % (tag, fn), flush=True)
            continue
        ids = old if keep == 'age77' else None
        info, n, _ = load_matrix(path, ids)
        out[tag] = (info, n)
        print('%-28s n=%4d species=%4d pool=%3d'
              % (tag, n, len(info), len(pool_of(info))), flush=True)
    return out


def pool_of(info):
    return [x for x in info if x['prev'] >= 0.20 and x['mean'] > 0]


def hstar_seeds(car_vals, pool, n_sp, n_ind):
    """四个种子各算一次，返回 (中位数, 最小值, 最大值, 原始值列表)。"""
    vals = [rf.h_cond(rf.families(car_vals, pool, n_sp, n_ind, seed=s)) for s in SEEDS]
    return statistics.median(vals), min(vals), max(vals), vals


def hdca_state(info, carr, ns_median, ns_single):
    """hdcA 家族状态。H 与种子无关；阈值分别用中位数表与单种子表各算一次。"""
    hit = match_hdca(info, carr)
    if not hit:
        return None
    n_ind = len(info[0]['vals'])
    tot = sum(x['mean'] for _, x, _ in hit) or 1.0
    ws = [x['mean'] / tot for _, x, _ in hit]
    cs = [c for _, _, c in hit]
    c_eff = sum(w * c for w, c in zip(ws, cs))
    H = sum(w * (c - c_eff) ** 2 for w, c in zip(ws, cs))
    thr_m, exact = interp(ns_median, len(hit))
    thr_s, _ = interp(ns_single, len(hit))
    p = [sum(x['vals'][i] for _, x, _ in hit) for i in range(n_ind)]
    c = [sum(x['vals'][i] * f for _, x, f in hit) for i in range(n_ind)]
    pos = [i for i, v in enumerate(p) if v > 0]
    return dict(matched=len(hit), c_eff=c_eff, H=H,
                thr_m=thr_m, thr_s=thr_s, ratio_m=H / thr_m, ratio_s=H / thr_s,
                exact=exact, n_ind=n_ind,
                rho_all=rf.spearman(p, c),
                rho_pos=rf.spearman([p[i] for i in pos], [c[i] for i in pos])
                if len(pos) > 5 else float('nan'))


def run_t5():
    carr, car_vals, _, _, _, _ = rf.load()
    loaded = load_all()

    rows = []
    print('\n=== 四种子中位数 ===', flush=True)
    for tag, (info, n) in loaded.items():
        pool = pool_of(info)
        ns_med, ns_single, detail = {}, {}, {}
        t0 = time.time()
        for n_sp in NS_LIST:
            if len(pool) < n_sp:
                continue
            med, lo, hi, vals = hstar_seeds(car_vals, pool, n_sp, n)
            rng = (hi - lo) / med * 100
            ns_med[n_sp] = med
            ns_single[n_sp] = vals[0]          # 种子 0 = 旧表口径
            detail[n_sp] = (med, lo, hi, rng, vals)
            print('  %-28s n_sp=%2d  median=%.5f  range=%4.1f%%  seed0=%.5f'
                  % (tag, n_sp, med, rng, vals[0]), flush=True)
        st = hdca_state(info, carr, ns_med, ns_single)
        spans = (max(ns_med.values()) / min(ns_med.values())) if ns_med else float('nan')
        spans_s = (max(ns_single.values()) / min(ns_single.values())) if ns_single else float('nan')
        rows.append((tag, n, len(info), len(pool), ns_med, ns_single, detail, st, spans, spans_s))
        print('  -> 用时 %.0fs  span(median)=%.2f  span(seed0)=%.2f'
              % (time.time() - t0, spans, spans_s), flush=True)

    out = os.path.join(BASE, 'cross_cohort9_seed_t5.tsv')
    with open(out, 'w', encoding='utf-8') as f:
        head = ['cohort', 'n_ind', 'species', 'pool']
        head += ['H%d_med' % k for k in NS_LIST] + ['H%d_seed0' % k for k in NS_LIST]
        head += ['span_med', 'span_seed0', 'matched', 'c_eff', 'H',
                 'thr_med', 'thr_seed0', 'ratio_med', 'ratio_seed0',
                 'rho_all', 'rho_carriers']
        f.write('\t'.join(head) + '\n')
        for tag, n, nsp, npool, ns_med, ns_single, detail, st, sp, sps in rows:
            row = [tag, str(n), str(nsp), str(npool)]
            row += ['%.5f' % ns_med[k] if k in ns_med else '' for k in NS_LIST]
            row += ['%.5f' % ns_single[k] if k in ns_single else '' for k in NS_LIST]
            row += ['%.2f' % sp, '%.2f' % sps]
            if st:
                row += [str(st['matched']), '%.3f' % st['c_eff'], '%.4f' % st['H'],
                        '%.5f' % st['thr_m'], '%.5f' % st['thr_s'],
                        '%.2f' % st['ratio_m'], '%.2f' % st['ratio_s'],
                        '%.3f' % st['rho_all'], '%.3f' % st['rho_pos']]
            f.write('\t'.join(row) + '\n')

    with open(os.path.join(BASE, 'cross_cohort9_seed_t5_detail.txt'), 'w', encoding='utf-8') as f:
        for tag, n, nsp, npool, ns_med, ns_single, detail, st, sp, sps in rows:
            f.write('%s  n=%d pool=%d\n' % (tag, n, npool))
            for k in NS_LIST:
                if k in detail:
                    med, lo, hi, rng, vals = detail[k]
                    f.write('  n_sp=%2d  %s  median=%.5f  [%.5f-%.5f]  range=%.1f%%\n'
                            % (k, ' '.join('%.5f' % v for v in vals), med, lo, hi, rng))
            f.write('\n')

    # --- 汇总：跨队列跨度、判定分布、Spearman ---
    print('\n=== 汇总 ===', flush=True)
    for k in NS_LIST:
        xs = [(tag, ns_med[k]) for tag, _, _, _, ns_med, _, _, _, _, _ in rows if k in ns_med]
        if len(xs) > 1:
            vs = [v for _, v in xs]
            print('n_sp=%2d  H* 中位数 %d 队列：%.5f-%.5f  跨度 %.2f 倍'
                  % (k, len(vs), min(vs), max(vs), max(vs) / min(vs)), flush=True)
    for label, key in (('中位数口径', 'ratio_m'), ('单种子口径', 'ratio_s')):
        pairs = [(r[0], r[7][key], r[7]['c_eff']) for r in rows if r[7]]
        nfail = sum(1 for _, v, _ in pairs if v > 1)
        npas = len(pairs) - nfail
        rho = rf.spearman([v for _, v, _ in pairs], [c for _, _, c in pairs])
        print('%s：fail %d / pass %d，Spearman(ratio, c_eff) = %.3f'
              % (label, nfail, npas, rho), flush=True)
        for tag, v, c in sorted(pairs, key=lambda x: x[1]):
            print('    %-28s ratio=%.2f  c_eff=%.3f  %s'
                  % (tag, v, c, 'FAIL' if v > 1 else 'pass'), flush=True)
    print('\n写出 %s' % out, flush=True)


def run_t6():
    import random
    carr, car_vals, _, _, _, _ = rf.load()
    rnd = random.Random(20261029)
    loaded = load_all()
    nm = min(v[1] for v in loaded.values())
    print('=== n=%d 对齐，%d 抽样 × %d 种子 ===' % (nm, N_DRAW, len(SEEDS)), flush=True)

    out_rows = []
    for tag, (info, n) in loaded.items():
        t0 = time.time()
        per_draw_med = {k: [] for k in NS_LIST}
        per_draw_seed0 = {k: [] for k in NS_LIST}
        for i in range(N_DRAW):
            idx = rnd.sample(range(n), nm)
            sub = subset(info, idx)
            pool = pool_of(sub)
            for k in NS_LIST:
                if len(pool) < k:
                    continue
                med, lo, hi, vals = hstar_seeds(car_vals, pool, k, nm)
                per_draw_med[k].append(med)
                per_draw_seed0[k].append(vals[0])
        row = {'tag': tag}
        for k in NS_LIST:
            if per_draw_med[k]:
                row[k] = (statistics.median(per_draw_med[k]),
                          min(per_draw_med[k]), max(per_draw_med[k]),
                          statistics.median(per_draw_seed0[k]))
        out_rows.append(row)
        line = '%-28s' % tag
        for k in NS_LIST:
            if k in row:
                m, lo, hi, s0 = row[k]
                line += '  H*(%d)=%.4f[%.4f-%.4f] seed0med=%.4f' % (k, m, lo, hi, s0)
        print(line + '  %.0fs' % (time.time() - t0), flush=True)

    out = os.path.join(BASE, 'cross_cohort9_seed_t6.tsv')
    with open(out, 'w', encoding='utf-8') as f:
        head = ['cohort']
        for k in NS_LIST:
            head += ['H%d_med' % k, 'H%d_lo' % k, 'H%d_hi' % k, 'H%d_seed0med' % k]
        f.write('\t'.join(head) + '\n')
        for row in out_rows:
            cells = [row['tag']]
            for k in NS_LIST:
                if k in row:
                    m, lo, hi, s0 = row[k]
                    cells += ['%.5f' % m, '%.5f' % lo, '%.5f' % hi, '%.5f' % s0]
                else:
                    cells += ['', '', '', '']
            f.write('\t'.join(cells) + '\n')
    print('\n写出 %s' % out, flush=True)


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 't5'
    t = time.time()
    if mode == 't5':
        run_t5()
    elif mode == 't6':
        run_t6()
    else:
        raise SystemExit('用法: python cross_cohort9_seed.py t5|t6')
    print('总用时 %.0f 秒' % (time.time() - t))
