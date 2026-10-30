# -*- coding: utf-8 -*-
"""
gen_supp_tables.py —— 生成 Supplementary Table S1 与 S2。

S1 是"四档阈值表"的完整版：Table 5 只给中位数，这里补上每一格的八种子范围
与极差百分比（正文 §2.4 与 Methods "Seed dependence" 引用的 35–123% / 41–127%
/ 12–52% / 6–35% 全部由这张表合成）。

S2 是 M11 敏感性网格：保真门槛 {0.85, 0.90, 0.95} × 保留率 {90%, 95%, 99%}，
ELDERMET、n_sp = 5，每格为八种子中位数，并给出该格的种子极值与 hdcA 倍数。

数据来源（全部为已落盘的运行结果，本脚本不重跑抽样）：
  cross_cohort9_seed_t5_detail.txt   九队列八种子明细（广州为旧的 122 份口径）
  heq_cd_only.txt                    广州 49 例 CD 口径，替换上一文件的广州行
  cross_cohort9_seed_t5_heqcd.tsv    十行汇总（n、物种数、pool、hdcA 状态）
  review_M11_grid_median.tsv         3×3 网格（八种子中位数）

输出：
  TableS1_threshold_detail.tsv       长表，40 行（十队列 × 四档）
  TableS1_threshold_detail.md        宽表，便于阅读与排版
  TableS2_M11_grid.tsv               9 行
  TableS2_M11_grid.md                3×3 面板 + hdcA 倍数面板
  supp_tables_check.txt              与 Table 5 逐格核对的结果

运行：python gen_supp_tables.py    （纯标准库，无需 numpy）
"""
import csv
import os
import re
import statistics
from decimal import Decimal, ROUND_HALF_UP

BASE = os.path.dirname(os.path.abspath(__file__))
NS = [3, 5, 10, 20]

# 明细文件里的队列名 → 汇总 TSV 里的队列名（顺序一致，显式写出便于审计）
NAME_MAP = [
    ('ELDERMET (IE, MPA3)', 'ELDERMET (Ireland)'),
    ('CRC9 pool (MPA3)', 'CRC collection, 9 cohorts'),
    ('HeQ2017 (CN, MPA3)', 'Guangzhou (China)'),          # 该行用 heq_cd_only.txt 的值
    ('SG90 (SG, MPA4)', 'SG90 (Singapore), all'),
    ('SG90 77+ (SG, MPA4)', 'SG90, aged ≥77'),
    ('Villages (HN, MPA4)', 'Villages (Honduras)'),
    ('Yachida2019 (JP, MPA4)', 'Screening cohort (Japan)'),
    ('Wave3 65+ (CN, mOTU)', 'Wave 3 (China), aged ≥65'),
    ('Zhu2020 (CN, mOTU)', "Xi'an cohort (China)"),
    ('CHD Feng2016 (CN, mOTU)', 'Coronary cohort (China)'),
]

PROFILER = {
    'ELDERMET (Ireland)': 'MetaPhlAn 3',
    'CRC collection, 9 cohorts': 'MetaPhlAn 3',
    'Guangzhou (China)': 'MetaPhlAn 3',
    'SG90 (Singapore), all': 'MetaPhlAn 4',
    'SG90, aged ≥77': 'MetaPhlAn 4',
    'Villages (Honduras)': 'MetaPhlAn 4',
    'Screening cohort (Japan)': 'MetaPhlAn 4',
    'Wave 3 (China), aged ≥65': 'mOTU',
    "Xi'an cohort (China)": 'mOTU',
    'Coronary cohort (China)': 'mOTU',
}


def parse_detail(path):
    """返回 {队列名: {n_sp: [八个值]}}"""
    out, cur = {}, None
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            if not line.startswith('  '):
                cur = line.split('  ')[0].strip()
                out[cur] = {}
                continue
            m = re.match(r'\s+n_sp=\s*(\d+)\s+(.*?)median=', line)
            if m and cur:
                out[cur][int(m.group(1))] = [float(x) for x in m.group(2).split()]
    return out


def parse_heq(path):
    d = {}
    with open(path, encoding='utf-8') as f:
        for line in f:
            m = re.match(r'\s+n_sp=\s*(\d+)\s+(.*?)median=', line)
            if m:
                d[int(m.group(1))] = [float(x) for x in m.group(2).split()]
    return d


def load_summary(path):
    with open(path, encoding='utf-8') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        return {r['cohort']: r for r in rdr}


def fmt_pct(x):
    return '%.0f' % x


def q4(x):
    """四舍五入到 4 位小数（half-up），避免 Python 默认的银行家舍入把
    0.09355 印成 0.0935、与 Table 5 的 0.0936 对不上。"""
    return str(Decimal(repr(float(x))).quantize(Decimal('0.0001'),
                                                rounding=ROUND_HALF_UP))


def main():
    detail = parse_detail(os.path.join(BASE, 'cross_cohort9_seed_t5_detail.txt'))
    heq = parse_heq(os.path.join(BASE, 'heq_cd_only.txt'))
    summ = load_summary(os.path.join(BASE, 'cross_cohort9_seed_t5_heqcd.tsv'))

    # ---- 组装 40 格 ----
    rows = []          # 长表
    wide = {}          # 宽表用
    for det_name, tsv_name in NAME_MAP:
        src = heq if det_name.startswith('HeQ') else detail[det_name]
        s = summ[det_name]
        H = float(s['H'])
        for nsp in NS:
            vals = src[nsp]
            # 中位数取汇总 TSV 的全精度值（明细文件只打印 5 位小数，重算会有末位差）
            med = float(s['H%d_med' % nsp])
            lo, hi = min(vals), max(vals)
            rng = (hi - lo) / med * 100
            rows.append(dict(
                cohort=tsv_name, n_ind=s['n_ind'], profiler=PROFILER[tsv_name],
                species=s['species'], pool=s['pool'], n_sp=nsp,
                hstar_median=q4(med), hstar_min=q4(lo),
                hstar_max=q4(hi), seed_range_pct='%.1f' % rng,
                n_seeds=len(vals),
                hdca_matched=s['matched'], c_eff=s['c_eff'], H='%.4f' % H,
                ratio='%.2f' % float(s['ratio_med']),
                rho_all=s['rho_all'], rho_carriers=s['rho_carriers'],
            ))
            wide.setdefault(tsv_name, {})[nsp] = (med, lo, hi, rng)

    # ---- 写 S1 长表 ----
    fields = ['cohort', 'n_ind', 'profiler', 'species', 'pool', 'n_sp',
              'hstar_median', 'hstar_min', 'hstar_max', 'seed_range_pct', 'n_seeds',
              'hdca_matched', 'c_eff', 'H', 'ratio', 'rho_all', 'rho_carriers']
    p1 = os.path.join(BASE, 'TableS1_threshold_detail.tsv')
    with open(p1, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    # ---- 写 S1 宽表（Markdown）----
    L = []
    L.append('**Supplementary Table S1. Conditional thresholds in ten cohort-settings, '
             'with the spread over eight seeds.**')
    L.append('')
    L.append('Every entry is `median (min–max; range as % of the median)` over the eight '
             'simulation seeds, at the family size given by the column. `H*` is the '
             'conditional threshold of §2.4: the largest `H` at which '
             'P(ρ ≥ 0.90 | H ≈ h) ≥ 0.95, with a moving window of width 0.02, '
             '`min_n` = 40 and 1,200 simulated families per level. Cohorts are the ten '
             'settings of Table 5: nine cohorts plus one age-restricted subset of SG90.')
    L.append('')
    L.append('| Cohort | *n* | Profiler | Species | Pool | *H*\\*(3) | *H*\\*(5) | *H*\\*(10) | *H*\\*(20) |')
    L.append('|---|---|---|---|---|---|---|---|---|')
    for det_name, tsv_name in NAME_MAP:
        s = summ[det_name]
        cells = []
        for nsp in NS:
            med, lo, hi, rng = wide[tsv_name][nsp]
            cells.append('%s (%s–%s; %s%%)' % (q4(med), q4(lo), q4(hi), fmt_pct(rng)))
        L.append('| %s | %s | %s | %s | %s | %s |' % (
            tsv_name, '{:,}'.format(int(s['n_ind'])), PROFILER[tsv_name],
            '{:,}'.format(int(s['species'])), s['pool'],
            ' | '.join(cells)))
    L.append('')

    def rng_of(nsp):
        sub = [float(r['seed_range_pct']) for r in rows if r['n_sp'] == nsp]
        return min(sub), max(sub)

    allr = [float(r['seed_range_pct']) for r in rows]
    parts = ['`n_sp` = %d, %.0f–%.0f%%' % (nsp, rng_of(nsp)[0], rng_of(nsp)[1])
             for nsp in NS]
    L.append('Pooled ranges across the ten settings, by family size: '
             + '; '.join(parts) + '. Over all %d cells the range is %.0f–%.0f%% with a '
             'median of %.1f%%. The `n_sp` = 3 column is weakly identified (it moves '
             '1.8-fold with the width of the moving window, against 1.1-fold at '
             '`n_sp` = 20; §2.4) and is reported as an order-of-magnitude anchor only. '
             'No comparison in the paper rests on it; the single use is the robustness '
             'bound in §2.4, where the *hdcA* family still fails at the largest of the '
             'eight seed values at that size.'
             % (len(allr), min(allr), max(allr), statistics.median(allr)))
    L.append('')
    L.append('The Guangzhou row uses the 49 samples that the curated metadata labels as '
             'cases (Methods §2.1.6d); its *H*\\*(3) exceeds its *H*\\*(5), which is why '
             'the span of that cohort is quoted as max/min across the four sizes rather '
             'than as a ratio of endpoints (Table 5, note).')
    p1m = os.path.join(BASE, 'TableS1_threshold_detail.md')
    with open(p1m, 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')

    # ---- S2 ----
    grid = []
    with open(os.path.join(BASE, 'review_M11_grid_median.tsv'), encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            grid.append(r)
    H_hdca = float(summ['ELDERMET (IE, MPA3)']['H'])

    p2 = os.path.join(BASE, 'TableS2_M11_grid.tsv')
    gfields = ['fidelity_threshold', 'keep_level', 'hstar_median', 'hstar_min',
               'hstar_max', 'seed_range_pct', 'hstar_range_factor',
               'hdca_ratio_median', 'hdca_ratio_worst_seed']
    out2 = []
    for r in grid:
        med = float(r['H_star_median']); lo = float(r['H_star_min']); hi = float(r['H_star_max'])
        out2.append({
            'fidelity_threshold': r['fidelity_threshold'],
            'keep_level': r['keep_rate'],
            'hstar_median': q4(med),
            'hstar_min': q4(lo),
            'hstar_max': q4(hi),
            'seed_range_pct': '%.1f' % ((hi - lo) / med * 100),
            'hstar_range_factor': '%.1f' % (hi / lo),
            'hdca_ratio_median': '%.2f' % (H_hdca / med),
            'hdca_ratio_worst_seed': '%.2f' % (H_hdca / hi),
        })
    with open(p2, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=gfields, delimiter='\t')
        w.writeheader()
        w.writerows(out2)

    meds = [float(r['hstar_median']) for r in out2]
    r_med = [float(r['hdca_ratio_median']) for r in out2]
    r_worst = [float(r['hdca_ratio_worst_seed']) for r in out2]

    M = []
    M.append('**Supplementary Table S2. Sensitivity of the threshold to the two '
             'arbitrary conventions of the criterion.**')
    M.append('')
    M.append('ELDERMET cohort, *n*<sub>sp</sub> = 5, *H*(*hdcA*) = %.4f. The fidelity '
             'cut-off is the ρ that a family must reach; the retention level is the '
             'probability with which it must do so. Each cell is the median of the eight '
             'seeds, as everywhere in the paper; the bracketed interval is the range '
             'over those seeds. `hdcA` × is *H*(*hdcA*) divided by that cell\'s median '
             'threshold.' % H_hdca)
    M.append('')
    M.append('| Fidelity cut-off | Retained 90% | Retained 95% | Retained 99% |')
    M.append('|---|---|---|---|')
    by = {(r['fidelity_threshold'], r['keep_level']): r for r in out2}
    for cut in ['0.85', '0.90', '0.95']:
        cells = []
        for keep in ['0.90', '0.95', '0.99']:
            r = by[(cut, keep)]
            cells.append('%.4f (%.4f–%.4f) ×%.1f' % (
                float(r['hstar_median']), float(r['hstar_min']),
                float(r['hstar_max']), float(r['hdca_ratio_median'])))
        M.append('| %s | %s |' % (cut, ' | '.join(cells)))
    M.append('')
    M.append('Across the nine cells the median threshold spans %.4f–%.4f, a %.1f-fold '
             'range, and the *hdcA* margin runs ×%s to ×%s — it exceeds unity in all '
             'nine. Taking, in each cell, the least favourable of the eight seeds '
             'instead of the median, the smallest margin in the grid is ×%s, so no '
             'choice of convention reverses the verdict; only the margin moves.'
             % (min(meds), max(meds), max(meds) / min(meds),
                '%.1f' % min(r_med), '%.1f' % max(r_med), '%.1f' % min(r_worst)))
    M.append('')
    M.append('The permissive corner (cut-off 0.95, retention 99%%) has both the lowest '
             'threshold and the widest seed range (%.4f–%.4f, %.1f-fold), which is why '
             'the grid under the median convention (0.0116–0.1087) is wider than the '
             'single-draw version first reported (0.021–0.117, 1.3–7.1-fold).'
             % (float(by[('0.95', '0.99')]['hstar_min']),
                float(by[('0.95', '0.99')]['hstar_max']),
                float(by[('0.95', '0.99')]['hstar_range_factor'])))
    p2m = os.path.join(BASE, 'TableS2_M11_grid.md')
    with open(p2m, 'w', encoding='utf-8') as f:
        f.write('\n'.join(M) + '\n')

    # ---- 与 Table 5 逐格核对 ----
    t5 = {
        'ELDERMET (Ireland)': [0.0240, 0.0378, 0.0732, 0.0936],
        'CRC collection, 9 cohorts': [0.0300, 0.0459, 0.0856, 0.1056],
        'Guangzhou (China)': [0.0584, 0.0413, 0.0669, 0.0883],
        'SG90 (Singapore), all': [0.0447, 0.0798, 0.1064, 0.0994],
        'SG90, aged ≥77': [0.0369, 0.0645, 0.1067, 0.0907],
        'Villages (Honduras)': [0.0285, 0.0443, 0.0766, 0.0921],
        'Screening cohort (Japan)': [0.0318, 0.0413, 0.0714, 0.1081],
        'Wave 3 (China), aged ≥65': [0.0270, 0.0390, 0.0834, 0.0940],
        "Xi'an cohort (China)": [0.0338, 0.0525, 0.1045, 0.1125],
        'Coronary cohort (China)': [0.0241, 0.0367, 0.0806, 0.0936],
    }
    lines, bad = [], 0
    for name, ref in t5.items():
        for i, nsp in enumerate(NS):
            med = wide[name][nsp][0]
            ok = q4(med) == q4(ref[i])
            if not ok:
                bad += 1
                lines.append('MISMATCH %-28s n_sp=%2d  S1=%s  Table5=%s'
                             % (name, nsp, q4(med), q4(ref[i])))
    lines.insert(0, 'Table S1 与 Table 5 逐格核对：%d 格，%d 处不一致'
                 % (len(t5) * 4, bad))
    for nsp in NS:
        sub = [r['seed_range_pct'] for r in rows if r['n_sp'] == nsp]
        lines.append('  n_sp=%2d 极差区间 %5.1f–%5.1f%%（十队列）'
                     % (nsp, min(float(x) for x in sub), max(float(x) for x in sub)))
    allr = [float(r['seed_range_pct']) for r in rows]
    lines.append('  40 格：%.1f–%.1f%%，中位 %.2f%%'
                 % (min(allr), max(allr), statistics.median(allr)))
    lines.append('S2 网格中位数 %.4f–%.4f（%.1f 倍）；hdcA 倍数 ×%.1f–×%.1f；'
                 '最不利单种子 ×%.1f'
                 % (min(meds), max(meds), max(meds) / min(meds),
                    min(r_med), max(r_med), min(r_worst)))
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(BASE, 'supp_tables_check.txt'), 'w', encoding='utf-8') as f:
        f.write(txt + '\n')
    print('\n已写出 %s / %s / %s / %s'
          % tuple(os.path.basename(p) for p in (p1, p1m, p2, p2m)))


if __name__ == '__main__':
    main()
