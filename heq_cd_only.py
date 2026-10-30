# -*- coding: utf-8 -*-
"""
heq_cd_only.py —— 广州队列只保留克罗恩病（CD）患者后的重算。

背景：HeQ_2017 的公开矩阵有 122 份样本，原稿按"全部为 CD 患者"描述。核对
curatedMetagenomicData 的整理元数据（`HeQ_2017_sample_metadata.tsv`）后确认：
49 份 control == "Case"（CD），53 份 "Study Control"（健康对照），20 份未标注。
原稿的描述与原文（49 CD + 54 对照 + 14 例治疗后复查）不符，故本脚本按
"仅保留 CD 患者"重算该行。

与 cross_cohort9_seed.py 完全同一套口径：同一 pool 规则（prev>=0.20 且 mean>0）、
同一 hdcA 匹配、同一八个种子、同一对数插值。只把该队列的样本集合换掉。

输出：heq_cd_only.txt
"""
import csv
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import review_final as rf
from cross_cohort9 import load_matrix, interp, NS_LIST
from cross_cohort9_seed import SEEDS, pool_of, hstar_seeds, hdca_state

BASE = os.path.dirname(os.path.abspath(__file__))
META = os.path.join(BASE, 'HeQ_2017_sample_metadata.tsv')
MATRIX = os.path.join(BASE, 'HeQ2017_samples_122x.csv')


def cd_ids():
    """49 个 `control == Case` 的样本。

    优先读整理元数据；该文件是第三方发布的原始整理表，随包并不分发，
    故缺失时退回随包附带的 `HeQ2017_CD_ids.txt`（同一组 ID 的落盘形式）。
    """
    if os.path.exists(META):
        with open(META, encoding='utf-8') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        return {r['sample_id'] for r in rows if r['control'] == 'Case'}
    with open(os.path.join(BASE, 'HeQ2017_CD_ids.txt'), encoding='utf-8') as f:
        return {x.strip() for x in f if x.strip()}


def main():
    carr, car_vals, _, _, _, _ = rf.load()
    ids = cd_ids()

    info, n, _ = load_matrix(MATRIX, ids)
    pool = pool_of(info)
    lines = []
    lines.append('=== 广州队列：仅 CD 患者 ===')
    lines.append('样本 %d 份（元数据 control == Case）| 物种 %d | pool %d'
                 % (n, len(info), len(pool)))
    lines.append('SEEDS = %d 个，与主表同一组' % len(SEEDS))

    ns_med, ns_single, detail = {}, {}, {}
    t0 = time.time()
    for n_sp in NS_LIST:
        if len(pool) < n_sp:
            lines.append('  n_sp=%2d  [pool 不足]' % n_sp)
            continue
        med, lo, hi, vals = hstar_seeds(car_vals, pool, n_sp, n)
        rng = (hi - lo) / med * 100
        ns_med[n_sp] = med
        ns_single[n_sp] = vals[0]
        detail[n_sp] = (med, lo, hi, rng, vals)
        lines.append('  n_sp=%2d  %s  median=%.5f  [%.5f-%.5f]  range=%.1f%%'
                     % (n_sp, ' '.join('%.5f' % v for v in vals), med, lo, hi, rng))

    st = hdca_state(info, carr, ns_med, ns_single)
    span = max(ns_med.values()) / min(ns_med.values()) if ns_med else float('nan')
    lines.append('')
    lines.append('span(median) = %.2f   span(seed0) = %.2f'
                 % (span, max(ns_single.values()) / min(ns_single.values())))
    if st:
        lines.append('hdcA: 命中 %d 个类群 | c_eff=%.3f | H=%.4f'
                     % (st['matched'], st['c_eff'], st['H']))
        lines.append('      thr_med(%.0f)=%.5f  ratio_med=%.2f | thr_seed0=%.5f  ratio_seed0=%.2f'
                     % (st['matched'], st['thr_m'], st['ratio_m'], st['thr_s'], st['ratio_s']))
        lines.append('      rho_all=%.3f  rho_carriers=%.3f  n_ind=%d'
                     % (st['rho_all'], st['rho_pos'], st['n_ind']))
    lines.append('')
    lines.append('TSV 行（与 cross_cohort9_seed_t5.tsv 同列）：')
    row = ['HeQ2017 (CN, MPA3)', str(n), str(len(info)), str(len(pool))]
    row += ['%.5f' % ns_med.get(k, float('nan')) if k in ns_med else '' for k in NS_LIST]
    row += ['%.5f' % ns_single.get(k, float('nan')) if k in ns_single else '' for k in NS_LIST]
    row += ['%.2f' % span,
            '%.2f' % (max(ns_single.values()) / min(ns_single.values()))]
    if st:
        row += [str(st['matched']), '%.3f' % st['c_eff'], '%.4f' % st['H'],
                '%.5f' % st['thr_m'], '%.5f' % st['thr_s'],
                '%.2f' % st['ratio_m'], '%.2f' % st['ratio_s'],
                '%.3f' % st['rho_all'], '%.3f' % st['rho_pos']]
    lines.append('\t'.join(row))
    lines.append('用时 %.0f 秒' % (time.time() - t0))

    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(BASE, 'heq_cd_only.txt'), 'w', encoding='utf-8') as f:
        f.write(txt + '\n')

    # ---- 并入 Table 5：只替换广州一行，其余九行原样搬运 ----
    src = os.path.join(BASE, 'cross_cohort9_seed_t5.tsv')
    dst = os.path.join(BASE, 'cross_cohort9_seed_t5_heqcd.tsv')
    with open(src, encoding='utf-8') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        fields = rdr.fieldnames
        rows = list(rdr)
    changed = 0
    for r in rows:
        if r['cohort'].startswith('HeQ2017'):
            for k, v in zip(fields, row):
                r[k] = v
            changed += 1
    assert changed == 1, '表中应有且仅有一行 HeQ2017'
    with open(dst, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    lines2 = ['', '已写出 %s（%d 行，仅广州一行被替换，其余九行与表格逐字节相同）'
              % (os.path.basename(dst), len(rows))]
    print('\n'.join(lines2))
    with open(os.path.join(BASE, 'heq_cd_only.txt'), 'a', encoding='utf-8') as f:
        f.write('\n'.join(lines2) + '\n')


if __name__ == '__main__':
    main()
