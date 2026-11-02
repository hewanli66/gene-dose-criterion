# -*- coding: utf-8 -*-
"""生成 Table S4：胆盐水解酶（Pfam PF02275）的株水平携带率目录。

输入：../阳性对照/domain_carriage.tsv（由 domain_carriage.py 产出）
输出：TableS4_bsh_carriage.tsv 与 .md

表的用途是让 H 可被读者独立复算：给出每个物种的 c_s、在 ELDERMET 中的平均
相对丰度、归一化权重，以及该物种对 H 的贡献。贡献列之和即 H。

运行环境：Python 3.13，仅标准库。
"""
import csv
import os

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, '..', '阳性对照', 'domain_carriage.tsv')
DOMAIN = 'PF02275'

rows = [r for r in csv.DictReader(open(SRC, encoding='utf-8'), delimiter='\t')
        if r['domain'] == DOMAIN]

# 只在 c_s>0 的物种（即 catalog 判定的携带者）上归一化，与 domain_carriage.py 一致
carriers = [r for r in rows if float(r['c_s']) > 0]
wsum = sum(float(r['mean_abundance']) for r in carriers)
c_eff = sum(float(r['mean_abundance']) / wsum * float(r['c_s']) for r in carriers)
H = sum(float(r['mean_abundance']) / wsum * (float(r['c_s']) - c_eff) ** 2
        for r in carriers)

out = []
for r in sorted(rows, key=lambda x: -float(x['mean_abundance'])):
    n_p, n_s = int(r['n_present']), int(r['n_strains'])
    c_s = float(r['c_s'])
    w = float(r['mean_abundance']) / wsum if c_s > 0 else 0.0
    share = w * (c_s - c_eff) ** 2
    # n_s=1 的物种无法给出携带率（只能 0 或 1），标记出来，不计入 H
    # 只有 1 株的物种给不出携带率（只能取 0 或 1），标为 n/a 且不计入 H
    if n_s < 3:
        out.append([r['species'], n_p, n_s, 'n/a',
                    '%.6f' % float(r['mean_abundance']), 'n/a', 'n/a',
                    'single strain; not estimable'])
        continue
    out.append([r['species'], n_p, n_s, '%.4f' % c_s,
                '%.6f' % float(r['mean_abundance']), '%.4f' % w,
                '%.4f' % share, ''])

tsv = os.path.join(BASE, 'TableS4_bsh_carriage.tsv')
with open(tsv, 'w', encoding='utf-8', newline='') as f:
    w_ = csv.writer(f, delimiter='\t')
    w_.writerow(['species', 'strains_carrying', 'strains_tested', 'c_s',
                 'mean_relative_abundance_ELDERMET', 'weight', 'contribution_to_H', 'note'])
    w_.writerows(out)

md = os.path.join(BASE, 'TableS4_bsh_carriage.md')
with open(md, 'w', encoding='utf-8') as f:
    f.write('**Table S4. Strain-level carriage of bile salt hydrolase (Pfam PF02275) '
            'in twelve species of the ELDERMET cohort.**\n\n')
    f.write('| Species | Carrying / tested | *c*<sub>*s*</sub> | Mean rel. abundance | '
            'Weight | Contribution to *H* |\n|---|---|---|---|---|---|\n')
    for r in out:
        if r[3] == 'n/a':
            f.write('| *%s* | %d/%d | n/a | %.5f | n/a | n/a |\n'
                    % (r[0], r[1], r[2], float(r[4])))
        else:
            f.write('| *%s* | %d/%d | %.4f | %.5f | %.4f | %.4f |\n'
                    % (r[0], r[1], r[2], float(r[3]), float(r[4]), float(r[5]), float(r[6])))
    f.write('\n*H* = %.4f, *c*<sub>eff</sub> = %.4f, *n*<sub>sp</sub> = %d. '
            'Weights are relative abundances renormalised over carriers; the contribution '
            'column is the summand of *H* and the column sums to *H*. A species represented '
            'by a single genome cannot yield a carriage rate (the estimate can only be 0 or 1) '
            'and is marked n/a; it contributes to neither *c*<sub>eff</sub> nor *H*.\n'
            % (H, c_eff, len(carriers)))

print('c_eff = %.4f   H = %.4f   n_sp = %d' % (c_eff, H, len(carriers)))
print('written:', tsv)
print('written:', md)
