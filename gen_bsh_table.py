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

# 敏感性检验：只有 c_s<1 的物种可能被基因组抽样改变，逐个置 0 看判决是否翻转。
# H* = 0.0757 来自 §5.4 在 n_sp = 11 处的 log-log 插值（与正文一致，此处只做比值）。
# 置零一个物种等于把它从携带者里去掉，n_sp 由 11 降为 10，故敏感性比值须用
# H*(10) = 0.0757 之外的另一格：Table 5 第 10 物种格给的 0.0732。
HSTAR = 0.0757
HSTAR_N10 = 0.0732
_cmap = {r['species']: float(r['c_s']) for r in carriers}
_wmap = {r['species']: float(r['mean_abundance']) / wsum for r in carriers}


def H_with(over):
    cc = dict(_cmap)
    cc.update(over)
    ce = sum(_wmap[k] * cc[k] for k in cc)
    return ce, sum(_wmap[k] * (cc[k] - ce) ** 2 for k in cc)

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
    _, h_fp = H_with({'Faecalibacterium prausnitzii': 0.0})
    _, h_pc = H_with({'Prevotella copri': 0.0})
    f.write('\n*Sensitivity to the two species whose carriage is below 1.* Only *P. copri* '
            '(10/12) and *F. prausnitzii* (3/12) can be moved by the choice of genomes; the '
            'other nine are 12/12. Setting *c*<sub>*s*</sub> to 0 for one species also removes '
            'it from the count of carrying species, so *n*<sub>sp</sub> falls from 11 to 10 and '
            'both ratios below use *H*\\* = %.4f rather than the %.4f of the main comparison. '
            '*F. prausnitzii*, which supplies 88%% of *H*, gives *H* = %.4f (*H*/*H*\\* = %.2f); '
            '*P. copri* gives *H* = %.4f (*H*/*H*\\* = %.2f).\n'
            % (HSTAR_N10, HSTAR, h_fp, h_fp / HSTAR_N10, h_pc, h_pc / HSTAR_N10))

print('c_eff = %.4f   H = %.4f   n_sp = %d' % (c_eff, H, len(carriers)))
print('written:', tsv)
print('written:', md)
