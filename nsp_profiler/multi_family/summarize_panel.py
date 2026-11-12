# -*- coding: utf-8 -*-
"""汇总多家族面板的判决分布，回答「判据能否在 hdcA/BSH 之外复现」。

输入：out/families_panel.tsv（run_panel_mp.py 产出）
输出：out/panel_summary.txt + out/panel_fail.tsv（fail 家族清单）

口径提醒：n_sp 很小时（<5）H* 落在 §5.4 插值的下限 0.0240，判决敏感；
正文 §3.3 已声明 n_sp<5 只作筛选不作裁定，本汇总同样分列。
"""
import csv
import math
import os
import statistics
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, 'out')


def pct(v, q):
    s = sorted(v)
    return s[min(len(s) - 1, int(q * len(s)))]


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else 'panel'
    path = os.path.join(OUT, 'families_%s.tsv' % tag)
    rows = list(csv.DictReader(open(path, encoding='utf-8'), delimiter='\t'))
    rows = [r for r in rows if int(r['n_sp']) > 0]
    print('家族 %d' % len(rows))

    ratio = [float(r['H_over_Hstar']) for r in rows]
    H = [float(r['H']) for r in rows]
    ceff = [float(r['c_eff']) for r in rows]
    lam = [float(r['lambda_det']) for r in rows]

    lines = []
    lines.append('家族总数 %d' % len(rows))
    lines.append('H      中位 %.4f  p90 %.4f  最大 %.4f' %
                 (statistics.median(H), pct(H, .9), max(H)))
    lines.append('H/H*   中位 %.3f  p90 %.3f  最大 %.3f' %
                 (statistics.median(ratio), pct(ratio, .9), max(ratio)))
    lines.append('c_eff  中位 %.4f  p10 %.4f  最小 %.4f' %
                 (statistics.median(ceff), pct(ceff, .1), min(ceff)))
    lines.append('λ_det  中位 %.1f  p10 %.1f  最小 %.1f' %
                 (statistics.median(lam), pct(lam, .1), min(lam)))

    for lo, name in ((0, '全部'), (5, 'n_sp≥5（正文可裁定的区间）')):
        sub = [r for r in rows if int(r['n_sp']) >= lo]
        if not sub:
            continue
        f = sum(1 for r in sub if float(r['H_over_Hstar']) > 1)
        lines.append('')
        lines.append('[%s] 家族 %d：fail(H/H*>1) %d（%.1f%%），pass %d'
                     % (name, len(sub), f, 100 * f / len(sub), len(sub) - f))
        # 双约束四象限（§2.5 的 λ_det ≥ 5；§2.6 的四态）
        q = {'ok': 0, 'low_lam': 0, 'low_fid': 0, 'both': 0}
        for r in sub:
            bad_h = float(r['H_over_Hstar']) > 1
            bad_l = float(r['lambda_det']) < 5
            if not bad_h and not bad_l:
                q['ok'] += 1
            elif bad_h and bad_l:
                q['both'] += 1
            elif bad_h:
                q['low_fid'] += 1
            else:
                q['low_lam'] += 1
        n = len(sub)
        lines.append('  双约束：可用 %d（%.1f%%）| 仅保真度失败 %d | 仅可检测性失败 %d | 两者皆失败 %d'
                     % (q['ok'], 100 * q['ok'] / n, q['low_fid'], q['low_lam'], q['both']))

    # 按 n_sp 分组
    lines.append('')
    lines.append('按携带物种数：')
    lines.append('%-6s %6s %10s %10s %10s' % ('n_sp', 'n', 'medH', 'medH/H*', 'fail%'))
    by = {}
    for r in rows:
        by.setdefault(int(r['n_sp']), []).append(r)
    for k in sorted(by):
        if len(by[k]) < 5:
            continue
        g = by[k]
        hh = [float(x['H']) for x in g]
        rr = [float(x['H_over_Hstar']) for x in g]
        f = sum(1 for x in rr if x > 1)
        lines.append('%-6d %6d %10.4f %10.2f %9.1f%%'
                     % (k, len(g), statistics.median(hh), statistics.median(rr),
                        100 * f / len(g)))

    # fail 清单
    fails = sorted([r for r in rows if float(r['H_over_Hstar']) > 1],
                   key=lambda r: -float(r['H_over_Hstar']))
    p_fail = os.path.join(OUT, 'panel_fail.tsv')
    with open(p_fail, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['pfam', 'n_sp', 'c_eff', 'H', 'H_star', 'H_over_Hstar',
                    'lambda_det', 'L_gene', 'carrier_abundance'])
        for r in fails:
            w.writerow([r['pfam'], r['n_sp'], r['c_eff'], r['H'], r['H_star'],
                        r['H_over_Hstar'], r['lambda_det'], r['L_gene'],
                        r['carrier_abundance']])
    lines.append('')
    lines.append('fail 家族 %d 个，前 15 名：' % len(fails))
    for r in fails[:15]:
        lines.append('  %-9s n_sp=%-2s c_eff=%s H=%s H/H*=%s λ_det=%s'
                     % (r['pfam'], r['n_sp'], r['c_eff'], r['H'],
                        r['H_over_Hstar'], r['lambda_det']))

    txt = '\n'.join(lines)
    with open(os.path.join(OUT, 'panel_summary.txt'), 'w', encoding='utf-8') as f:
        f.write(txt + '\n')
    print(txt)
    print('\n写出\n  %s\n  %s' % (os.path.join(OUT, 'panel_summary.txt'), p_fail))


if __name__ == '__main__':
    main()
