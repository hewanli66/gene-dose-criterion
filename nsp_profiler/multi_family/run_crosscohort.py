# -*- coding: utf-8 -*-
"""跨队列外部验证：同一批 c_s（株水平携带率），换一个队列的丰度权重。

为什么要这一步
--------------
原面板的 538 个家族判决只用了一个队列（ELDERMET，189 位老年人）的丰度。
判据的输入有两半：c_s 来自基因组（队列无关），p_s 来自队列（队列相关）。
换队列只动 p_s 这一半，检验的是「判决对队列组成的敏感性」——正是 §2.9 对
hdcA 做的跨九队列检验，现在对 538 个家族一次做完。

数据
----
c_s       : out/carriage_panel.tsv（2026-11-11 面板扫描产出，132 株 × 家族）
ELDERMET  : 同文件的 mean_abundance 列
HMP2      : experiments-p2-p3/humann3/hmp2_h3/taxonomic_profiles_3.tsv.gz
            MetaPhlAn 3，1,638 个样本，单位 %，取物种行全样本算术平均

口径（与 run_panel_mp.py 一致，只有一处必须改）
------------------------------------------------
n_sp = 携带者中**平均丰度非零**的物种数。正文 §2.6 就是这样定的（hdcA「only
five of the 20 contribute non-zero abundance, so the appropriate threshold is
that for n_sp = 5」）。原脚本用 c_s>0 计数，在 ELDERMET 上等价（top 物种丰度
全非零），跨队列必须改成非零丰度，否则一个在本队列中缺席的物种会撑大 n_sp，
把 H* 抬高（更宽松）却对 c_eff/H 毫无贡献。

用法：python run_crosscohort.py [--min-ab 0.001] [--tag hmp2]
"""
import argparse
import csv
import gzip
import math
import os

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))
HMP2 = os.path.join(ROOT, 'experiments-p2-p3', 'humann3', 'hmp2_h3',
                    'taxonomic_profiles_3.tsv.gz')
PANEL = os.path.join(BASE, 'out', 'carriage_panel.tsv')
FAM = os.path.join(BASE, 'out', 'families_panel.tsv')
OUT = os.path.join(BASE, 'out')

HSTAR = {3: 0.0240, 5: 0.0378, 10: 0.0732, 20: 0.0936}
READS = 6.3e6      # HMP2 深度未取到，沿用 ELDERMET 的 6.3 M，仅影响 λ 绝对值
EPS = 0.71


def hstar_at(n):
    ks = sorted(HSTAR)
    if n <= ks[0]:
        return HSTAR[ks[0]]
    if n >= ks[-1]:
        return HSTAR[ks[-1]]
    for a, b in zip(ks, ks[1:]):
        if a <= n <= b:
            t = (math.log(n) - math.log(a)) / (math.log(b) - math.log(a))
            return math.exp(math.log(HSTAR[a]) * (1 - t) + math.log(HSTAR[b]) * t)
    return HSTAR[ks[-1]]


def hmp2_abundance():
    """物种名（下划线）-> 全样本平均相对丰度（%）。"""
    out = {}
    with gzip.open(HMP2, 'rt', encoding='utf-8', errors='replace') as f:
        r = csv.reader(f, delimiter='\t')
        next(r)
        for row in r:
            if '|s__' not in row[0]:
                continue
            nm = row[0].split('|s__')[-1]
            vals = []
            for v in row[1:]:
                try:
                    vals.append(float(v))
                except ValueError:
                    pass
            out[nm] = sum(vals) / len(vals) if vals else 0.0
    return out


def judge(rows, ab, min_ab):
    """rows: [species, c_s, L_bp, G_bp]；返回判决元组或 None。"""
    keep = [r for r in rows if ab.get(r[0], 0.0) > min_ab]
    if not keep:
        return None
    wsum = sum(ab[r[0]] for r in keep)
    ceff = sum(ab[r[0]] / wsum * r[1] for r in keep)
    H = sum(ab[r[0]] / wsum * (r[1] - ceff) ** 2 for r in keep)
    nsp = len(keep)
    hs = hstar_at(nsp)
    lam = sum(READS * ab[r[0]] / 100.0 * r[1] * r[2] / r[3]
              for r in keep if r[3])
    return nsp, ceff, H, hs, H / hs, lam, lam * EPS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-ab', type=float, default=0.001,
                    help='物种视为「在本队列中存在」的平均丰度下限（%%）')
    ap.add_argument('--tag', default='hmp2')
    a = ap.parse_args()

    ab2 = {k.replace('_', ' '): v for k, v in hmp2_abundance().items()}
    byfam = {}
    with open(PANEL, encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            byfam.setdefault(r['pfam'], []).append(
                [r['species'], float(r['c_s']),
                 float(r['median_hit_len']) * 3, float(r['genome_bp'])])

    eld = {}
    with open(FAM, encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            eld[r['pfam']] = r

    out = []
    for pf in sorted(byfam):
        e = eld.get(pf)
        if not e:
            continue
        got = judge(byfam[pf], ab2, a.min_ab)
        if not got:
            continue
        nsp, ceff, H, hs, ratio, lam, lam_det = got
        e_nsp = int(e['n_sp'])
        e_ratio = float(e['H_over_Hstar'])
        e_det = float(e['lambda_det'])
        v_e = 'pass' if (e_ratio <= 1 and e_det >= 5) else (
            'fail_fid' if e_ratio > 1 and e_det >= 5 else (
                'fail_det' if e_ratio <= 1 else 'fail_both'))
        v_2 = 'pass' if (ratio <= 1 and lam_det >= 5) else (
            'fail_fid' if ratio > 1 and lam_det >= 5 else (
                'fail_det' if ratio <= 1 else 'fail_both'))
        out.append([pf, e_nsp, e['c_eff'], e['H'], e['H_over_Hstar'],
                    '%.2f' % e_det, v_e, nsp, '%.4f' % ceff, '%.4f' % H,
                    '%.3f' % ratio, '%.2f' % lam_det, v_2,
                    'same' if v_e == v_2 else 'flip'])

    p = os.path.join(OUT, 'crosscohort_%s.tsv' % a.tag)
    with open(p, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['pfam', 'n_sp_eldermet', 'c_eff_eldermet', 'H_eldermet',
                    'ratio_eldermet', 'lambda_det_eldermet', 'verdict_eldermet',
                    'n_sp_hmp2', 'c_eff_hmp2', 'H_hmp2', 'ratio_hmp2',
                    'lambda_det_hmp2', 'verdict_hmp2', 'change'])
        w.writerows(out)

    n = len(out)
    flip = [r for r in out if r[-1] == 'flip']
    print('共同家族 %d，判决翻转 %d（%.1f%%）' % (n, len(flip), 100 * len(flip) / n))
    from collections import Counter
    print('ELDERMET 判决:', Counter(r[6] for r in out).most_common())
    print('HMP2     判决:', Counter(r[12] for r in out).most_common())
    print('转移矩阵（行=ELDERMET，列=HMP2）:')
    vs = ['pass', 'fail_fid', 'fail_det', 'fail_both']
    print('%-12s' % '' + ''.join('%-11s' % v for v in vs))
    for v in vs:
        print('%-12s' % v + ''.join(
            '%-11d' % sum(1 for r in out if r[6] == v and r[12] == v2) for v2 in vs))
    k = [r for r in out if r[6] == 'fail_fid']
    k2 = [r for r in k if r[12] == 'fail_fid']
    print('\nELDERMET 保真度失败的 %d 个家族，在 HMP2 仍失败 %d 个（%.0f%%）'
          % (len(k), len(k2), 100 * len(k2) / max(1, len(k))))
    print('写入 %s' % p)


if __name__ == '__main__':
    main()
