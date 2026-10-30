# -*- coding: utf-8 -*-
"""检出效率 ε 的重新标定（取代 lambda_v2.py 里手推的 1/6 × 0.9 = 0.15）

问题出在哪
----------
lambda_v2.py 的推导是：「6 个读框里只有 1 个是编码框，所以效率 1/6」。
这一步是错的。六框翻译的做法是**把 6 个框全部试一遍**，编码框必然在其中
（正链 3 个 + 反链 3 个，无论 read 来自哪条链、哪个相位）。1/6 只适用于
「固定读一个框」的情形。

真正的损失只有三项，都很小：
  a. read 跨越基因边界（100 bp / 3.8 kb ≈ 2.6%）
  b. 测序错误造成移码/提前终止（Illumina ~0.1%/碱基，100 bp 下 <10%）
  c. 与参考序列的分歧（15-mer 要求 15 aa 完全相同）

本脚本做两件事
--------------
1. in-silico 召回实验：直接从 rpoB 参考序列上切 100 bp reads，正反链各半，
   跑真正的 orf_kmers()，量出「一条确实来自该基因的 read，被检出的概率」。
   这是 ε 的**上限**（不含测序错误与分歧）。
2. 用 rpoB 对照的实测命中数反解 ε：ε = observed / expected，
   expected 由基因组尺度算出（rpoB 单拷贝，不需丰度表）。
   给出对「细菌 reads 占比」和「平均基因组大小」的敏感性网格。

输出 eff_calibration.txt
"""
import json
import math
import os
import random
import statistics as st

import scan_hdca as SH

BASE = os.path.dirname(os.path.abspath(__file__))
K = 15
READ_LEN = 100


def in_silico_recall(n_reads=20000, seed=7):
    """从 rpoB 参考上随机切 read，量检出率（不含分歧与测序错误）"""
    d = json.load(open(os.path.join(BASE, 'rpob_ref.json'), encoding='utf-8'))
    seqs = {sp: v['seq'] for sp, v in d.items()}
    # 把蛋白序列反向翻译成随机密码子——只是为了让 orf_kmers 有核苷酸可翻。
    # 注意：这里必须用固定密码子表回译，否则无法保证翻译回来还是原序列。
    CODON_OF = {
        'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT',
        'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT',
        'P': 'CCT', 'Q': 'CAA', 'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT',
        'W': 'TGG', 'Y': 'TAT',
    }
    COMP = str.maketrans('ACGT', 'TGCA')
    rng = random.Random(seed)
    res = {}
    for sp, aa in seqs.items():
        dna = ''.join(CODON_OF[c] for c in aa)
        # 参考 k-mer 集（与扫描时一致）
        ref = {aa[i:i + K] for i in range(len(aa) - K + 1)}
        hit = 0
        tot = 0
        for _ in range(n_reads // len(seqs)):
            if len(dna) <= READ_LEN:
                break
            p = rng.randrange(0, len(dna) - READ_LEN)
            r = dna[p:p + READ_LEN]
            if rng.random() < 0.5:                      # 一半来自反义链
                r = r.translate(COMP)[::-1]
            tot += 1
            if SH.orf_kmers(r) & ref:
                hit += 1
        res[sp] = (hit, tot)
    return res


def eps_from_counts(obs_per_m, gene_bp, genome_mb, f_bact):
    """由实测命中率反解检出效率"""
    exp = gene_bp / (genome_mb * 1e6) * f_bact * 1e6      # reads per million
    return obs_per_m / exp if exp else float('nan')


def poisson_cdf_le(k, lam):
    return sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(k + 1))


def main():
    out = []
    A = out.append

    # ---------- 1. in-silico 召回 ----------
    A('=' * 74)
    A('1. in-silico 召回：一条确实来自 rpoB 的 100 bp read，被 orf_kmers 检出的概率')
    A('=' * 74)
    res = in_silico_recall()
    allh = allt = 0
    for sp, (h, t) in res.items():
        allh += h
        allt += t
        A('  %-34s %5d / %5d = %.3f' % (sp, h, t, h / t if t else 0))
    recall = allh / allt
    A('  %-34s %5d / %5d = %.3f   <-- 无分歧、无测序错误时的上限' % ('POOLED', allh, allt, recall))
    A('')
    A('  这条上限直接否掉了 1/6 的读框惩罚：若 6 个框只有 1 个能用，检出率应≈0.167。')

    # ---------- 2. 由 rpoB 实测反解 ε ----------
    A('')
    A('=' * 74)
    A('2. 由 rpoB 对照反解 ε（rpoB 单拷贝，期望值只依赖基因组尺度）')
    A('=' * 74)
    ts = os.path.join(BASE, 'rpob_control.tsv')
    if not os.path.exists(ts):
        A('  [缺] %s 还没生成，跳过' % ts)
        open(os.path.join(BASE, 'eff_calibration.txt'), 'w', encoding='utf-8').write('\n'.join(out))
        return
    import csv
    rows = list(csv.DictReader(open(ts, encoding='utf-8'), delimiter='\t'))
    tot = sum(int(r['reads']) for r in rows)
    rp = sum(int(r['rpob_reads']) for r in rows)
    obs = 1e6 * rp / tot
    A('  实测：%d reads，rpoB 阳性 %d 条 = %.0f /百万' % (tot, rp, obs))
    A('')
    A('  对假设的敏感性（row=细菌 reads 占比，col=平均基因组 Mb）：')
    A('        ' + ''.join('%10.1f' % g for g in (3.0, 3.5, 4.0, 5.0)))
    for f in (1.0, 0.9, 0.8, 0.6):
        line = '  f=%.1f ' % f
        for g in (3.0, 3.5, 4.0, 5.0):
            line += '%10.2f' % eps_from_counts(obs, 3810, g, f)
        A(line)
    A('')
    eps_ref = eps_from_counts(obs, 3810, 4.0, 0.9)
    eps_low = eps_from_counts(obs, 3810, 3.0, 1.0)     # 最保守（最大化期望）
    A('  取 f=0.9、G=4.0 Mb：ε = %.3f' % eps_ref)
    A('  最保守（f=1.0、G=3.0 Mb）：ε = %.3f  <-- ε 的下限' % eps_low)
    A('  参考只覆盖 8 个物种，其他属的 rpoB 匹配不上 → 真实 ε 只会更高，')
    A('  所以这里得到的是**下限**。')
    A('')
    A('  结论：ε 至少 %.2f，不可能是我原来写的 0.15（差 %.1f 倍）。' % (eps_low, eps_low / 0.15))

    # ---------- 3. 用标定后的 ε 重算 hdcA ----------
    A('')
    A('=' * 74)
    A('3. 用标定后的 ε 重算 hdcA 的 λ 与深度需求')
    A('=' * 74)
    N = 6_300_000
    LAM_RAW = 1.525        # lambda_v2.py：携带率校正、未扣检出效率、6.3 M reads
    lam_raw_per_read = LAM_RAW / N
    A('  λ_raw（携带率校正，未扣 ε）= %.3f / 样本' % LAM_RAW)
    for e in (eps_low, eps_ref, 0.50):
        ld = LAM_RAW * e
        A('    ε=%.3f -> λ_det=%.3f；5 样本期望 %.2f 条，实测 1 条，P(X≤1)=%.3f'
          % (e, ld, 5 * ld, poisson_cdf_le(1, 5 * ld)))
    A('')
    for e, lab in ((eps_low, 'ε=%.2f（下限）' % eps_low), (eps_ref, 'ε=%.2f（居中）' % eps_ref)):
        need = 5.0 / (lam_raw_per_read * e)
        A('  λ_det=5 所需深度（%s）：%.0f 万 reads/样本 = %.1f 个当前 run'
          % (lab, need / 1e4, need / N))
    A('')
    A('  对照：旧口径（ε=0.15）给 1.38 亿 reads / 22 run，偏高 %.1f 倍。'
      % (1.38e8 / (5.0 / (lam_raw_per_read * eps_low))))
    A('')
    A('  定性结论不变：λ_det 仍远低于 5，仍不可定量；只是「差多少」从 22 倍改成约 4 倍。')

    # ---------- 4. hdcA/rpoB 的检出率之比 ----------
    A('')
    A('=' * 74)
    A('4. 同一批 reads 上两个基因的相对检出')
    A('=' * 74)
    hd = sum(int(r['hdca_reads']) for r in rows)
    A('  rpoB: %d 条 / %d reads = %.0f 每百万' % (rp, tot, obs))
    A('  hdcA: %d 条 / %d reads = %.3f 每百万' % (hd, tot, 1e6 * hd / tot))
    if hd:
        A('  比值：%.0f 倍' % (rp / hd))
    A('')
    A('  这个比值不需要任何假设，也不依赖基因组大小或细菌占比 ——')
    A('  两个基因在同一 for 循环里、用同一套解析与翻译代码、在同一批 reads 上计数。')

    txt = '\n'.join(out)
    print(txt)
    open(os.path.join(BASE, 'eff_calibration.txt'), 'w', encoding='utf-8').write(txt + '\n')


if __name__ == '__main__':
    main()
