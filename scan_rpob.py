# -*- coding: utf-8 -*-
"""rpoB 阳性对照：与 hdcA 在同一批 reads、同一条代码路径上并行扫描。

为什么需要这个对照
------------------
已有的 16S 引物对照只验证「核苷酸精确匹配」这条路径；hdcA 走的却是
「六框翻译 → 按终止子切 ORF → 15-mer 匹配」另一条路径。16S 命中再多，
也不能证明后者没坏，所以需要一条只依赖基因组尺度的对照。

rpoB（RNA 聚合酶 β 亚基）是细菌单拷贝核心基因：每个基因组恰好一份，
不受操纵子拷贝数、不依赖丰度表，因此是天然的量尺 —— 它把「检出能力」
从「我们相信流程没问题」变成「同一批 reads 里，流程实际找到了多少条」。

配对设计
--------
两个基因家族在同一 for 循环里计数，共享同一份 reads、同一套解析与翻译代码、
同一次运行。因此两组数字之间不存在批次、深度、样本构成的差异，
唯一差别就是基因本身。

用法
----
  python scan_rpob.py --limit 500000            # 先小规模试跑
  python scan_rpob.py reads/*.gz                # 全量，逐文件一行结果
"""
import argparse
import csv
import json
import os
import statistics as st
import sys
import time
from collections import Counter

import scan_hdca as SH

BASE = os.path.dirname(os.path.abspath(__file__))
K = 15
MIN_SPECIES_SHARE = 2       # 保守核心：至少在 2 个物种中出现


def load_rpob():
    p = os.path.join(BASE, 'rpob_ref.json')
    d = json.load(open(p, encoding='utf-8'))
    per = {}
    for sp, v in d.items():
        seq = v['seq']
        per[sp] = {seq[i:i + K] for i in range(len(seq) - K + 1)}
    return per


def build_index(per):
    """k-mer -> 携带该 k-mer 的物种集合"""
    idx = {}
    for sp, ks in per.items():
        for k in ks:
            idx.setdefault(k, set()).add(sp)
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fastq', nargs='+')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--out', default='rpob_control.tsv')
    args = ap.parse_args()

    # --- 参考：hdcA（复用既有库）与 rpoB（新抓的 8 个全长序列）---
    per_sp, type_of, conserved = SH.load_ref()
    hd_kmer2sp = {}
    for sp, ks in per_sp.items():
        for k in ks:
            hd_kmer2sp.setdefault(k, set()).add(sp)
    hd_cons = conserved['pyruvoyl'] | conserved['PLP']

    rp_per = load_rpob()
    rp_idx = build_index(rp_per)
    union = set(rp_idx)
    core = {k for k, v in rp_idx.items() if len(v) >= MIN_SPECIES_SHARE}

    print('hdcA ref : %d species, %d kmers' % (len(per_sp), len(hd_kmer2sp)), flush=True)
    print('rpoB ref : %d species, %d kmers union, %d kmers in >=%d species'
          % (len(rp_per), len(union), len(core), MIN_SPECIES_SHARE), flush=True)

    rows = []
    for path in args.fastq:
        sid = os.path.basename(path).split('_')[0].replace('.gz', '')
        t0 = time.time()
        n = 0
        pcr = 0
        hd_reads = 0
        rp_reads = 0
        rp_core_reads = 0
        rp_sp = Counter()
        rp_nk = []            # 每条阳性 read 命中的 rpoB k-mer 数
        hd_nk = []
        for seq in SH.fastq_records(path, args.limit or None):
            n += 1
            # 16S 引物对照
            rc = SH.revcomp(seq)
            for p in SH.P515_FULL:
                if p in seq or p in rc:
                    pcr += 1
                    break
            else:
                for p in SH.P806_FULL:
                    if p in seq or p in rc:
                        pcr += 1
                        break

            ks = SH.orf_kmers(seq)
            if not ks:
                continue

            # rpoB
            hit_r = ks & union
            if hit_r:
                rp_reads += 1
                rp_nk.append(len(hit_r))
                if hit_r & core:
                    rp_core_reads += 1
                touched = set()
                for k in hit_r:
                    touched |= rp_idx[k]
                for sp in touched:
                    rp_sp[sp] += 1

            # hdcA（同一份 ks，不重复翻译）
            touched_h = set()
            for k in ks:
                s = hd_kmer2sp.get(k)
                if s:
                    touched_h |= s
            if touched_h or (ks & hd_cons):
                hd_reads += 1
                hd_nk.append(len(touched_h) if touched_h else 1)

            if n % 500000 == 0:
                el = time.time() - t0
                print('  %-12s %9d reads  %.0fs  %.0f r/s  rpoB=%d  hdcA=%d 16S=%d'
                      % (sid, n, el, n / el, rp_reads, hd_reads, pcr), flush=True)

        el = time.time() - t0
        row = dict(
            sample=sid, reads=n, pcr_16S=pcr, pcr_per_M=1e6 * pcr / n if n else 0,
            rpob_reads=rp_reads, rpob_per_M=1e6 * rp_reads / n if n else 0,
            rpob_core_reads=rp_core_reads,
            rpob_kmer_median=st.median(rp_nk) if rp_nk else 0,
            rpob_kmer_max=max(rp_nk) if rp_nk else 0,
            hdca_reads=hd_reads,
            rpob_top=(';'.join('%s(%d)' % (a, b) for a, b in rp_sp.most_common(3))
                      if rp_sp else ''),
            seconds=round(el, 1),
        )
        rows.append(row)
        print('  DONE %-12s %d reads  rpoB %d (%.0f/M)  hdcA %d  16S %d (%.0f/M)  %.0fs'
              % (sid, n, rp_reads, row['rpob_per_M'], hd_reads, pcr,
                 row['pcr_per_M'], el), flush=True)

    out = os.path.join(BASE, args.out)
    with open(out, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print('\nwrote %s' % out)

    # --- 汇总与期望值 ---
    tot = sum(r['reads'] for r in rows)
    rp = sum(r['rpob_reads'] for r in rows)
    hd = sum(r['hdca_reads'] for r in rows)
    print('\n==== SUMMARY ====')
    print('total reads      : %d' % tot)
    print('rpoB positive    : %d  (%.1f per million)' % (rp, 1e6 * rp / tot))
    print('hdcA positive    : %d  (%.2f per million)' % (hd, 1e6 * hd / tot))
    if hd:
        print('ratio rpoB/hdcA  : %.0f x' % (rp / hd))
    print()
    print('NOTE: the rpoB/hdcA ratio above needs no assumptions and is the primary result.')
    print('      eps must NOT be derived from frame geometry (six-frame translation always')
    print('      contains the coding frame; in-silico recall = 1.000). Use eff_calibration.py,')
    print('      which back-solves eps from this table: 0.48 (lower bound) to 0.94 per sample.')


if __name__ == '__main__':
    main()
