# -*- coding: utf-8 -*-
"""多家族面板扫描（多进程版）——对上千个 Pfam 家族在 12 物种 × 12 株上跑域搜索。

为什么要多进程
--------------
单进程版 run_multi_family.py 的顺序是「外层 proteome、内层家族」，每个 proteome
只加载一次序列，看起来省，但实际上 1,410 个 HMM 的 hmmsearch 占了 95% 时间，
16 核全闲着；实测 49 秒才起步。本版按 proteome 分片到 8 个进程，
每个进程在初始化时把 HMM 读一次（spawn 下每进程独立），再顺序处理自己分到的株。

口径与 run_multi_family.py 完全一致（已用 BSH 逐位复现原结果）：
    c_eff = Σ w_s c_s,  H = Σ w_s (c_s − c_eff)²,  H* 按 n_sp 对数插值,
    ratio = H/H*,  λ_raw = N·Σ(a_s c_s)·L_gene/G_s,  λ_det = 0.71·λ_raw

用法：python run_panel_mp.py --domains-file candidates_meta.txt --tag panel --procs 8
"""
import argparse
import csv
import math
import multiprocessing as mp
import os
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))
POS = os.path.join(ROOT, '阳性对照')
HMMDIR = os.path.join(BASE, 'hmm')
OUT = os.path.join(BASE, 'out')
MATRIX = os.path.join(ROOT, '数据', 'ELDERMET_species_189x836.csv')
HSTAR = {3: 0.0240, 5: 0.0378, 10: 0.0732, 20: 0.0936}
READS = 6.3e6
EPS = 0.71

_HMMS = None
_PFS = None


def _init(pfs):
    """worker 初始化：读入全部 HMM（每进程一份）。"""
    global _HMMS, _PFS
    import pyhmmer
    from pyhmmer.plan7 import HMMFile
    os.chdir(POS)                      # pyhmmer 只吃 ASCII 相对路径
    _PFS = pfs
    _HMMS = {}
    for pf in pfs:
        p = os.path.join('..', 'experiments-p2-p3', 'multi-family', 'hmm', pf + '.hmm.gz')
        if not os.path.exists(p):
            continue
        try:
            with HMMFile(p) as hf:
                _HMMS[pf] = hf.read()
        except Exception:
            pass


def _scan_one(job):
    """处理一个 proteome，返回 (species, upid, n_prot, total_aa, counts, medlens)。"""
    import pyhmmer
    from pyhmmer.easel import Alphabet, DigitalSequenceBlock, SequenceFile
    species, upid = job
    import os as _os
    rel = 'proteomes/%s.fasta.gz' % upid
    if not _os.path.exists(rel):
        return species, upid, 0, 0, {}, {}
    abc = Alphabet.amino()
    seqs, total_aa = [], 0
    with SequenceFile(rel, digital=True, alphabet=abc) as sf:
        for s in sf:
            seqs.append(s)
            total_aa += len(s)
    block = DigitalSequenceBlock(abc, seqs)
    counts, medlens = {}, {}
    for pf, hmm in _HMMS.items():
        n = 0
        L = []
        for hits in pyhmmer.hmmsearch([hmm], block, E=1e-5, cpus=1):
            for h in hits:
                if h.included:
                    n += 1
                    if h.length:
                        L.append(int(h.length))
        if n:
            counts[pf] = n
            medlens[pf] = sorted(L)[len(L) // 2] if L else 0
    return species, upid, len(seqs), total_aa, counts, medlens


def top_species(n):
    with open(MATRIX, encoding='utf-8-sig') as f:
        r = csv.reader(f)
        hdr = next(r)
        rows = list(r)
    total = len(rows)
    out = []
    for j, name in enumerate(hdr[1:], 1):
        v = [float(x[j]) / 100 for x in rows]
        out.append((sum(v) / total, name))
    out.sort(reverse=True)
    return out[:n]


def cached_proteomes(species_key, limit):
    path = os.path.join(POS, 'domain_hits.tsv')
    ids = []
    with open(path, encoding='utf-8') as f:
        r = csv.reader(f, delimiter='\t')
        next(r)
        for row in r:
            if row[0].replace('_', ' ') == species_key.replace('_', ' '):
                ids.append(row[1])
    seen, out = set(), []
    for i in ids:
        if i not in seen:
            seen.add(i)
            out.append(i)
    return out[:limit]


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--domains-file', required=True)
    ap.add_argument('--species', type=int, default=12)
    ap.add_argument('--per', type=int, default=12)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--tag', default='panel')
    a = ap.parse_args()

    pfs = [l.strip() for l in open(a.domains_file, encoding='utf-8')
           if l.strip() and not l.startswith('#')]
    pfs = [p for p in pfs if os.path.exists(os.path.join(HMMDIR, p + '.hmm.gz'))]
    print('家族 %d（有 HMM）' % len(pfs), flush=True)

    sp = top_species(a.species)
    jobs, ab = [], {}
    for mean_ab, name in sp:
        ids = cached_proteomes(name, a.per)
        ab[name] = mean_ab
        for u in ids:
            jobs.append((name, u))
    print('proteome 任务 %d，进程 %d' % (len(jobs), a.procs), flush=True)

    t0 = time.time()
    results = []
    with mp.Pool(a.procs, initializer=_init, initargs=(pfs,)) as pool:
        for i, r in enumerate(pool.imap_unordered(_scan_one, jobs, chunksize=1), 1):
            results.append(r)
            if i % 10 == 0:
                el = time.time() - t0
                print('  %d/%d  %.0fs  预计剩余 %.0fs' %
                      (i, len(jobs), el, el / i * (len(jobs) - i)), flush=True)

    # 汇总 c_s
    agg = {}
    for species, upid, nprot, aa, counts, medlens in results:
        if nprot == 0:
            continue
        d = agg.setdefault(species, {'n': 0, 'G': [], 'present': {}, 'lens': {}})
        d['n'] += 1
        d['G'].append(aa * 3)
        for pf, n in counts.items():
            d['present'][pf] = d['present'].get(pf, 0) + 1
        for pf, L in medlens.items():
            d['lens'].setdefault(pf, []).append(L)

    cs_rows = []
    for species, d in agg.items():
        G = sorted(d['G'])[len(d['G']) // 2]
        for pf, n in d['present'].items():
            L = sorted(d['lens'].get(pf, [0]))
            cs_rows.append([species.replace('_', ' '), pf, n, d['n'],
                            n / d['n'], ab[species], L[len(L) // 2], G])

    fam_rows = []
    byfam = {}
    for r in cs_rows:
        byfam.setdefault(r[1], []).append(r)
    for pf, rows in byfam.items():
        rows = [r for r in rows if r[4] > 0]
        if not rows:
            continue
        wsum = sum(r[5] for r in rows)
        ceff = sum(r[5] / wsum * r[4] for r in rows)
        H = sum(r[5] / wsum * (r[4] - ceff) ** 2 for r in rows)
        nsp = len(rows)
        hs = hstar_at(nsp)
        Ls = [int(r[6]) for r in rows if r[6]]
        L = (Ls[len(Ls) // 2] if Ls else 0) * 3
        lam = sum(READS * r[5] * r[4] * L / r[7] for r in rows if r[7])
        fam_rows.append([pf, nsp, '%.4f' % ceff, '%.4f' % H, '%.4f' % hs,
                         '%.3f' % (H / hs), '%.2f' % lam, '%.2f' % (lam * EPS),
                         L, '%.4f' % wsum])

    os.makedirs(OUT, exist_ok=True)
    p1 = os.path.join(OUT, 'carriage_%s.tsv' % a.tag)
    p2 = os.path.join(OUT, 'families_%s.tsv' % a.tag)
    with open(p1, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['species', 'pfam', 'n_present', 'n_strains', 'c_s',
                    'mean_abundance', 'median_hit_len', 'genome_bp'])
        w.writerows(cs_rows)
    with open(p2, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['pfam', 'n_sp', 'c_eff', 'H', 'H_star', 'H_over_Hstar',
                    'lambda_raw', 'lambda_det', 'L_gene', 'carrier_abundance'])
        w.writerows(fam_rows)

    nfail = sum(1 for r in fam_rows if float(r[5]) > 1)
    print('\n家族 %d，有携带者 %d；H/H*>1（fail）%d，pass %d'
          % (len(pfs), len(fam_rows), nfail, len(fam_rows) - nfail))
    print('用时 %.0f 分钟' % ((time.time() - t0) / 60))
    print('写出\n  %s\n  %s' % (p1, p2))


if __name__ == '__main__':
    main()
