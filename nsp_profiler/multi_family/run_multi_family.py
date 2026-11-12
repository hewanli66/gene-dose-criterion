# -*- coding: utf-8 -*-
"""多家族扩展：对一批 Pfam 家族算 c_s、c_eff、H、H*、λ，检验判据能否复现。

与 ../../阳性对照/domain_carriage.py 的关系
-------------------------------------------
那个脚本负责「建目录」（解析 taxid、下载菌株蛋白组、跑域搜索、写 c_s 表）。
本脚本只做家族层面的重算，并复用它的 proteome 缓存（167 个 fasta，含 12 物种
× 12 株），因此新增一个家族的成本只是下载 HMM（几十 KB）+ 秒级域搜索。

路径约定
--------
pyhmmer 的 C 层按 ASCII 处理路径，含中文的绝对路径会抛 FileNotFoundError。
工作目录切到 `阳性对照/`（路径含中文），此后所有传给 pyhmmer 的路径都用相对路径：
    proteomes/UPxxxxxx.fasta.gz                      （已缓存）
    ../experiments-p2-p3/multi-family/hmm/PFxxxxx.hmm.gz  （ASCII 全路径）
结果文件走 Python 的 open()，用绝对路径无碍。

判据口径与正文一致
------------------
    c_eff = Σ w_s c_s          （w_s 按携带物种的丰度归一）
    H     = Σ w_s (c_s-c_eff)²
    H*    = §5.4 的 ELDERMET 门槛，按 n_sp 对数插值
    ratio = H / H*             ratio>1 → fail（物种代理不可靠）
    λ_raw = N · Σ (a_s c_s) · L_gene / G_s ,  λ_det = 0.71 · λ_raw

用法
----
    python run_multi_family.py --domains PF02275,PF00884
    python run_multi_family.py --domains-file candidates.txt --per 12
"""
import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))          # .../experiments-p2-p3/multi-family
ROOT = os.path.dirname(os.path.dirname(BASE))              # 项目根
POS = os.path.join(ROOT, '阳性对照')
HMMDIR = os.path.join(BASE, 'hmm')
OUT = os.path.join(BASE, 'out')
MATRIX = os.path.join(ROOT, '数据', 'ELDERMET_species_189x836.csv')

HSTAR = {3: 0.0240, 5: 0.0378, 10: 0.0732, 20: 0.0936}   # §5.4，ELDERMET
READS = 6.3e6
EPS = 0.71


def get_hmm(pf, retries=4):
    os.makedirs(HMMDIR, exist_ok=True)
    path = os.path.join(HMMDIR, pf + '.hmm.gz')
    if os.path.exists(path) and os.path.getsize(path) > 2000:
        return path
    url = 'https://www.ebi.ac.uk/interpro/api/entry/pfam/%s?annotation=hmm' % pf
    for k in range(retries):
        r = subprocess.run(['curl', '-sL', '-m', '120', '-o', path, url],
                           capture_output=True)
        if os.path.exists(path) and os.path.getsize(path) > 2000:
            return path
        time.sleep(1.5 * (k + 1))
    return None


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
    """从 domain_hits.tsv 取该物种已缓存的 proteome 列表（按原 index 排序）。"""
    path = os.path.join(POS, 'domain_hits.tsv')
    if not os.path.exists(path):
        return []
    ids = []
    with open(path, encoding='utf-8') as f:
        r = csv.reader(f, delimiter='\t')
        next(r)
        for row in r:
            if row[0].replace('_', ' ') == species_key.replace('_', ' '):
                ids.append(row[1])
    # 去重保序
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
    ap.add_argument('--domains', default='PF02275')
    ap.add_argument('--domains-file', default=None)
    ap.add_argument('--species', type=int, default=12)
    ap.add_argument('--per', type=int, default=12)
    ap.add_argument('--tag', default='run')
    ap.add_argument('--quiet', action='store_true', help='不逐株打印（大家族面板用）')
    a = ap.parse_args()

    pfs = [x.strip() for x in a.domains.split(',') if x.strip()]
    if a.domains_file:
        pfs += [l.strip() for l in open(a.domains_file, encoding='utf-8')
                if l.strip() and not l.startswith('#')]
    pfs = list(dict.fromkeys(pfs))

    import pyhmmer
    from pyhmmer.easel import Alphabet, DigitalSequenceBlock, SequenceFile
    from pyhmmer.plan7 import HMMFile

    os.chdir(POS)                      # pyhmmer 只吃 ASCII 相对路径
    hmms = {}
    for pf in pfs:
        p = get_hmm(pf)
        if not p:
            print('HMM 下载失败 %s' % pf, flush=True)
            continue
        rel = os.path.relpath(p, POS).replace('\\', '/')
        with HMMFile(rel) as hf:
            hmms[pf] = hf.read()
        if not a.quiet:
            print('HMM %-9s %-8s M=%d' % (pf, hmms[pf].name, hmms[pf].M), flush=True)

    sp = top_species(a.species)
    cs_rows = []
    for mean_ab, name in sp:
        ids = cached_proteomes(name, a.per)
        if not ids:
            print('SKIP %-40s 无缓存蛋白组' % name, flush=True)
            continue
        present = {pf: 0 for pf in hmms}
        lens = {pf: [] for pf in hmms}
        genome_bp, n_ok = [], 0
        for upid in ids:
            rel = 'proteomes/%s.fasta.gz' % upid
            if not os.path.exists(rel):
                continue
            abc = Alphabet.amino()
            seqs, total_aa = [], 0
            with SequenceFile(rel, digital=True, alphabet=abc) as sf:
                for s in sf:
                    seqs.append(s)
                    total_aa += len(s)
            block = DigitalSequenceBlock(abc, seqs)
            n_ok += 1
            genome_bp.append(total_aa * 3)
            for pf, hmm in hmms.items():
                n = 0
                L = []
                for hits in pyhmmer.hmmsearch([hmm], block, E=1e-5, cpus=2):
                    for h in hits:
                        if h.included:
                            n += 1
                            if h.length:
                                L.append(int(h.length))
                if n:
                    present[pf] += 1
                    if L:
                        lens[pf].append(sorted(L)[len(L) // 2])
        if not n_ok:
            continue
        G = sorted(genome_bp)[len(genome_bp) // 2]
        for pf in hmms:
            Lm = sorted(lens[pf])[len(lens[pf]) // 2] if lens[pf] else 0
            cs_rows.append([name.replace('_', ' '), pf, present[pf], n_ok,
                            present[pf] / n_ok, mean_ab, Lm, G])
        if not a.quiet:
            print('  %-38s n_ok=%-3d %s' % (name, n_ok,
                  ' '.join('%s=%d/%d' % (pf, present[pf], n_ok) for pf in hmms)), flush=True)

    fam = []
    for pf in hmms:
        rows = [r for r in cs_rows if r[1] == pf and r[4] > 0]
        if not rows:
            fam.append([pf, '', 0, '', '', '', '', '', '', '', 'no_carrier'])
            continue
        wsum = sum(r[5] for r in rows)
        ceff = sum(r[5] / wsum * r[4] for r in rows)
        H = sum(r[5] / wsum * (r[4] - ceff) ** 2 for r in rows)
        nsp = len(rows)
        hs = hstar_at(nsp)
        L = (sorted(int(r[6]) for r in rows if r[6]) or [0])
        L = L[len(L) // 2] * 3
        lam = 0.0
        for r in rows:
            if r[7]:
                lam += READS * r[5] * r[4] * L / r[7]
        fam.append([pf, hmms[pf].name, nsp, '%.4f' % ceff, '%.4f' % H,
                    '%.4f' % hs, '%.2f' % (H / hs), '%.2f' % lam,
                    '%.2f' % (lam * EPS), L, ''])
    os.makedirs(OUT, exist_ok=True)
    cs_path = os.path.join(OUT, 'carriage_%s.tsv' % a.tag)
    fm_path = os.path.join(OUT, 'families_%s.tsv' % a.tag)
    with open(cs_path, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['species', 'pfam', 'n_present', 'n_strains', 'c_s',
                    'mean_abundance', 'median_hit_len', 'genome_bp'])
        w.writerows(cs_rows)
    with open(fm_path, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['pfam', 'hmm_name', 'n_sp', 'c_eff', 'H', 'H_star', 'H_over_Hstar',
                    'lambda_raw', 'lambda_det', 'L_gene', 'note'])
        w.writerows(fam)

    print('\n%-9s %-26s %3s %8s %8s %6s %9s  %s' %
          ('pfam', 'name', 'nsp', 'c_eff', 'H', 'ratio', 'lam_det', 'verdict'))
    for r in fam:
        if r[2] == 0:
            print('%-9s %-26s   -        -        -     -         -  no_carrier' %
                  (r[0], str(r[1])[:26]))
            continue
        v = 'FAIL' if float(r[6]) > 1 else 'pass'
        print('%-9s %-26s %3d %8s %8s %6s %9s  %s' %
              (r[0], str(r[1])[:26], r[2], r[3], r[4], r[6], r[8], v))
    print('\n写出 %s\n     %s' % (cs_path, fm_path))


if __name__ == '__main__':
    main()
