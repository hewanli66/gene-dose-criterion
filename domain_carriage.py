# -*- coding: utf-8 -*-
"""用「Pfam 域 + 完整蛋白组 + pyhmmer」测株水平携带率 c_s。

为什么这样做：UniProt 的 Pfam 交叉引用只覆盖已审校条目（全库仅数十条），不足以
代表菌株；NCBI 从分析机不可达；KEGG 对肠道厌氧菌几乎没有株深度。剩下唯一自洽的
做法是下载每个菌株的完整蛋白组，自己跑域检索判定有无。

输出（都在本目录）：
    domain_hits.tsv     每个（物种，蛋白组）的逐域命中数
    domain_carriage.tsv 每个（物种，域）的 c_s、株数
    domain_family.tsv   每个域的 c_eff、H、H*、λ、判决

用法：
    python domain_carriage.py --domains PF02275,PF02286,PF05985 --species 18 --per 16
"""
import argparse
import csv
import json
import math
import os
import random
import subprocess
import time
import urllib.parse

import pyhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, SequenceFile
from pyhmmer.plan7 import HMMFile

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, 'proteomes')
HMMDIR = os.path.join(BASE, 'hmm')
MATRIX = os.path.join(BASE, '..', '数据', 'ELDERMET_species_189x836.csv')

HSTAR = {3: 0.0240, 5: 0.0378, 10: 0.0732, 20: 0.0936}   # §5.4，ELDERMET
READS = 6.3e6


def curl(url, out=None, timeout=240, retries=4):
    for k in range(retries):
        cmd = ['curl', '-sL', '-m', str(timeout)]
        if out:
            cmd += ['-o', out]
        cmd.append(url)
        r = subprocess.run(cmd, capture_output=True)
        # 阈值 10 字节：UniProt 对无条目的蛋白组返回 20 字节空 gzip，那是有效响应而非失败
        ok = r.returncode == 0 and (out is None or
                                    (os.path.exists(out) and os.path.getsize(out) > 10))
        if ok:
            return r.stdout
        time.sleep(1.5 * (k + 1))
    return None


def get_hmm(pf):
    os.makedirs('hmm', exist_ok=True)
    path = os.path.join('hmm', pf + '.hmm.gz')   # 相对路径：pyhmmer 不接受含中文的绝对路径
    if not os.path.exists(path) or os.path.getsize(path) < 2000:
        curl('https://www.ebi.ac.uk/interpro/api/entry/pfam/%s?annotation=hmm' % pf, path)
    return path


def top_species(n):
    f = open(MATRIX, encoding='utf-8-sig')
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


def resolve_taxid(name):
    """学名（含过时的旧名）-> species 级 taxid。

    必要性：UniProt 已按新分类改名（Prevotella copri -> Segatella copri,
    Bacteroides dorei -> Phocaeicola dorei），ELDERMET 矩阵用的是旧名，
    直接按 organism_name 查 proteomes 会返回 0 条而被误判为「无蛋白组」。
    taxonomy 接口的 otherNames 保留旧名，故用自由文本检索再筛 rank。
    """
    raw = curl('https://rest.uniprot.org/taxonomy/search?query='
               + urllib.parse.quote(name) + '&format=json&size=10')
    if not raw:
        return None
    try:
        rs = json.loads(raw).get('results', [])
    except ValueError:
        return None
    for r in rs:                       # 先要完全同名的 species
        if r.get('rank') == 'species' and \
           (r.get('scientificName') or '').lower() == name.lower():
            return r.get('taxonId')
    for r in rs:                       # 退而求其次：任一 species
        if r.get('rank') == 'species':
            return r.get('taxonId')
    return None


def proteomes(taxname, limit, seed=20261030):
    """返回候选蛋白组队列：reference 优先，其余随机洗牌。"""
    q = ' '.join(taxname.split('_')[:2]).replace('_unclassified', '')
    if q.lower().endswith('unclassified') or q.lower().endswith('sp.'):
        q = q.rsplit(' ', 1)[0]        # 属级未定种：只用属名检索
    url = None
    taxid = resolve_taxid(q)
    if taxid:
        url = ('https://rest.uniprot.org/proteomes/search?query=organism_id:%d'
               '&format=json&size=500' % taxid)
    raw = curl(url) if url else None
    if not raw:
        url = ('https://rest.uniprot.org/proteomes/search?query='
               + urllib.parse.quote('organism_name:"%s"' % q) + '&format=json&size=500')
        raw = curl(url)
    if not raw:
        return q, []
    try:
        d = json.loads(raw)
    except ValueError:
        return q, []
    ref, other = [], []
    for r in d.get('results', []):
        if r.get('proteomeType') in ('Redundant proteome', 'Excluded proteome'):
            continue
        (ref if 'Reference' in (r.get('proteomeType') or '') else other).append(
            (r['id'], r.get('proteinCount') or 0))
    rnd = random.Random(seed)
    rnd.shuffle(other)
    ref.sort()
    return q, ref + other


def fetch_fasta(upid):
    """下载蛋白组 FASTA，返回**相对路径**（相对 BASE）。

    pyhmmer 的 C 层按 ASCII 处理路径，含中文的绝对路径会抛 FileNotFoundError，
    因此下载用绝对路径、返回用相对路径（调用方已 chdir 到 BASE）。
    """
    os.makedirs(CACHE, exist_ok=True)
    abspath = os.path.join(CACHE, upid + '.fasta.gz')
    if os.path.exists(abspath) and os.path.getsize(abspath) > 1000:
        return 'proteomes/' + upid + '.fasta.gz'
    if os.path.exists(abspath):
        os.remove(abspath)
    # UniParc 而非 UniProtKB：后者只导入了每个物种的参考株（organism_id 查询证实，
    # B. uniformis 全库仅 4,303 条 ≈ 单株蛋白数），其余菌株只有元数据、没有序列。
    # UniParc 持有全部菌株的序列并按 proteome 建了索引，返回条数与 proteinCount 吻合。
    # 记录在「策略」一节，勿改回 uniprotkb。
    curl('https://rest.uniprot.org/uniparc/stream?query=proteome:%s'
         '&format=fasta&compressed=true' % upid, abspath)
    if os.path.exists(abspath) and os.path.getsize(abspath) > 1000:
        return 'proteomes/' + upid + '.fasta.gz'
    return None


def scan(path, hmms):
    """返回 ({pf: 命中蛋白数}, {pf: 命中蛋白的中位长度}, 蛋白数, 氨基酸总量)

    注意 pyhmmer 的 Hit 暴露的是 .length（目标蛋白长度），不是 .target_length；
    写错字段会静默得到 0，进而使 λ 恒为 0——这个错误在单样本试跑时暴露过。
    """
    abc = Alphabet.amino()
    seqs = []
    total_aa = 0
    with SequenceFile(path, digital=True, alphabet=abc) as sf:
        for s in sf:
            seqs.append(s)
            total_aa += len(s)
    block = DigitalSequenceBlock(abc, seqs)
    counts, lens = {}, {}
    for pf, hmm in hmms.items():
        n = 0
        L = []
        for hits in pyhmmer.hmmsearch([hmm], block, E=1e-5, cpus=2):
            for h in hits:
                if h.included:
                    n += 1
                    tl = int(h.length or 0)
                    if tl:
                        L.append(tl)
        counts[pf] = n
        lens[pf] = sorted(L)[len(L) // 2] if L else 0
    return counts, lens, len(seqs), total_aa


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
    ap.add_argument('--domains', default='PF02275,PF02286,PF05985')
    ap.add_argument('--species', type=int, default=18)
    ap.add_argument('--per', type=int, default=16)
    a = ap.parse_args()

    os.chdir(BASE)   # pyhmmer 的 C 层不处理非 ASCII 绝对路径，统一用相对路径
    pfs = [x.strip() for x in a.domains.split(',') if x.strip()]
    hmms = {}
    for pf in pfs:
        with HMMFile(get_hmm(pf)) as hf:
            hmms[pf] = hf.read()
        print('HMM %s  %s  M=%d' % (pf, hmms[pf].name, hmms[pf].M), flush=True)

    species = top_species(a.species)
    hit_rows, cs_rows = [], []
    for mean_ab, name in species:
        q, ids = proteomes(name, a.per)
        if not ids:
            print('SKIP %-40s 无蛋白组' % q, flush=True)
            continue
        n_ok = 0
        genome_bp = []
        present = {pf: 0 for pf in pfs}
        medlen = {pf: [] for pf in pfs}
        for upid, np_ in ids:
            if n_ok >= a.per:
                break
            fa = fetch_fasta(upid)
            if not fa:      # 该蛋白组在 UniProtKB 里没有条目（空 gzip），换下一个
                continue
            try:
                counts, lens, nprot, total_aa = scan(fa, hmms)
            except Exception as e:
                print('  ERR', upid, e, flush=True)
                continue
            n_ok += 1
            genome_bp.append(total_aa * 3)
            for pf in pfs:
                if counts[pf]:
                    present[pf] += 1
                if lens[pf]:
                    medlen[pf].append(lens[pf])
            hit_rows.append([q, upid, n_ok, nprot] + [counts[pf] for pf in pfs])
            print('  %-40s %-12s n=%-5d %s' % (q, upid, nprot,
                  ' '.join('%s=%d' % (pf, counts[pf]) for pf in pfs)), flush=True)
        if n_ok == 0:
            continue
        G_med = sorted(genome_bp)[len(genome_bp) // 2] if genome_bp else 0
        for pf in pfs:
            L = sorted(medlen[pf])[len(medlen[pf]) // 2] if medlen[pf] else 0
            cs_rows.append([q, pf, present[pf], n_ok,
                            '%.4f' % (present[pf] / n_ok), '%.6f' % mean_ab, L, G_med])

    with open(os.path.join(BASE, 'domain_hits.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['species', 'proteome', 'index', 'n_proteins'] + pfs)
        w.writerows(hit_rows)
    with open(os.path.join(BASE, 'domain_carriage.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['species', 'domain', 'n_present', 'n_strains', 'c_s', 'mean_abundance',
                    'median_protein_len', 'genome_bp'])
        w.writerows(cs_rows)

    # 家族层面：carriers = c_s>0 的物种；w 按丰度归一
    fam_rows = []
    for pf in pfs:
        rows = [r for r in cs_rows if r[1] == pf and float(r[4]) > 0]
        if not rows:
            continue
        wsum = sum(float(r[5]) for r in rows)
        ceff = sum(float(r[5]) / wsum * float(r[4]) for r in rows)
        H = sum(float(r[5]) / wsum * (float(r[4]) - ceff) ** 2 for r in rows)
        nsp = len(rows)
        hs = hstar_at(nsp)
        L = sorted(int(r[6]) for r in rows)[len(rows) // 2] * 3
        # λ = N · Σ(a_s c_s) · L_gene / G_s，G_s 由蛋白数×1kb 估
        lam = 0.0
        for r in rows:
            G = int(r[7])
            if not G:
                continue
            lam += READS * float(r[5]) * float(r[4]) * L / G
        fam_rows.append([pf, nsp, '%.4f' % ceff, '%.4f' % H, '%.4f' % hs,
                         '%.2f' % (H / hs), '%.2f' % lam, L, wsum])
        print('\n%s  n_sp=%d  c_eff=%.3f  H=%.4f  H*=%.4f  H/H*=%.2f  λ=%.1f  L=%d'
              % (pf, nsp, ceff, H, hs, H / hs, lam, L), flush=True)

    with open(os.path.join(BASE, 'domain_family.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['domain', 'n_species', 'c_eff', 'H', 'H_star', 'ratio', 'lambda',
                    'L_gene', 'carrier_abundance'])
        w.writerows(fam_rows)


if __name__ == '__main__':
    main()
