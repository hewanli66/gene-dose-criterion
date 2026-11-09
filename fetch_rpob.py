# -*- coding: utf-8 -*-
"""抓取 rpoB（RNA 聚合酶 β 亚基）蛋白序列，作为蛋白级搜索路径的阳性对照。

为什么要这个对照
----------------
现有 16S 引物对照只验证了「核苷酸精确匹配」这一条代码路径；
而 hdcA 用的是「六框翻译 → 按终止子切 ORF → 15-mer 匹配」另一条路径。
16S 命中再多，也不能证明后者没坏。零命中首先会被怀疑成翻译/切框那一步的
bug，所以需要一条只依赖基因组尺度的对照。

rpoB 是细菌单拷贝核心基因，每个基因组恰好一份，因此它的期望命中数
可以直接由基因组尺度算出（不依赖任何物种丰度表），是天然的量尺。

数据来源
--------
- taxid：ENA taxonomy REST（www.ebi.ac.uk）
- 序列：EBI Proteins API 镜像（同源 UniProt，境内可达）

输出 rpob_ref.json：{species: {accession, length, seq}}
"""
import json, os, sys, time, urllib.parse, urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
UA = {'User-Agent': 'Mozilla/5.0 (research; rpob-control)'}
JSONH = dict(UA, Accept='application/json')

# 覆盖 ELDERMET 中主要门类的肠道常见菌，每门尽量 2 种以上，避免只用变形菌门。
# 名录给多了无妨：抓不到的自动跳过，只保留 ≥MIN_LEN 的全长条目。
SPECIES = [
    # Bacteroidota
    'Bacteroides thetaiotaomicron', 'Bacteroides uniformis', 'Bacteroides fragilis',
    'Bacteroides caccae', 'Phocaeicola vulgatus', 'Parabacteroides distasonis',
    'Alistipes putredinis', 'Prevotella copri', 'Segatella copri',
    # Bacillota
    'Faecalibacterium prausnitzii', 'Ruminococcus bromii', 'Agathobacter rectalis',
    'Eubacterium rectale', 'Roseburia intestinalis', 'Anaerostipes hadrus',
    'Clostridium perfringens', 'Clostridioides difficile', 'Lactobacillus fermentum',
    # Actinomycetota
    'Bifidobacterium longum', 'Bifidobacterium adolescentis', 'Eggerthella lenta',
    'Collinsella aerofaciens',
    # Pseudomonadota
    'Escherichia coli', 'Klebsiella pneumoniae',
]

# 105 aa 之类的条目是结构域片段，会把对照的灵敏度压低，只收全长
MIN_LEN = 900


def get(url, headers=JSONH, tries=3, timeout=30):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            last = e
            time.sleep(1.5 * (i + 1))
    raise last


def taxid(name):
    u = 'https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/' + urllib.parse.quote(name)
    d = json.loads(get(u))
    if not d:
        return None
    for e in d:
        if e.get('formalName') == 'true':
            return e['taxId']
    return d[0]['taxId']


def rpob(tid):
    """返回该 taxid 下最长的 rpoB 序列（片段条目会被 MIN_LEN 过滤掉）"""
    u = ('https://www.ebi.ac.uk/proteins/api/proteins?gene=rpoB&taxid=%s&size=100' % tid)
    d = json.loads(get(u))
    best = None
    for e in d:
        seq = e.get('sequence', {}).get('sequence', '')
        if not seq or len(seq) < MIN_LEN:
            continue
        if any(c not in 'ACDEFGHIKLMNPQRSTVWY' for c in seq):
            continue
        if best is None or len(seq) > len(best[1]):
            best = (e.get('accession', '?'), seq)
    return best


def main():
    out = {}
    for name in SPECIES:
        try:
            tid = taxid(name)
            if not tid:
                print('  [skip] no taxid: %s' % name)
                continue
            got = rpob(tid)
            if not got:
                print('  [skip] no rpoB: %s (taxid %s)' % (name, tid))
                continue
            acc, seq = got
            out[name] = dict(taxid=tid, accession=acc, length=len(seq), seq=seq)
            print('  %-34s taxid %-7s %-10s %4d aa' % (name, tid, acc, len(seq)))
        except Exception as e:
            print('  [err] %s: %s' % (name, e))
        time.sleep(0.4)

    p = os.path.join(BASE, 'rpob_ref.json')
    json.dump(out, open(p, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    L = [v['length'] for v in out.values()]
    print('\n%d/%d species -> %s' % (len(out), len(SPECIES), p))
    if L:
        print('length: min %d  median %d  max %d aa' % (min(L), sorted(L)[len(L) // 2], max(L)))


if __name__ == '__main__':
    main()
