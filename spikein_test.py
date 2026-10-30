# -*- coding: utf-8 -*-
"""
spike-in 对照：把合成的 hdcA reads 混进真实数据，验证流程能检出。
这是排除「流程系统性 bug 导致永不报阳性」的唯一办法。
"""
import os, random, gzip, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

BASE = os.path.dirname(os.path.abspath(__file__))
random.seed(11)

CODON = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'], 'K': ['AAA', 'AAG'],
    'L': ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
    'M': ['ATG'], 'N': ['AAT', 'AAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'W': ['TGG'], 'Y': ['TAT', 'TAC'],
}


def back_translate(pep):
    return ''.join(random.choice(CODON[a]) for a in pep if a in CODON)


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def get_seqs():
    """从 xlsx 取 E. lenta 与 S. thermophilus 的 hdcA 蛋白序列"""
    import zipfile, re
    from xml.etree import ElementTree as ET
    NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    ns = {'m': NS}
    z = zipfile.ZipFile(os.path.join(BASE, '12864_2021_8004_MOESM1_ESM.xlsx'))
    ss = []
    for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', ns):
        ss.append(''.join(t.text or '' for t in si.iter('{%s}t' % NS)))

    def ci(ref):
        c = re.match(r'([A-Z]+)', ref).group(1)
        n = 0
        for ch in c:
            n = n * 26 + (ord(ch) - 64)
        return n - 1
    R = {}
    for row in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).iter('{%s}row' % NS):
        v = {}
        for c in row.findall('m:c', ns):
            vv = c.find('m:v', ns)
            t = c.get('t')
            val = ss[int(vv.text)] if (t == 's' and vv is not None) else (vv.text if vv is not None else '')
            v[ci(c.get('r'))] = (val or '')
        R[int(row.get('r'))] = v
    want = {'Eggerthella lenta', 'Streptococcus thermophilus', 'Morganella morganii'}
    out = {}
    for r in sorted(R):
        if r < 3:
            continue
        v = R[r]
        sp = v.get(0, '').strip()
        if sp not in want:
            continue
        txt = v.get(6, '').strip()
        seq = ''.join(l.strip() for l in txt.split('\n')[1:]).replace(' ', '').upper()
        clean = []
        for ch in seq:
            if ch in 'ACDEFGHIKLMNPQRSTVWY':
                clean.append(ch)
            else:
                break
        out[sp] = ''.join(clean)
    return out


def make_spike_reads(pep, n=10, read_len=100):
    dna = back_translate(pep)
    reads = []
    for i in range(n):
        # 随机起点，随机链
        L = read_len
        if len(dna) <= L:
            s = dna
        else:
            p = random.randint(0, len(dna) - L)
            s = dna[p:p + L]
        if random.random() < 0.5:
            s = rc(s)
        reads.append(s)
    return reads


def main():
    src = os.path.join(BASE, 'reads', 'ERR3986546_1.part.fastq.gz')
    out = os.path.join(BASE, 'reads', 'spikein_test.fastq')
    seqs = get_seqs()
    print('取得参考序列:', {k: len(v) for k, v in seqs.items()})

    spiked = []
    for sp, pep in seqs.items():
        for s in make_spike_reads(pep, 10):
            spiked.append((sp, s))
    print('合成 reads 数:', len(spiked))

    # 写混合 fastq：真实 reads 10 万条 + 10 条 spike（pyr）
    # 注意：只对 pyr 型做 spike，PLP 的 Morganella 放在最后单独测
    op = gzip.open if src.endswith('.gz') else open
    f = op(src, 'rt', errors='replace')
    n = 0
    with open(out, 'w') as w:
        while n < 100000:
            try:
                h = f.readline()
                if not h:
                    break
                s = f.readline()
                q1 = f.readline()
                q2 = f.readline()
            except Exception:
                break
            w.write(h)
            w.write(s)
            w.write(q1)
            w.write(q2)
            n += 1
        # 在中间插入 spike
        for i, (sp, s) in enumerate(spiked):
            w.write('@spike_%d_%s\n' % (i, sp.replace(' ', '_')))
            w.write(s + '\n')
            w.write('+\n')
            w.write('I' * len(s) + '\n')
    try:
        f.close()
    except Exception:
        pass
    print('写出 %s : 真实 %d 条 + spike %d 条' % (out, n, len(spiked)))


if __name__ == '__main__':
    main()
