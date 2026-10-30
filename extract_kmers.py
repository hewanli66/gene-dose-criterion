# -*- coding: utf-8 -*-
"""
从 Table S1 的 96 条 HdcA 蛋白序列中提取类别特异的保守肽段(k-mer)。
用途：在没有 HMMER 的 Windows 环境下，用「六框翻译 + 保守肽段精确匹配」
      代替 hmmsearch，回答「PLP 型 hdcA 在肠道宏基因组里到底存不存在」。

设计要点：
 1. 分两类提取（pyruvoyl / PLP），互不混用，避免交叉命中。
 2. 只保留在同类 >= MIN_SUPPORT 条序列中出现的 k-mer —— 种间变异大，
    过严会漏，过松会假阳性。
 3. 输出两类 k-mer 的「交叉重叠」检查，确认两套互不干扰。
"""
import zipfile, re, os, json, itertools
from collections import Counter, defaultdict
from xml.etree import ElementTree as ET

NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
ns = {'m': NS}
BASE = os.path.dirname(os.path.abspath(__file__))
XL = os.path.join(BASE, '12864_2021_8004_MOESM1_ESM.xlsx')

K = 15          # 肽段长度(氨基酸)
MIN_SUPPORT_FRAC = 0.25   # 至少在 25% 同类序列中出现
STD_AA = set('ACDEFGHIKLMNPQRSTVWY')


def colidx(ref):
    c = re.match(r'([A-Z]+)', ref).group(1)
    n = 0
    for ch in c:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def load_seqs():
    z = zipfile.ZipFile(XL)
    ss = []
    for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', ns):
        ss.append(''.join(t.text or '' for t in si.iter('{%s}t' % NS)))
    R = {}
    for row in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).iter('{%s}row' % NS):
        v = {}
        for c in row.findall('m:c', ns):
            vv = c.find('m:v', ns)
            t = c.get('t')
            val = ss[int(vv.text)] if (t == 's' and vv is not None) else (vv.text if vv is not None else '')
            v[colidx(c.get('r'))] = (val or '')
        R[int(row.get('r'))] = v
    out = []
    for r in sorted(R):
        if r < 3:
            continue
        v = R[r]
        sp = v.get(0, '').strip()
        if not sp:
            continue
        txt = v.get(6, '').strip()
        if not txt.startswith('>'):
            continue
        ls = txt.split('\n')
        seq = ''.join(l.strip() for l in ls[1:]).replace(' ', '').upper()
        # 清理非标准字符：保留到第一个非标准字符之前
        clean = []
        for ch in seq:
            if ch in STD_AA:
                clean.append(ch)
            else:
                break
        out.append({'species': sp, 'type': v.get(3, '').strip(),
                    'raw_len': len(seq), 'seq': ''.join(clean)})
    return out


def kmers(seq, k=K):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def main():
    seqs = load_seqs()
    print('loaded %d sequences' % len(seqs))

    groups = {'pyruvoyl': [], 'PLP': []}
    for s in seqs:
        if s['type'] in groups:
            groups[s['type']].append(s)
        elif s['type'] == 'pyruvoyl/PLP':
            # 双型物种不入任一组，避免污染
            print('  skipped dual-type:', s['species'])

    results = {}
    for gname, glist in groups.items():
        n = len(glist)
        need = max(2, int(round(n * MIN_SUPPORT_FRAC)))
        cnt = Counter()
        for s in glist:
            for km in kmers(s['seq']):
                cnt[km] += 1
        # 只保留达到支持度的
        keep = {km: c for km, c in cnt.items() if c >= need}
        # 去掉被更短集合包含的冗余：保留出现次数高的
        results[gname] = keep
        print(f"\n{gname} : {n} sequences, {len(cnt)} distinct {K}-mers, "
              f"{len(keep)} kept (support>={need}, i.e. {100.0*MIN_SUPPORT_FRAC:.0f}%)")

        # 该组内覆盖度：有多少序列能被至少一个保留 k-mer 命中
        covered = 0
        for s in glist:
            if any(km in keep for km in kmers(s['seq'])):
                covered += 1
        print('   sequences covered by kept k-mers: %d/%d (%.0f%%)' % (covered, n, 100.0 * covered / n))

        # 展示最有代表性的几个（支持度最高的）
        top = sorted(keep.items(), key=lambda x: -x[1])[:5]
        for km, c in top:
            print('   %s  support=%d/%d' % (km, c, n))

    # 交叉检查
    a = set(results['pyruvoyl'])
    b = set(results['PLP'])
    inter = a & b
    print('\ncross-overlap between the two k-mer sets: %d' % len(inter))
    if inter:
        print('   examples:', list(inter)[:8])

    # 保存
    out = {
        'k': K,
        'pyruvoyl': sorted(results['pyruvoyl'].items(), key=lambda x: -x[1]),
        'PLP': sorted(results['PLP'].items(), key=lambda x: -x[1]),
    }
    p = os.path.join(BASE, 'conserved_kmers.json')
    with open(p, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1)
    print('\nwritten', p)

    # 同时导出人类可读 TSV
    p2 = os.path.join(BASE, 'conserved_kmers.tsv')
    with open(p2, 'w', encoding='utf-8-sig', newline='') as f:
        f.write('type\tkmer\tsupport\n')
        for t in ['pyruvoyl', 'PLP']:
            for km, c in out[t]:
                f.write('%s\t%s\t%d\n' % (t, km, c))
    print('written', p2)


if __name__ == '__main__':
    main()
