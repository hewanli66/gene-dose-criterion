# -*- coding: utf-8 -*-
"""motu_chd2matrix.py —— 把 Feng 2016 冠心病队列的 mOTU prof 表转成 species × sample 矩阵。

输入 zen/motu_chd.txt（Zenodo 4991531 / mOTU.prof.upload.txt）：
  行为 mOTU，列为样本，首行是样本名，首列是行名。
  行名有两种：
    * 已定名物种（Bacteroides_uniformis）—— 保留
    * 未定名 motu_linkage_group_510 —— 丢弃（无参考基因组，无物种名，
      与 Zhu2020 / Wave3 两个 mOTU 矩阵的口径一致：那两个矩阵也只含已定名物种）
  丢弃后按列重新归一化到 1。

输出 CHD_species_102xNNN.csv（SampleID × 物种，与其它队列矩阵同格式）。
"""
import csv
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, 'zen', 'motu_chd.txt')
OUT = os.path.join(BASE, 'CHD_species_102x.csv')


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return 0.0


def norm(name):
    """mOTU 行名 → 物种名：下划线换空格，去掉 |[C] clade 后缀。"""
    s = name.strip()
    if s.endswith('[C]'):
        s = s[:-3].strip()
    s = s.replace('|', ' ').replace('_', ' ')
    return ' '.join(s.split())


def main():
    with open(SRC, encoding='utf-8', errors='replace') as f:
        rd = csv.reader(f, delimiter='\t')
        hdr = next(rd)
        samples = [s.strip() for s in hdr[1:] if s.strip()]
        acc = {}
        n_lg = 0
        lg_abund = 0.0
        total = 0.0
        for row in rd:
            if len(row) < 2 or not row[0].strip():
                continue
            name = norm(row[0])
            vals = [fnum(x) for x in row[1:1 + len(samples)]]
            vals += [0.0] * (len(samples) - len(vals))
            total += sum(vals)
            if not name or name.lower().startswith('motu linkage group'):
                n_lg += 1
                lg_abund += sum(vals)
                continue
            if name.lower().startswith('unassigned') or name.lower() in ('unknown', 'other'):
                n_lg += 1
                lg_abund += sum(vals)
                continue
            vec = acc.get(name)
            if vec is None:
                vec = [0.0] * len(samples)
                acc[name] = vec
            for k, v in enumerate(vals):
                vec[k] += v

    # 列归一
    for k in range(len(samples)):
        tot = sum(v[k] for v in acc.values())
        if tot > 0:
            for v in acc.values():
                v[k] /= tot

    keys = sorted(acc)
    with open(OUT, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f)
        w.writerow(['SampleID'] + keys)
        for k, s in enumerate(samples):
            w.writerow([s] + ['%.6g' % acc[key][k] for key in keys])

    print('samples=%d  species=%d  dropped LG rows=%d (%.1f%% of abundance)'
          % (len(samples), len(keys), n_lg, 100 * lg_abund / total if total else 0))
    print('→ %s' % OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
