# -*- coding: utf-8 -*-
"""
motu2matrix.py —— 把 mOTU 丰度表折叠成 species × sample 矩阵（供跨队列重算使用）。

输入（Wave 3 队列，Front Aging Neurosci 2022;14:820108 的 Supplementary Table 2）：
  - 工作表 `Our`：7,726 行 mOTU × 86 列样本；行标签为
    `k__Bacteria|p__...|g__X|s__Y` 形式的 mOTU 共识分类
  - Supplementary Table 1：样本分组与年龄（HL 90–103 / HE 65–75 / HY 24–48）

处理：
  1. 取 `|s__` 之后的物种名；同名 mOTU 相加（mOTU 是一种更细的分类单元，
     同一物种可能对应多个 mOTU，必须求和才等价于 species-level 丰度）
  2. 丢弃非细菌/古菌（k__Viruses、k__Eukaryota）与无法归种的条目
  3. 按列重新归一化到 1（原表列和约 0.98，去掉不可归种的部分后需要补回）
  4. 只保留指定分组的样本（默认 HL + HE = 年龄 ≥65）

运行环境：Python 3.8+，仅标准库。
"""
import csv
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
sys.path.insert(0, os.path.join(ROOT, '工具'))
import xlsx2csv  # noqa: E402

RAW = os.path.join(BASE, '_raw')
KEEP_GROUPS = ('HL', 'HE')


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return 0.0


def read_sheet(xlsx, sheet, out_csv):
    """用项目里的 xlsx2csv.load() 把工作表导出成 CSV。"""
    rows = xlsx2csv.load(xlsx, sheet)
    with open(out_csv, 'w', newline='', encoding='utf-8-sig') as f:
        csv.writer(f).writerows(rows)
    return out_csv


def build(motu_csv, meta_csv, out_csv, groups=KEEP_GROUPS, verbose=True):
    # --- 元数据 ---
    meta = {}
    with open(meta_csv, encoding='utf-8-sig') as f:
        rows = list(csv.reader(f))
    hdr_i = next(i for i, r in enumerate(rows) if r and r[0].strip() == 'Group')
    hdr = rows[hdr_i]
    ci_g, ci_s, ci_a = hdr.index('Group'), hdr.index('Sample_ID'), hdr.index('Age')
    for r in rows[hdr_i + 1:]:
        if not r or not r[0].strip():
            continue
        try:
            age = int(float(r[ci_a]))
        except (TypeError, ValueError):
            age = None
        meta[r[ci_s].strip()] = dict(group=r[ci_g].strip(), age=age)

    # --- mOTU 表 ---
    with open(motu_csv, encoding='utf-8-sig') as f:
        r = list(csv.reader(f))
    # 第 1 行：标题；第 2 行：Group；第 3 行：Sample；第 4 行：列名
    grp_row, smp_row, hdr_row = r[1], r[2], r[3]
    samples = [x.strip() for x in smp_row[2:]]
    sel = [j for j, s in enumerate(samples) if meta.get(s, {}).get('group') in groups]
    if verbose:
        from collections import Counter
        c = Counter(meta.get(s, {}).get('group') for s in samples)
        print('全部样本分组:', dict(c))
        print('保留 %s → %d 例' % ('+'.join(groups), len(sel)))

    tax_of = {}
    acc = {}
    for row in r[4:]:
        if len(row) < 3 or not row[0].strip():
            continue
        tax = row[1]
        if not tax.startswith('k__Bacteria') and not tax.startswith('k__Archaea'):
            continue
        if '|s__' not in tax:
            continue
        sp = tax.split('|s__')[-1].strip()
        if not sp:
            continue
        low = sp.lower()
        # 剔除未定种单元（unknown X / X sp. / gen.），口径与 MetaPhlAn 的已定名物种一致：
        # ELDERMET 矩阵只含 836 个已定名物种，未分类部分不在矩阵内。
        if low.startswith('unknown') or ' sp.' in low or 'gen.' in low or 'unclassified' in low:
            continue
        # 去掉 mOTU 的 clade 后缀，保留可归种的名字
        if sp.endswith(' [C]'):
            sp = sp[:-4].strip()
        key = sp
        vec = acc.get(key)
        if vec is None:
            vec = [0.0] * len(sel)
            acc[key] = vec
            tax_of[key] = tax
        for k, j in enumerate(sel):
            vec[k] += fnum(row[2 + j])

    # --- 列归一 ---
    n = len(sel)
    for k in range(n):
        tot = sum(v[k] for v in acc.values())
        if tot > 0:
            for v in acc.values():
                v[k] /= tot

    with open(out_csv, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f)
        w.writerow(['SampleID'] + sorted(acc))
        for k, j in enumerate(sel):
            w.writerow([samples[j]] + ['%.6g' % acc[s][k] for s in sorted(acc)])
    if verbose:
        print('species=%d  samples=%d  → %s' % (len(acc), n, out_csv))
    return dict(n_species=len(acc), n_samples=n,
                ids=[samples[j] for j in sel])


def main():
    raw_dir = RAW
    motu = os.path.join(raw_dir, 'table2_Our.csv')
    meta = os.path.join(raw_dir, 'table1_meta.csv')
    for p in (motu, meta):
        if not os.path.exists(p):
            print('缺少 %s' % p)
            print('先把 Supplementary Table 1/2 导出到 数据/_raw/：')
            print('  python 工具/xlsx2csv.py Table_2.xlsx "Our" 数据/_raw/table2_Our.csv')
            print('  python 工具/xlsx2csv.py Table_1.xlsx "Sheet1" 数据/_raw/table1_meta.csv')
            return 1
    build(motu, meta, os.path.join(BASE, 'Wave3_species_59x.csv'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
