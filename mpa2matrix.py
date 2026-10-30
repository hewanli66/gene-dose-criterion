# -*- coding: utf-8 -*-
"""
mpa2matrix.py —— 把 MetaPhlAn 输出（clade × 样本，百分数）转成 样本 × 物种 CSV。

约定与 ELDERMET / CRC9 / Villages 等既有矩阵一致：
  · 首列列名 SampleID，样本名为行
  · 物种名去掉 ``s__`` 前缀，保留下划线（如 Bacteroides_vulgatus）
  · 值统一归一到 1（原表是百分数时除以 100）

用法：
  python mpa2matrix.py <in.tsv> <out.csv> [--keep-ids ids.txt]
"""
import csv
import sys


def clade_to_species(path):
    """返回 (samples, species_names, matrix[sample_idx][sp_idx])。"""
    with open(path, encoding='utf-8', errors='replace') as f:
        rd = csv.reader(f, delimiter='\t')
        # cMD 合并表首行是 "#mpa_v30_CHOCOPhlAn_201901" 版本注释，跳过
        hdr = next(rd)
        while hdr and hdr[0].startswith('#'):
            hdr = next(rd)
        samples = [h.strip() for h in hdr[1:]]
        sp_rows = []
        for r in rd:
            if not r or not r[0].strip():
                continue
            last = r[0].split('|')[-1].strip()
            # 只取末端为 s__ 的行：MetaPhlAn 4 在物种行下面还有 t__ SGB 行，
            # 两者相加会把每个物种算两遍（列和 200%，踩过的坑）
            if not last.startswith('s__'):
                continue
            sp_rows.append((last[3:].strip(), r[1:]))
    names = [x[0] for x in sp_rows]
    mat = [[0.0] * len(names) for _ in samples]
    for j, (_, vals) in enumerate(sp_rows):
        for i in range(len(samples)):
            v = vals[i] if i < len(vals) else ''
            try:
                mat[i][j] = float(v)
            except (TypeError, ValueError):
                mat[i][j] = 0.0
    return samples, names, mat


def main():
    src, out = sys.argv[1], sys.argv[2]
    keep = None
    if '--keep-ids' in sys.argv:
        keep = {l.strip() for l in open(sys.argv[sys.argv.index('--keep-ids') + 1],
                                        encoding='utf-8') if l.strip()}
    samples, names, mat = clade_to_species(src)

    # 判断单位：物种层列和若在 95~105 之间说明是百分数
    colsum = [sum(row) for row in mat[:50]]
    scale = 0.01 if colsum and sum(colsum) / len(colsum) > 2 else 1.0

    idx = [i for i, s in enumerate(samples) if keep is None or s in keep]
    with open(out, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(['SampleID'] + names)
        for i in idx:
            w.writerow([samples[i]] + ['%.6g' % (v * scale) for v in mat[i]])
    print('%s -> %s : %d 样本 × %d 物种（scale=%.2f）'
          % (src, out, len(idx), len(names), scale))


if __name__ == '__main__':
    main()
