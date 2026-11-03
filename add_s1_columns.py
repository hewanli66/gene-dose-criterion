# -*- coding: utf-8 -*-
"""给 Table S1（Additional file 7）补两列：contributing_species 与 ratio_nsp。

M1-b：ratio 列是按 matched 口径（目录匹配物种数）算的，因此同一队列的
n_sp = 3/5/10/20 四行完全相同。读者会误读为"比值随家族规模变化"。
本脚本补上 contributing（队列内平均丰度 > 0 的匹配物种数）与该比值的
实际 n_sp 口径（ratio_nsp），并把两列写进 TSV。幂等。

运行：python 工具/add_s1_columns.py
"""
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, '数据')
PACK = os.path.join(ROOT, '投稿包')
TSV = os.path.join(PACK, 'TableS1_threshold_detail.tsv')

# 队列 → (矩阵文件, 是否样本为行, 元数据过滤)
MATRICES = {
    'ELDERMET (Ireland)':        ('ELDERMET_species_189x836.csv', None),
    'CRC collection, 9 cohorts': ('CRC9_species_1262x917.csv', None),
    'Guangzhou (China)':         ('HeQ2017_samples_122x.csv', None),
    'SG90 (Singapore), all':     ('SG90_species_abundances.csv', None),
    'SG90, aged ≥77':            ('SG90_species_abundances.csv', ('SG90_metadata.csv', 77)),
    'Villages (Honduras)':       ('Village_species_1787x1910.csv', None),
    'Screening cohort (Japan)':  ('Yachida2019_samples_615x1835.csv', None),
    'Wave 3 (China), aged ≥65':  ('Wave3_species_59x.csv', None),
    "Xi'an cohort (China)":      ('Zhu2020_species_171x359.csv', None),
    'Coronary cohort (China)':   ('CHD_species_102x.csv', None),
}


def norm(s):
    return s.strip().lower().replace('_', ' ').replace('-', ' ')


def load_catalogue():
    out = {}
    with open(os.path.join(DATA, 'TableS1_species_fixed.csv'), encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            fr = r['strains_fraction_as_published'].strip()
            if '/' in fr:
                a, b = fr.split('/')
                v = float(a) / float(b)
            else:
                v = float(fr)
            out[norm(r['species'])] = v
    return out


def load_means(fn, filt):
    rows = list(csv.reader(open(os.path.join(DATA, fn), encoding='utf-8-sig')))
    hdr, body = rows[0], [r for r in rows[1:] if r]
    rows_are_samples = norm(hdr[0]) not in ('species', 'taxon', 'clade')
    keep = None
    if filt:
        meta_fn, age_min = filt
        age = {}
        with open(os.path.join(DATA, meta_fn), encoding='utf-8-sig') as f:
            for m in csv.DictReader(f):
                try:
                    age[m['Subject ID'].strip()] = float(m['Age (in years)'])
                except (ValueError, KeyError):
                    pass
        keep = {k for k, v in age.items() if v >= age_min}

    if not rows_are_samples:
        # 行=物种，列=样本
        sums, n = {}, 0
        idx = list(range(1, len(hdr)))
        for r in body:
            vals = []
            for i in idx:
                if i < len(r) and r[i] not in ('', 'NA'):
                    try:
                        vals.append(float(r[i]))
                    except ValueError:
                        pass
            if not vals:
                continue
            n = max(n, len(vals))
            sums[norm(r[0])] = sum(vals)
        return {k: v / n for k, v in sums.items()} if n else {}

    # 行=样本，列=物种
    sp = [norm(c) for c in hdr[1:]]
    acc = [0.0] * len(sp)
    n = 0
    for r in body:
        sid = r[0].strip()
        if keep is not None and sid not in keep:
            continue
        n += 1
        for j, cell in enumerate(r[1:]):
            if cell not in ('', 'NA'):
                try:
                    acc[j] += float(cell)
                except ValueError:
                    pass
    if n == 0:
        return {}
    return {s: a / n for s, a in zip(sp, acc)}


def main():
    cat = load_catalogue()
    contrib, matched = {}, {}
    for name, (fn, filt) in MATRICES.items():
        means = load_means(fn, filt)
        hit = [k for k, m in means.items() if k in cat and m > 0]
        contrib[name] = len(hit)
        matched[name] = sum(1 for k in means if k in cat)
        print('%-28s matched=%2d contributing=%2d' % (name, matched[name], contrib[name]))

    rows = list(csv.DictReader(open(TSV, encoding='utf-8'), delimiter='\t'))
    fields = list(rows[0].keys())
    for col in ('contributing_species', 'ratio_nsp'):
        if col not in fields:
            fields.append(col)
    for r in rows:
        name = r['cohort']
        r['contributing_species'] = str(contrib[name])
        r['ratio_nsp'] = r['hdca_matched']       # ratio 就是在这个 n_sp 上算的
    with open(TSV, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print('written:', TSV)
    # 交叉核对：ratio 是否等于 H / H*(ratio_nsp)
    thr = {}
    for r in rows:
        thr.setdefault(r['cohort'], {})[int(r['n_sp'])] = float(r['hstar_median'])
    bad = []
    for r in rows:
        c, ns = r['cohort'], int(r['ratio_nsp'])
        if ns in thr[c]:
            calc = float(r['H']) / thr[c][ns]
            if abs(calc - float(r['ratio'])) > 0.011:
                bad.append((c, r['n_sp'], r['ratio'], round(calc, 3)))
    print('ratio 与 H/H*(ratio_nsp) 不符的格:', bad if bad else '无')


if __name__ == '__main__':
    main()
