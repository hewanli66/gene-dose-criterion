# -*- coding: utf-8 -*-
"""
cross_cohort8.py —— 每条流程各补一个独立队列后的六队列 H* 重算。

队列（按物种定量流程分组）：
  MetaPhlAn 3 : ELDERMET（爱尔兰，189 例） + CRC-9（1262 例） + HeQ2017（中国广东，122 例）
  MetaPhlAn 4 : SG90（新加坡华人，516 例） + Villages（洪都拉斯，1787 例） + Yachida2019（日本，615 例）
  mOTU        : Wave3（中国长寿带，59 例） + Zhu2020（西安，171 例）

同一流程内两个队列互相验证，跨流程再验证一次；这样"H* 是队列属性而不是流程假象"
这句话才有依据。

输出：cross_cohort8.tsv（表格）+ 屏幕摘要
运行：python cross_cohort8.py
"""
import csv
import os
import random
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import review_final as rf

BASE = os.path.dirname(os.path.abspath(__file__))
NS_LIST = [3, 5, 10, 20]
N_DRAW = 12
SEED = 20261029


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return 0.0


def load_matrix(path, keep_ids=None):
    """样本 × 物种 CSV（首列 SampleID）。"""
    with open(path, encoding='utf-8-sig') as f:
        rd = csv.reader(f)
        hdr = next(rd)
        cols = hdr[1:]
        rows = [r for r in rd if r and r[0].strip()
                and (keep_ids is None or r[0].strip() in keep_ids)]
    n = len(rows)
    info = []
    for j, c in enumerate(cols):
        v = [fnum(r[j + 1]) for r in rows]
        info.append(dict(name=c, mean=sum(v) / n,
                         prev=sum(1 for x in v if x > 0) / n, vals=v))
    return info, n, [r[0].strip() for r in rows]


def subset(info, idx):
    out = []
    n = len(idx)
    for x in info:
        v = [x['vals'][i] for i in idx]
        out.append(dict(name=x['name'], mean=sum(v) / n,
                        prev=sum(1 for y in v if y > 0) / n, vals=v))
    return out


def summarise(tag, info, car_vals, quiet=False):
    n_ind = len(info[0]['vals']) if info else 0
    pool = [x for x in info if x['prev'] >= 0.20 and x['mean'] > 0]
    out = dict(tag=tag, n_ind=n_ind, species=len(info), pool=len(pool), ns={})
    if not quiet:
        print('%-26s n_ind=%4d  species=%4d  pool=%3d'
              % (tag, n_ind, len(info), len(pool)))
    for n_sp in NS_LIST:
        if len(pool) < n_sp:
            continue
        fam = rf.families(car_vals, pool, n_sp, n_ind)
        out['ns'][n_sp] = rf.h_cond(fam)
        if not quiet:
            print('    n_sp=%2d families=%4d H*=%.5f' % (n_sp, len(fam), out['ns'][n_sp]))
    return out


def match_hdca(info, carr):
    by = {x['name'].lower().replace('_', ' '): x for x in info}
    return [(sp, by[sp.lower().replace('_', ' ')], d['frac'])
            for sp, d in carr.items() if sp.lower().replace('_', ' ') in by]


def interp(ns, n_sp):
    import math
    ks = sorted(ns)
    if n_sp in ns:
        return ns[n_sp], True
    if n_sp <= ks[0]:
        return ns[ks[0]], False
    if n_sp >= ks[-1]:
        return ns[ks[-1]], False
    lo = max(k for k in ks if k <= n_sp)
    hi = min(k for k in ks if k >= n_sp)
    t = (math.log(n_sp) - math.log(lo)) / (math.log(hi) - math.log(lo))
    return math.exp(math.log(ns[lo]) + t * (math.log(ns[hi]) - math.log(ns[lo]))), False


def hdca_state(tag, info, carr, res, quiet=False):
    hit = match_hdca(info, carr)
    if not hit:
        return None
    n_ind = len(info[0]['vals'])
    tot = sum(x['mean'] for _, x, _ in hit) or 1.0
    ws = [x['mean'] / tot for _, x, _ in hit]
    cs = [c for _, _, c in hit]
    c_eff = sum(w * c for w, c in zip(ws, cs))
    H = sum(w * (c - c_eff) ** 2 for w, c in zip(ws, cs))
    thr, exact = interp(res['ns'], len(hit))
    # 实测保真度：全体 / 仅携带者
    p = [sum(x['vals'][i] for _, x, _ in hit) for i in range(n_ind)]
    c = [sum(x['vals'][i] * f for _, x, f in hit) for i in range(n_ind)]
    rho_all = rf.spearman(p, c)
    pos = [i for i, v in enumerate(p) if v > 0]
    rho_pos = rf.spearman([p[i] for i in pos], [c[i] for i in pos]) if len(pos) > 5 else float('nan')
    st = dict(tag=tag, matched=len(hit), c_eff=c_eff, H=H, thr=thr, exact=exact,
              ratio=H / thr, rho_all=rho_all, rho_pos=rho_pos,
              carrier_frac=len(pos) / n_ind,
              top_w=max(ws))
    if not quiet:
        print('%-26s hdcA 命中=%2d (携带者 %3d/%3d=%.0f%%)  c_eff=%.3f  H=%.4f  H*(%d)=%.5f  H/H*=%.2f  rho=%.2f/%.2f'
              % (tag, len(hit), len(pos), n_ind, 100 * st['carrier_frac'],
                 c_eff, H, len(hit), thr, st['ratio'], rho_all, rho_pos))
    return st


def main():
    carr, car_vals, _, _, _, _ = rf.load()
    rnd = random.Random(SEED)

    spec = []
    p = lambda *a: os.path.join(BASE, *a)
    spec.append(('ELDERMET (IE, MPA3)', p('ELDERMET_species_189x836.csv'), None))
    spec.append(('CRC9 pool (MPA3)', p('CRC9_samples_1262x917.csv'), None))
    spec.append(('SG90 (SG, MPA4)', p('SG90_species_abundances.csv'), None))
    spec.append(('Villages (HN, MPA4)', p('Village_samples_1787x1910.csv'), None))
    spec.append(('Wave3 65+ (CN, mOTU)', p('Wave3_species_59x.csv'), None))
    spec.append(('Zhu2020 (CN, mOTU)', p('Zhu2020_samples_171x359.csv'), None))
    spec.append(('HeQ2017 (CN, MPA3)', p('HeQ2017_samples_122x.csv'), None))
    spec.append(('Yachida2019 (JP, MPA4)', p('Yachida2019_samples_615x1835.csv'), None))

    loaded = {}
    print('=== 各队列 H* ===')
    for tag, path, keep in spec:
        if not os.path.exists(path):
            print('%-26s [缺文件] %s' % (tag, os.path.basename(path)))
            continue
        info, n, ids = load_matrix(path, keep)
        loaded[tag] = (info, n, ids)
        summarise(tag, info, car_vals)

    print('\n=== hdcA 家族状态 ===')
    states = []
    res = {}
    for tag, (info, n, ids) in loaded.items():
        res[tag] = summarise(tag, info, car_vals, quiet=True)
        st = hdca_state(tag, info, carr, res[tag])
        if st:
            states.append(st)

    # --- 样本量匹配 ---
    if loaded:
        nm = min(v[1] for v in loaded.values())
        print('\n=== 样本量匹配：全部抽到 n=%d（%d 次）===' % (nm, N_DRAW))
        for tag, (info, n, ids) in loaded.items():
            vals = []
            for _ in range(N_DRAW):
                idx = rnd.sample(range(n), nm)
                r = summarise(tag, subset(info, idx), car_vals, quiet=True)
                vals.append(r['ns'])
            line = '%-26s' % tag
            for n_sp in NS_LIST:
                xs = [v[n_sp] for v in vals if n_sp in v]
                if xs:
                    line += '  H*(%d)=%.5f [%.5f-%.5f]' % (n_sp, statistics.median(xs), min(xs), max(xs))
            print(line)

    # --- 输出 ---
    with open(os.path.join(BASE, 'cross_cohort8.tsv'), 'w', encoding='utf-8') as f:
        f.write('cohort\tn_ind\tspecies\tpool\t' +
                '\t'.join('H%d' % k for k in NS_LIST) +
                '\tmatched\tceff\tH\tthr\tratio\trho_all\trho_carriers\tcarrier_frac\ttop_w\n')
        for tag, (info, n, ids) in loaded.items():
            r = res[tag]
            st = next((s for s in states if s['tag'] == tag), None)
            row = [tag, str(n), str(len(info)), str(r['pool'])]
            row += ['%.5f' % r['ns'][k] if k in r['ns'] else '' for k in NS_LIST]
            if st:
                row += [str(st['matched']), '%.3f' % st['c_eff'], '%.4f' % st['H'],
                        '%.5f' % st['thr'], '%.2f' % st['ratio'],
                        '%.3f' % st['rho_all'], '%.3f' % st['rho_pos'],
                        '%.3f' % st['carrier_frac'], '%.3f' % st['top_w']]
            f.write('\t'.join(row) + '\n')
    print('\n写出 cross_cohort8.tsv')


if __name__ == '__main__':
    main()
