# -*- coding: utf-8 -*-
"""
cross_cohort.py —— 在第二个独立队列（SG90，新加坡华人）上重算 H* 阈值。

目的：回应审稿意见 M8——“三个子集同属 ELDERMET，证明的是队列内稳健，
不是跨队列可迁移”。

队列：
  A. ELDERMET（爱尔兰，189 例，836 物种谱）—— 原队列，用于复现 0.0343/0.0426/0.0769/0.0923
  B. SG90 全部（新加坡，516 例，382 物种谱）
  C. SG90 八旬老人子集（年龄 ≥77，221 例）—— 与 ELDERMET 年龄段最可比

方法：与 review_final.py 完全相同的模拟（同 seed、同 1,200 族/规模、同条件阈值定义），
只是把丰度矩阵换成另一队列，并把物种抽样池换成该队列 prevalence ≥20% 的物种
（ELDERMET 149 个；SG90 全部 69 个；SG90 老年 63 个）。

注意：carriage 率 c_s 仍从 Mou 2021 的 96 个观测值里抽——携带率属于基因目录，
不随队列变化；队列只改变丰度结构与个体组成。

运行环境：Python 3.8+，仅标准库。约 3–5 分钟。
输入：数据/ELDERMET_species_189x836.csv、数据/SG90_species_abundances.csv、
      数据/SG90_metadata.csv、数据/TableS1_species_fixed.csv
输出：数据/cross_cohort.txt、数据/cross_cohort.tsv
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import review_final as rf

BASE = os.path.dirname(os.path.abspath(__file__))
NS_LIST = [3, 5, 10, 20]
AGE_CUT = 77


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return 0.0


def load_matrix(path, keep_ids=None):
    """读 species × sample 矩阵；keep_ids 为 None 时全取。返回 info, n_ind, cols。"""
    with open(path, encoding='utf-8-sig') as f:
        rd = csv.reader(f)
        hdr = next(rd)
        cols = hdr[1:]
        rows = [r for r in rd if keep_ids is None or r[0] in keep_ids]
    n = len(rows)
    info = []
    for j, c in enumerate(cols):
        v = [fnum(r[j + 1]) for r in rows]
        info.append(dict(name=c, mean=sum(v) / n,
                         prev=sum(1 for x in v if x > 0) / n, vals=v))
    return info, n, cols


def pool_of(info, thr=0.20):
    return [x for x in info if x['prev'] >= thr and x['mean'] > 0]


def summarise(tag, info, n_ind, car_vals):
    pool = pool_of(info)
    out = dict(tag=tag, n_ind=n_ind, pool=len(pool), ns={}, cum={})
    line = '%-26s n_ind=%4d  species=%3d  pool(prev>=20%%)=%3d' % (
        tag, n_ind, len(info), len(pool))
    print(line)
    if len(pool) < max(NS_LIST):
        print('  !! 池小于最大 n_sp=%d，仅算可行规模' % max(NS_LIST))
    for n_sp in NS_LIST:
        if len(pool) < n_sp:
            continue
        fam = rf.families(car_vals, pool, n_sp, n_ind)
        hc = rf.h_cond(fam)
        hu = rf.h_cum(fam)
        # 该规模下保真率的中位数（描述这个队列的 ρ–H 关系形状）
        mid = fam[len(fam) // 2][1]
        out['ns'][n_sp] = hc
        out['cum'][n_sp] = hu
        print('   n_sp=%2d  families=%4d  H*_cond=%.5f  H*_cum=%.5f  median_rho=%.3f'
              % (n_sp, len(fam), hc, hu, mid))
    return out


def match_hdca(info, carr):
    """hdcA 目录物种在本队列中能匹配到的那些（名称下划线/空格归一）。"""
    by_key = {}
    for x in info:
        by_key[x['name'].lower().replace('_', ' ')] = x
    hit = []
    for sp in carr:
        k = sp.lower().replace('_', ' ')
        if k in by_key:
            hit.append((sp, by_key[k], carr[sp]['frac']))
    return hit


def interp_hstar(ns, n_sp):
    """按正文口径在 log H* vs log n_sp 上线性插值；超出表范围则取端点。"""
    import math
    ks = sorted(ns)
    if n_sp <= ks[0]:
        return ns[ks[0]], ks[0]
    if n_sp >= ks[-1]:
        return ns[ks[-1]], ks[-1]
    lo = max(k for k in ks if k <= n_sp)
    hi = min(k for k in ks if k >= n_sp)
    if lo == hi:
        return ns[lo], lo
    t = (math.log(n_sp) - math.log(lo)) / (math.log(hi) - math.log(lo))
    v = math.exp(math.log(ns[lo]) + t * (math.log(ns[hi]) - math.log(ns[lo])))
    return v, None   # None 表示插值而来


def hdca_state(tag, info, n_ind, carr, hstar):
    hit = match_hdca(info, carr)
    if not hit:
        print('%-26s hdcA: 无匹配物种' % tag)
        return None
    tot = sum(x['mean'] for _, x, _ in hit) or 1.0
    ws = [x['mean'] / tot for _, x, _ in hit]
    cs = [c for _, _, c in hit]
    c_eff = sum(w * c for w, c in zip(ws, cs))
    H = sum(w * (c - c_eff) ** 2 for w, c in zip(ws, cs))
    n_sp = len(hit)
    thr, exact = interp_hstar(hstar['ns'], n_sp)
    ratio = H / thr if thr else float('nan')
    nz = sum(1 for _, x, _ in hit if x['mean'] > 0)
    tag_thr = ('H*(%d)=%.5f' % (n_sp, thr)) if exact else ('H*(%d, interp)=%.5f' % (n_sp, thr))
    print('%-26s hdcA: matched=%2d  non-zero=%-2d  c_eff=%.3f  H=%.4f  %s  H/H*=%.2f'
          % (tag, n_sp, nz, c_eff, H, tag_thr, ratio))
    return dict(tag=tag, n_sp=n_sp, nz=nz, c_eff=c_eff, H=H, thr=thr, ratio=ratio)


def main():
    carr, car_vals, _, _, _, _ = rf.load()
    print('carriage 观测值 n=%d\n' % len(car_vals))

    met = {}
    with open(os.path.join(BASE, 'SG90_metadata.csv'), encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            met[r['Subject ID']] = fnum(r['Age (in years)'])
    old_ids = {k for k, a in met.items() if a >= AGE_CUT}

    e_info, e_n, _ = load_matrix(os.path.join(BASE, 'ELDERMET_species_189x836.csv'))
    s_info, s_n, _ = load_matrix(os.path.join(BASE, 'SG90_species_abundances.csv'))
    o_info, o_n, _ = load_matrix(os.path.join(BASE, 'SG90_species_abundances.csv'), old_ids)

    res = {}
    for tag, info, n in (('ELDERMET (189)', e_info, e_n),
                         ('SG90 all (%d)' % s_n, s_info, s_n),
                         ('SG90 age>=%d (%d)' % (AGE_CUT, o_n), o_info, o_n)):
        res[tag] = summarise(tag, info, n, car_vals)

    print('\n--- hdcA 家族在各队列的状态 ---')
    states = []
    for tag, info, n in (('ELDERMET (189)', e_info, e_n),
                         ('SG90 all (%d)' % s_n, s_info, s_n),
                         ('SG90 age>=%d (%d)' % (AGE_CUT, o_n), o_info, o_n)):
        st = hdca_state(tag, info, n, carr, res[tag])
        if st:
            states.append(st)

    # 跨队列 H* 比值
    print('\n--- 跨队列比较（同 n_sp） ---')
    e = res['ELDERMET (189)']
    for tag in list(res)[1:]:
        r = res[tag]
        for n_sp in NS_LIST:
            if n_sp in r['ns'] and n_sp in e['ns']:
                print('  n_sp=%2d  ELDERMET=%.5f  %-18s=%.5f  ratio=%.2f'
                      % (n_sp, e['ns'][n_sp], tag, r['ns'][n_sp],
                         r['ns'][n_sp] / e['ns'][n_sp]))

    span_e = max(e['ns'].values()) / min(e['ns'].values())
    print('\n  ELDERMET 内部跨度 = %.2f 倍' % span_e)
    for tag in list(res)[1:]:
        r = res[tag]
        print('  %-18s 内部跨度 = %.2f 倍' % (tag, max(r['ns'].values()) / min(r['ns'].values())))

    with open(os.path.join(BASE, 'cross_cohort.tsv'), 'w', encoding='utf-8') as f:
        f.write('cohort\tn_ind\tpool\tn_sp\tH_star_cond\tH_star_cum\n')
        for tag, r in res.items():
            for n_sp in NS_LIST:
                if n_sp in r['ns']:
                    f.write('%s\t%d\t%d\t%d\t%.5f\t%.5f\n'
                            % (tag, r['n_ind'], r['pool'], n_sp,
                               r['ns'][n_sp], r['cum'][n_sp]))
    print('\nwritten cross_cohort.tsv')


if __name__ == '__main__':
    main()
