# -*- coding: utf-8 -*-
"""把 v2 -> v3 的判决变化拆成「家族规模 n_sp 变了」与「族内强度分布 H 变了」两半。

判据是 H / H*(n_sp) > 1 则 fail。v3 相对于 v2 同时改了两样东西：
  (a) n_sp  —— 哪些物种被判为该家族的携带者
  (b) H     —— 在同一批物种内部，相对强度分布齐不齐
只有把 H*(·) 固定在同一个 n_sp 上，才能分清翻转来自哪一半。

三个量：
  ratio2   = H2 / H*(n2)           v2 的实际判决量
  ratio3   = H3 / H*(n3)           v3 的实际判决量
  ratio3n2 = H3 / H*(n2)           反事实：H 用 v3 的，阈值仍按 v2 的家族规模取

于是
  翻转(ratio3, ratio2)    总翻转
  翻转(ratio3n2, ratio2)  仅由 H 变化引起的翻转
  翻转(ratio3, ratio3n2)  仅由 n_sp 变化引起的翻转
两者不是简单相加（同一个家族可能两半都翻，方向相反会抵消），所以都报出来，并报
「只有 n_sp 翻」「只有 H 翻」「两者同向翻」「两者反向相互抵消」四类的计数。

用法：python decompose_v2v3.py [--rel 0.001] [--min-p 0.1] [--min-sp 3]
"""
import argparse
import glob
import math
import os
import re
import statistics
from collections import Counter

from pair_v2v3 import hstar_at, load_species, read_fam, stats, spearman

BASE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(BASE, 'work')
OUT = os.path.join(BASE, 'out')


def base_name(s):
    return re.sub(r'_(P|TR)$', '', s)


def collect(a):
    sp = load_species()
    v3files = {os.path.basename(p).replace('_genefamilies_cpm.tsv', ''): p
               for p in glob.glob(os.path.join(WORK, '*_genefamilies_cpm.tsv'))}
    v2files = {os.path.basename(p).replace('_genefamilies.tsv', ''): p
               for p in glob.glob(os.path.join(WORK, '*_genefamilies.tsv'))
               if not p.endswith('_cpm.tsv')}
    v3by = {base_name(k): k for k in v3files}
    pairs = [(v3by[b], b) for b in sorted(set(v3by) & set(v2files))]

    recs = []
    for sid3, sid2 in pairs:
        ptab = sp.get(sid3) or sp.get(sid2)
        if not ptab:
            continue
        f3, f2 = read_fam(v3files[sid3]), read_fam(v2files[sid2])
        for fam in set(f3) & set(f2):
            s3 = stats(f3[fam], ptab, a.rel, a.min_p, a.min_sp)
            s2 = stats(f2[fam], ptab, a.rel, a.min_p, a.min_sp)
            if not s3 or not s2:
                continue
            recs.append(dict(sample=sid3, fam=fam,
                             n2=s2[0], n3=s3[0], H2=s2[2], H3=s3[2],
                             ratio2=s2[2] / hstar_at(s2[0]),
                             ratio3=s3[2] / hstar_at(s3[0]),
                             ratio3n2=s3[2] / hstar_at(s2[0])))
    return recs


def report(recs, tag=''):
    n = len(recs)
    print('=' * 72)
    print('参数 %s  家族-样本对 n=%d' % (tag, n))
    same = [r for r in recs if r['n2'] == r['n3']]
    print('n_sp 在 v2/v3 间相同：%d / %d = %.1f%%' % (len(same), n, 100 * len(same) / n))
    dn = Counter(r['n3'] - r['n2'] for r in recs)
    print('Δn_sp 分布：' + '  '.join('%+d:%d(%.0f%%)' % (k, v, 100 * v / n)
                                    for k, v in sorted(dn.items()) if abs(k) <= 4))
    print('n_sp 中位：v2 = %.0f，v3 = %.0f' %
          (statistics.median(r['n2'] for r in recs),
           statistics.median(r['n3'] for r in recs)))

    # 判决
    f2 = [r['ratio2'] > 1 for r in recs]
    f3 = [r['ratio3'] > 1 for r in recs]
    f3n2 = [r['ratio3n2'] > 1 for r in recs]
    tot = sum(1 for x, y in zip(f2, f3) if x != y)
    by_h = sum(1 for x, y in zip(f2, f3n2) if x != y)
    by_n = sum(1 for x, y in zip(f3n2, f3) if x != y)
    print('总翻转 %d (%.1f%%)；仅 H 变引起的 %d (%.1f%%)；仅 n_sp 变引起的 %d (%.1f%%)'
          % (tot, 100 * tot / n, by_h, 100 * by_h / n, by_n, 100 * by_n / n))
    cls = Counter()
    for x, y, z in zip(f2, f3n2, f3):
        if x == y == z:
            cls['未翻'] += 1
        elif x != y and y == z:
            cls['仅H'] += 1
        elif x == y and y != z:
            cls['仅n_sp'] += 1
        elif x != y and y != z:
            cls['两段各翻一次'] += 1
    print('判决轨迹分类：' + '  '.join('%s:%d(%.1f%%)' % (k, v, 100 * v / n)
                                       for k, v in cls.most_common()))

    # 翻转是否集中在边界
    print('\n按 v2 的判决量分箱（|log ratio2|，越接近 1 越在边界）：')
    print('  %-14s %6s %8s %8s' % ('ratio2 区间', 'n', '总翻转%', '仅n_sp%'))
    bins = [(0, 0.5), (0.5, 0.8), (0.8, 1.0), (1.0, 1.25), (1.25, 2.0), (2.0, 1e9)]
    for lo, hi in bins:
        sub = [i for i, r in enumerate(recs) if lo <= r['ratio2'] < hi]
        if not sub:
            continue
        t = sum(1 for i in sub if f2[i] != f3[i])
        nn = sum(1 for i in sub if f3n2[i] != f3[i])
        print('  %-14s %6d %7.1f%% %7.1f%%'
              % ('%g-%g' % (lo, hi) if hi < 1e8 else '>2', len(sub),
                 100 * t / len(sub), 100 * nn / len(sub)))

    # Δn_sp 分组下的 H 变化
    print('\n按 Δn_sp 分组看 H 与 H 的相对变化：')
    print('  %-8s %6s %10s %12s %10s' % ('Δn_sp', 'n', 'ΔH中位', 'ΔH/H2中位%', 'ratio中位变化'))
    groups = [('<=-1', -99, -1), ('0', 0, 0), ('+1', 1, 1), ('+2', 2, 2), ('>=+3', 3, 99)]
    for name, lo, hi in groups:
        sub = [r for r in recs if lo <= r['n3'] - r['n2'] <= hi]
        if len(sub) < 5:
            continue
        dh = [r['H3'] - r['H2'] for r in sub]
        rel = [100 * (r['H3'] - r['H2']) / r['H2'] for r in sub if r['H2'] > 0]
        print('  %-8s %6d %10.4f %11.1f%% %10.3f'
              % (name, len(sub), statistics.median(dh),
                 statistics.median(rel),
                 statistics.median([r['ratio3'] - r['ratio2'] for r in sub])))

    # H* 的队列间不确定性：把阈值整体缩放后翻转率怎么变
    print('\n阈值整体缩放（模拟换队列用 H*）后的总翻转率：')
    for k in (0.79, 1.0, 1.26):
        ff = 0
        for r in recs:
            a2 = r['H2'] / (hstar_at(r['n2']) * k) > 1
            a3 = r['H3'] / (hstar_at(r['n3']) * k) > 1
            ff += (a2 != a3)
        print('  H* x %.2f -> %.1f%%' % (k, 100 * ff / n))
    return recs


def extra(recs):
    """n_sp 是不是通过抬高 H 才起作用的？按 Δn_sp 分层看翻转。"""
    n = len(recs)
    print('\n' + '=' * 72)
    print('按 Δn_sp 分层的翻转（若 n_sp 是真驱动，这一层应显著分化）')
    print('  %-8s %6s %8s %10s %10s %10s'
          % ('Δn_sp', 'n', '翻转%', 'fail->pass', 'pass->fail', 'H2中位'))
    for name, lo, hi in [('<=-1', -99, -1), ('0', 0, 0), ('+1', 1, 1),
                         ('+2', 2, 2), ('>=+3', 3, 99)]:
        sub = [r for r in recs if lo <= r['n3'] - r['n2'] <= hi]
        if not sub:
            continue
        fp = sum(1 for r in sub if r['ratio2'] > 1 and r['ratio3'] <= 1)
        pf = sum(1 for r in sub if r['ratio2'] <= 1 and r['ratio3'] > 1)
        print('  %-8s %6d %7.1f%% %10d %10d %10.4f'
              % (name, len(sub), 100 * (fp + pf) / len(sub), fp, pf,
                 statistics.median(r['H2'] for r in sub)))

    print('\n整体翻转方向：fail->pass %d，pass->fail %d'
          % (sum(1 for r in recs if r['ratio2'] > 1 and r['ratio3'] <= 1),
             sum(1 for r in recs if r['ratio2'] <= 1 and r['ratio3'] > 1)))

    print('\n限制在 H2 不太小的家族上（H2 很小时 ratio 数值不稳，翻转无意义）')
    print('  %-14s %6s %8s %10s' % ('条件', 'n', '翻转%', '仅n_sp%'))
    for floor in (0.0, 0.002, 0.005, 0.01, 0.02):
        sub = [r for r in recs if r['H2'] >= floor]
        if not sub:
            continue
        t = sum(1 for r in sub if (r['ratio2'] > 1) != (r['ratio3'] > 1))
        nn = sum(1 for r in sub if (r['ratio3n2'] > 1) != (r['ratio3'] > 1))
        print('  H2 >= %-8g %6d %7.1f%% %9.1f%%'
              % (floor, len(sub), 100 * t / len(sub), 100 * nn / len(sub)))

    print('\nratio2 分箱内 H2 的量级（说明高 ratio 箱是否只是 H 极小的噪声）')
    print('  %-12s %6s %10s %10s %8s' % ('ratio2', 'n', 'H2中位', 'ratio3中位', '翻转%'))
    for lo, hi in [(0, 0.5), (0.5, 0.8), (0.8, 1.0), (1.0, 1.25),
                   (1.25, 2.0), (2.0, 5.0), (5.0, 1e9)]:
        sub = [r for r in recs if lo <= r['ratio2'] < hi]
        if not sub:
            continue
        t = sum(1 for r in sub if (r['ratio2'] > 1) != (r['ratio3'] > 1))
        print('  %-12s %6d %10.5f %10.3f %7.1f%%'
              % ('%g-%g' % (lo, hi) if hi < 1e8 else '>5', len(sub),
                 statistics.median(r['H2'] for r in sub),
                 statistics.median(r['ratio3'] for r in sub),
                 100 * t / len(sub)))


def crosstab(recs):
    """控制「离判决边界多远」之后，Δn_sp 是否还影响翻转方向。

    翻转概率强烈依赖 ratio2 距 1 的距离（见上）。若 Δn_sp≥+1 的家族恰好都挤在
    边界附近，分层翻转率差异就只是构成效应。这里在每个 ratio2 箱内比较。
    """
    print('\n' + '=' * 72)
    print('控制边界距离：同一 ratio2 箱内，Δn_sp 组之间的翻转率与方向')
    print('  %-12s %-8s %6s %7s %8s %8s'
          % ('ratio2', 'Δn_sp', 'n', '翻转%', 'f->p', 'p->f'))
    for lo, hi in [(0, 0.8), (0.8, 1.0), (1.0, 1.25), (1.25, 2.0), (2.0, 1e9)]:
        for name, a, b in [('0', 0, 0), ('>=+1', 1, 99), ('<=-1', -99, -1)]:
            sub = [r for r in recs
                   if lo <= r['ratio2'] < hi and a <= r['n3'] - r['n2'] <= b]
            if len(sub) < 10:
                continue
            fp = sum(1 for r in sub if r['ratio2'] > 1 and r['ratio3'] <= 1)
            pf = sum(1 for r in sub if r['ratio2'] <= 1 and r['ratio3'] > 1)
            print('  %-12s %-8s %6d %6.1f%% %8d %8d'
                  % ('%g-%g' % (lo, hi) if hi < 1e8 else '>2', name, len(sub),
                     100 * (fp + pf) / len(sub), fp, pf))
        print('')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rel', type=float, default=0.001)
    ap.add_argument('--min-p', type=float, default=0.1)
    ap.add_argument('--min-sp', type=int, default=3)
    a = ap.parse_args()
    recs = report(collect(a), 'rel=%g min_p=%g' % (a.rel, a.min_p))
    extra(recs)
    crosstab(recs)
    print('\n分样本：')
    print('  %-14s %6s %8s %10s %10s %8s %8s'
          % ('sample', 'n', 'nsp同%', 'H2中位', 'H3中位', '翻转%', '仅nsp%'))
    for s in sorted({r['sample'] for r in recs}):
        sub = [r for r in recs if r['sample'] == s]
        f2 = [r['ratio2'] > 1 for r in sub]
        f3 = [r['ratio3'] > 1 for r in sub]
        f3n2 = [r['ratio3n2'] > 1 for r in sub]
        print('  %-14s %6d %7.0f%% %10.4f %10.4f %7.1f%% %7.1f%%'
              % (s, len(sub),
                 100 * sum(1 for r in sub if r['n2'] == r['n3']) / len(sub),
                 statistics.median(r['H2'] for r in sub),
                 statistics.median(r['H3'] for r in sub),
                 100 * sum(1 for x, y in zip(f2, f3) if x != y) / len(sub),
                 100 * sum(1 for x, y in zip(f3n2, f3) if x != y) / len(sub)))
    print('\nρ(H2,H3) 全体 = %.3f' % spearman([r['H2'] for r in recs],
                                              [r['H3'] for r in recs]))

    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, 'decompose_v2v3.tsv')
    cols = ['sample', 'fam', 'n2', 'n3', 'H2', 'H3', 'ratio2', 'ratio3', 'ratio3n2']
    with open(p, 'w', encoding='utf-8', newline='') as f:
        f.write('\t'.join(cols) + '\n')
        for r in recs:
            f.write('\t'.join(('%.6g' % r[c]) if isinstance(r[c], float) else str(r[c])
                              for c in cols) + '\n')
    print('写入 %s' % p)


if __name__ == '__main__':
    main()
