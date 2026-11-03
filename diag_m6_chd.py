# -*- coding: utf-8 -*-
"""N-F 诊断：36 个"c_eff 放行但 ρ<0.90"的子族为何全部来自冠状动脉队列。

对每个队列、k=3/4，枚举真实 hdcA 携带物种的子族，统计：
  - c_eff > 0.84 的子族数（这是 c_eff 判据可能放行的前提）
  - 其中 ρ < 0.90 的个数（真正的不一致案例）
  - 该队列 c_eff 的最大值，以及"最不齐"的亚族里 0<c<1 物种的权重占比

只读 数据/ 与 投稿包/，不写任何文件。
"""
import itertools
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m6_discriminant as m6
import cross_cohort as cc

N_SUB = 250
KS = (3, 4, 5, 6, 8, 10)


def main():
    carr, _ = m6.load_catalogue()
    hstar = m6.load_hstar()
    rnd = random.Random(20261103)

    print('%-28s %-2s %6s %6s %6s %7s %6s' %
          ('cohort', 'k', 'n_sub', 'ceff>84', 'rho<90', 'max_ce', 'n_hit'))
    print('-' * 78)
    for name, fn in m6.COHORTS:
        coh = m6.load_cohort(fn)
        if coh is None:
            print('%-28s  (matrix missing)' % name)
            continue
        hit = cc.match_hdca(coh['info'], carr)
        name2j = {x: j for j, x in enumerate(coh['names_all'])}
        usable = [(x['name'], frac) for _, x, frac in hit
                  if x['name'] in name2j and coh['mean_all'][name2j[x['name']]] > 0]
        if len(usable) < 3:
            print('%-28s  (only %d usable carriers)' % (name, len(usable)))
            continue
        cols = [name2j[nm] for nm, _ in usable]
        fr = [f for _, f in usable]
        sub_all = coh['Aall'][:, cols]
        mn = coh['mean_all'][cols]
        m = len(usable)
        for k in KS:
            if k > m:
                continue
            combos = list(itertools.combinations(range(m), k))
            if len(combos) > N_SUB:
                combos = rnd.sample(combos, N_SUB)
            n_hi = n_bad = 0
            max_ce = 0.0
            worst = None
            for cb in combos:
                cb = list(cb)
                ws_raw = mn[cb]
                if ws_raw.sum() <= 0:
                    continue
                ws = ws_raw / ws_raw.sum()
                cs = np.array([fr[i] for i in cb])
                ce, H = m6.moments(ws, cs)
                if ce > 0.84:
                    n_hi += 1
                    max_ce = max(max_ce, ce)
                    sub = sub_all[:, cb]
                    rho = m6.spear(sub.sum(1), sub @ cs)
                    if rho < 0.90:
                        n_bad += 1
                        if worst is None or ce > worst[0]:
                            worst = (ce, [(usable[i][0], fr[i], w)
                                          for i, w in zip(cb, ws)])
            print('%-28s %-2d %6d %6d %6d %7.3f %6d' %
                  (name, k, len(combos), n_hi, n_bad, max_ce, m))
            if worst and worst[0] > 0.84:
                detail = '; '.join('%s c=%.2f w=%.2f' % (s, c, w)
                                   for s, c, w in worst[1])
                print('      worst c_eff=%.3f : %s' % (worst[0], detail))


if __name__ == '__main__':
    main()
