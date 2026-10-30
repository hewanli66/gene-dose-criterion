# -*- coding: utf-8 -*-
"""
diag_carriage_frac.py —— 复核正文里"carriers present in X% of individuals"这类句子。

c_eff 是"按丰度加权的平均携带率"，与"个体中有多少比例至少有任一携带物种"
（下面记 frac_pos）不是同一个量。正文混用两处措辞，故把两个都算出来对账。

只算三个队列：广州（49 例 CD）、广州（122 份，旧口径）、冠心病（102 例）。
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import review_final as rf
from cross_cohort9 import load_matrix
from cross_cohort9_seed import hdca_state, pool_of
from cross_cohort9 import match_hdca

BASE = os.path.dirname(os.path.abspath(__file__))


def frac_pos(mat, ids=None):
    carr, car_vals, _, _, _, _ = rf.load()
    info, n, _ = load_matrix(mat, ids) if ids else load_matrix(mat)
    st = hdca_state(info, carr, {20: 0.09}, {20: 0.09})
    # 重新算 p（个体各携带物种丰度之和）
    h = match_hdca(info, carr)
    vals = [sum(x['vals'][i] for _, x, _ in h) for i in range(n)]
    pos = [v for v in vals if v > 0]
    return len(pos), n, len(pos) / n, st['c_eff']


def cd_ids():
    with open(os.path.join(BASE, 'HeQ_2017_sample_metadata.tsv'), encoding='utf-8') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    return {r['sample_id'] for r in rows if r['control'] == 'Case'}


for tag, mat, ids in [
    ('Guangzhou 49 CD', 'HeQ2017_samples_122x.csv', cd_ids()),
    ('Guangzhou 122 all', 'HeQ2017_samples_122x.csv', None),
    ('CHD Feng2016 102', 'CHD_species_102x.csv', None),
]:
    k, n, fr, ce = frac_pos(os.path.join(BASE, mat), ids)
    print('%-18s  n=%4d  pos=%4d  frac_pos=%.3f  c_eff=%.3f' % (tag, n, k, fr, ce))
