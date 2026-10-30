# -*- coding: utf-8 -*-
"""
diag_ratio_spread.py —— 用已有的八种子明细（不重算）核对正文两处极差口径：

  1) 各队列 n_sp=5 处的种子极差（正文写 41–127%）
  2) 5→20 比值自身的种子极差（正文写 43–92% of midpoint）
  3) 全表 40 格极差的上下界与中位数（Table 5 注：6–127%，中位 47%）

数据来源：
  cross_cohort9_seed_t5_detail.txt   九个队列（广州为旧的 122 份口径）
  heq_cd_only.txt                    广州 49 例 CD 口径，替换上一文件的广州行
"""
import os
import re
import statistics

BASE = os.path.dirname(os.path.abspath(__file__))


def parse_detail(path):
    """返回 {cohort: {n_sp: [八个值]}}"""
    out, cur, buf = {}, None, ''
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            if not line.startswith('  '):
                cur = line.split('  ')[0].strip()
                out[cur] = {}
                continue
            m = re.match(r'\s+n_sp=\s*(\d+)\s+(.*?)median=', line)
            if m and cur:
                vals = [float(x) for x in m.group(2).split()]
                out[cur][int(m.group(1))] = vals
    return out


def parse_heq(path):
    d = {}
    with open(path, encoding='utf-8') as f:
        for line in f:
            m = re.match(r'\s+n_sp=\s*(\d+)\s+(.*?)median=', line)
            if m:
                d[int(m.group(1))] = [float(x) for x in m.group(2).split()]
    return d


data = parse_detail(os.path.join(BASE, 'cross_cohort9_seed_t5_detail.txt'))
heq = parse_heq(os.path.join(BASE, 'heq_cd_only.txt'))
for k in list(data):
    if k.startswith('HeQ'):
        del data[k]
data['HeQ2017 (CN, MPA3) 49CD'] = heq

print('=== 1. 各格极差 %d 格 ===' % (len(data) * 4))
cells = []
for c, d in data.items():
    for nsp in (3, 5, 10, 20):
        v = d[nsp]
        rng = (max(v) - min(v)) / statistics.median(v) * 100
        cells.append((c, nsp, rng))
        print('  %-28s n_sp=%2d  %6.1f%%  [%.5f-%.5f]' % (c, nsp, rng, min(v), max(v)))

for nsp in (3, 5, 10, 20):
    sub = [r for _, n, r in cells if n == nsp]
    print('  n_sp=%2d  十队列极差区间 %.0f–%.0f%%' % (nsp, min(sub), max(sub)))
allv = [r for _, _, r in cells]
print('  全 %d 格：%.0f–%.0f%%  中位 %.0f%%'
      % (len(allv), min(allv), max(allv), statistics.median(allv)))

print()
print('=== 2. 5→20 比值自身的种子极差 ===')
spreads = []
for c, d in data.items():
    r = [a / b for a, b in zip(d[20], d[5])]
    mid = (max(r) + min(r)) / 2
    sp = (max(r) - min(r)) / mid * 100
    spreads.append(sp)
    print('  %-28s ratio %.2f–%.2f  mid %.2f  spread %5.1f%%' % (c, min(r), max(r), mid, sp))
print('  十队列：%.0f–%.0f%%（中位 %.0f%%）' % (min(spreads), max(spreads), statistics.median(spreads)))

print()
print('=== 3. 方向核查：H20 > H5 的种子数 ===')
for c, d in data.items():
    k = sum(1 for a, b in zip(d[20], d[5]) if a > b)
    med = statistics.median(d[20]) > statistics.median(d[5])
    print('  %-28s 八个种子中 %d 个成立，中位数成立=%s' % (c, k, med))
