# -*- coding: utf-8 -*-
"""下载与已取得 HUMAnN3 样本配对的 HUMAnN2 包，用于同一样本的 v2/v3 对比。

依据（2026-11-11 实测）：HMP2 同时提供
  2018-05-04/func_profiles/<SID>_humann2.tar.bz2   HUMAnN2（旧）
  func_profile3_WGS/<SID>_humann3.tar.bz2          HUMAnN3（新）
两套样本 ID 的基础名重叠 1,317 个，v3 名字带 _P 后缀、v2 不带，
故配对规则：v3 样本名去掉 _P/_TR 后缀即 v2 的 SID。

坑（2026-11-11 15:00 更正）：v2 路径必须带 `2018-05-04/` 前缀。早先漏了这一段，
全部 404（137 字节错误体），被误判为「Globus 端点故障」——端点从未故障。
"""
import os
import re
import subprocess
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
DEST = os.path.join(BASE, 'hmp2_h3')
IDS = os.path.join(BASE, 'v2_ids.txt')
V2 = ('https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/'
      '2018-05-04/func_profiles/%s_humann2.tar.bz2')
MIN = 300000


def main():
    # 只处理 v2 与 v3 的真实交集（1,317 个）。v3 有 1,638 个样本，其中 321 个
    # 在 v2 目录里没有对应包（实测 CSM5FZ3N/R/T 即属此类）——不是下载失败。
    have2 = set()
    if os.path.exists(IDS):
        have2 = {l.strip() for l in open(IDS) if l.strip()}
    v3 = sorted(f[:-len('_humann3.tar.bz2')] for f in os.listdir(DEST)
                if f.endswith('_humann3.tar.bz2'))
    bases = []
    for s in v3:
        b = re.sub(r'_(P|TR)$', '', s)
        if have2 and b not in have2:
            continue
        if b not in bases:
            bases.append(b)
    print('v3 样本 %d，其中可配对 %d' % (len(v3), len(bases)), flush=True)
    ok = fail = cached = 0
    for i, sid in enumerate(bases, 1):
        path = os.path.join(DEST, '%s_humann2.tar.bz2' % sid)
        if os.path.exists(path) and os.path.getsize(path) > MIN:
            cached += 1
            continue
        good = False
        for k in range(2):
            if os.path.exists(path):
                os.remove(path)
            subprocess.run(['curl', '-sL', '-m', '240', '-o', path, V2 % sid],
                           capture_output=True)
            if os.path.exists(path) and os.path.getsize(path) > MIN:
                t = subprocess.run(['tar', '-tjf', path], capture_output=True)
                if t.returncode == 0:
                    good = True
                    break
                os.remove(path)
        if good:
            ok += 1
            print('[%d/%d] %-14s %.2f MB' % (i, len(bases), sid,
                  os.path.getsize(path) / 1e6), flush=True)
        else:
            fail += 1
            print('[%d/%d] %-14s FAIL' % (i, len(bases), sid), flush=True)
    print('完成 ok=%d cached=%d fail=%d' % (ok, cached, fail))


if __name__ == '__main__':
    main()
