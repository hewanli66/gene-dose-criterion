# -*- coding: utf-8 -*-
"""批量下载 HMP2 的 HUMAnN3 样本包（并发 + 短超时 + 完整性校验）。

URL（2026-11-11 实测可达）：
  https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/
      func_profile3_WGS/<SID>_humann3.tar.bz2        约 1.3 MB/样本
  .../HMP2/MGX/2018-05-04/taxonomic_profiles_3.tsv.gz  MetaPhlAn 3 物种表

坑：Globus 对个别样本会挂住不返回，串行下载会整体停摆（第一版 -m 300 × 5 次
重试，4 分钟零进展）。故改并发 4 + 单次 90 秒上限 + 失败即换下一个。
"""
import concurrent.futures as cf
import csv
import gzip
import os
import random
import re
import subprocess
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
HTML = os.path.join(BASE, 'mgx_products.html')
SPTAB = os.path.join(BASE, 'hmp2_h3', 'taxonomic_profiles_3.tsv.gz')
DEST = os.path.join(BASE, 'hmp2_h3')
URL = ('https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/'
       'func_profile3_WGS/%s_humann3.tar.bz2')
MIN_BYTES = 200000


def species_samples():
    with gzip.open(SPTAB, 'rt', encoding='utf-8', errors='replace') as f:
        hdr = next(csv.reader(f, delimiter='\t'))
    return {h.replace('_profile', '') for h in hdr[1:] if h.strip()}


def one(sid):
    path = os.path.join(DEST, '%s_humann3.tar.bz2' % sid)
    if os.path.exists(path) and os.path.getsize(path) > MIN_BYTES:
        return sid, True, 'cached'
    for k in range(3):
        if os.path.exists(path):
            os.remove(path)
        subprocess.run(['curl', '-sL', '-m', '90', '-o', path, URL % sid],
                       capture_output=True)
        if os.path.exists(path) and os.path.getsize(path) > MIN_BYTES:
            t = subprocess.run(['tar', '-tjf', path], capture_output=True)
            if t.returncode == 0:
                return sid, True, 'ok'
            # 半截文件：tar 报损坏
            if os.path.exists(path):
                os.remove(path)
    return sid, False, 'fail'


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261111
    html = open(HTML, encoding='utf-8', errors='replace').read()
    allsids = sorted(set(re.findall(r'func_profile3_WGS/([A-Z0-9_]+)_humann3\.tar\.bz2', html)))
    have = species_samples()
    cand = [s for s in allsids if s in have]
    rnd = random.Random(seed)
    rnd.shuffle(cand)
    picked = cand[:n]
    print('候选 %d，本次取 %d' % (len(cand), len(picked)), flush=True)
    os.makedirs(DEST, exist_ok=True)
    ok = fail = cached = 0
    with cf.ThreadPoolExecutor(6) as ex:
        for i, (sid, good, msg) in enumerate(ex.map(one, picked), 1):
            if good:
                ok += 1
                if msg == 'cached':
                    cached += 1
            else:
                fail += 1
            print('[%d/%d] %-14s %s' % (i, len(picked), sid, msg), flush=True)
    print('完成 ok=%d（其中缓存 %d）fail=%d' % (ok, cached, fail), flush=True)


if __name__ == '__main__':
    main()
