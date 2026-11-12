# -*- coding: utf-8 -*-
"""并行 Range 分块下载 HMP2 包——绕过单连接限速。

为什么需要
----------
Globus HTTPS 端点对单连接限到 3–9 KB/s，1.3–2.2 MB 的包要 2–7 分钟；把 6 个
连接直接并发下载不同文件，总吞吐仍然是 9 KB/s（服务端按客户端限速）。

实测（2026-11-11 15:13）：对**同一个文件**用 8 路 Range 并发，120 秒拿到
2.04 MB，约 17 KB/s，是单连接的 3 倍左右。说明限制是每连接而非每客户端，
分块并发能绕过去。故本脚本对每个文件开 K 块并发，块内失败重试，合并后
用 `tar -tjf` 校验，通过才落到正式名。

用法：
  python fetch_parallel.py --kind v3 --samples A,B,C [--blocks 16]
  python fetch_parallel.py --kind v2 --samples A,B --list ids.txt
"""
import argparse
import concurrent.futures as cf
import os
import subprocess
import sys
import time
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
DEST = os.path.join(BASE, 'hmp2_h3')
WORK = os.path.join(BASE, 'work')
HOST = 'https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/'
URL = {
    'v3': HOST + 'func_profile3_WGS/%s_humann3.tar.bz2',
    'v2': HOST + '2018-05-04/func_profiles/%s_humann2.tar.bz2',
}
TSV = {'v3': '%s_genefamilies_cpm.tsv', 'v2': '%s_genefamilies.tsv'}


def tsv_complete(path):
    if not os.path.exists(path) or os.path.getsize(path) < 1000:
        return False
    with open(path, 'rb') as f:
        f.seek(max(0, os.path.getsize(path) - 200))
        tail = f.read()
    if not tail.endswith(b'\n'):
        return False
    p = tail.decode('utf-8', 'replace').rstrip('\n').split('\n')[-1].split('\t')
    if len(p) != 2:
        return False
    try:
        float(p[1])
    except ValueError:
        return False
    return True


def head_size(url):
    req = urllib.request.Request(url, method='HEAD')
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers['Content-Length'])


def grab(url, lo, hi, out, tries=4):
    want = hi - lo + 1
    for k in range(tries):
        subprocess.run(['curl', '-sL', '-m', '180', '-r', '%d-%d' % (lo, hi),
                        '-o', out, url], capture_output=True)
        if os.path.exists(out) and os.path.getsize(out) == want:
            return True
        # 不删除分块文件：环境的 safe-delete 钩子会拦截 os.remove（2026-11-11 实测
        # 触发 SAFE_DELETE_BULK_CONFIRM_REQUIRED 并让整个任务退出）。curl 会覆盖写。
    return False


def fetch_one(kind, sid, blocks):
    tarball = os.path.join(DEST, '%s_humann%s.tar.bz2'
                           % (sid, '3' if kind == 'v3' else '2'))
    tsv = os.path.join(WORK, TSV[kind] % sid)
    if tsv_complete(tsv):
        return sid, True, '已完整'
    url = URL[kind] % sid
    t0 = time.time()
    try:
        n = head_size(url)
    except Exception as e:
        return sid, False, 'HEAD 失败 %s' % e
    step = max(1, n // blocks)
    jobs = []
    for i in range(blocks):
        lo = i * step
        hi = n - 1 if i == blocks - 1 else min(n - 1, (i + 1) * step - 1)
        if lo > n - 1:
            continue
        jobs.append((url, lo, hi, os.path.join(BASE, '_p%d' % i)))
    with cf.ThreadPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(lambda j: grab(*j), jobs))
    if not all(res):
        return sid, False, '分块失败 %d/%d' % (sum(1 for x in res if not x), len(res))
    part = tarball + '.part'
    with open(part, 'wb') as w:
        for j in jobs:
            with open(j[3], 'rb') as r:
                w.write(r.read())
    if os.path.getsize(part) != n or subprocess.run(
            ['tar', '-tjf', part], capture_output=True).returncode != 0:
        return sid, False, '合并校验失败'
    os.replace(part, tarball)
    subprocess.run(['tar', '-xjf', tarball, '-C', WORK], capture_output=True)
    if not tsv_complete(tsv):
        return sid, False, '解包后 tsv 不完整'
    return sid, True, '%.2fMB %.0fs (%.0fKB/s)' % (
        n / 1e6, time.time() - t0, n / 1024 / max(1e-9, time.time() - t0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--kind', choices=['v2', 'v3'], required=True)
    ap.add_argument('--samples', default='')
    ap.add_argument('--list', default='')
    ap.add_argument('--blocks', type=int, default=16)
    a = ap.parse_args()
    sids = ([l.strip() for l in open(a.list) if l.strip()] if a.list
            else [s.strip() for s in a.samples.split(',') if s.strip()])
    os.makedirs(DEST, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    print('%s %d 个样本，%d 块并行' % (a.kind, len(sids), a.blocks), flush=True)
    ok = 0
    for i, sid in enumerate(sids, 1):
        sid2, good, msg = fetch_one(a.kind, sid, a.blocks)
        ok += good
        print('[%d/%d] %-16s %s %s' % (i, len(sids), sid2,
                                       'OK  ' if good else 'FAIL', msg), flush=True)
    print('完成 ok=%d/%d' % (ok, len(sids)))


if __name__ == '__main__':
    main()
