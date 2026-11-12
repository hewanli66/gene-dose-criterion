# -*- coding: utf-8 -*-
"""批量下载 Pfam HMM（InterPro 接口），限速 + 重试 + 完整性校验。

坑：InterPro 对并发请求返回 429。故并发 3、失败退避重试；下载后校验文件头
含 HMMER3 标记，避免把错误页当成 HMM（pyhmmer 会静默读到脏数据）。
"""
import concurrent.futures as cf
import gzip
import os
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
HMMDIR = os.path.join(BASE, 'hmm')
URL = 'https://www.ebi.ac.uk/interpro/api/entry/pfam/%s?annotation=hmm'


def valid(path):
    """HMM 是 gzip 压缩的，直接读原始字节看不到 HMMER3 标记。
    第一版因此把全部下载误判为失败（2026-11-11）。必须解压后验。"""
    if not os.path.exists(path) or os.path.getsize(path) < 1500:
        return False
    try:
        with gzip.open(path, 'rb') as f:
            if b'HMMER3' in f.read(200):
                return True
    except (OSError, EOFError):
        pass
    try:
        with open(path, 'rb') as f:
            return b'HMMER3' in f.read(200)
    except OSError:
        return False


def one(pf):
    path = os.path.join(HMMDIR, pf + '.hmm.gz')
    if valid(path):
        return pf, True, 'cached'
    for k in range(6):
        r = subprocess.run(['curl', '-sL', '-m', '90', '-o', path, URL % pf],
                           capture_output=True)
        if valid(path):
            return pf, True, 'ok'
        time.sleep(1.5 + 2.0 * k)
    return pf, False, 'fail'


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, 'candidates_meta.txt')
    pfs = [l.strip() for l in open(src, encoding='utf-8')
           if l.strip() and not l.startswith('#')]
    os.makedirs(HMMDIR, exist_ok=True)
    print('待下载 %d 个' % len(pfs), flush=True)
    ok = fail = 0
    todo = [p for p in pfs if not valid(os.path.join(HMMDIR, p + '.hmm.gz'))]
    print('其中已缓存 %d，需下 %d' % (len(pfs) - len(todo), len(todo)), flush=True)
    with cf.ThreadPoolExecutor(3) as ex:
        for i, (pf, good, msg) in enumerate(ex.map(one, todo), 1):
            if good:
                ok += 1
            else:
                fail += 1
                print('FAIL', pf, flush=True)
            if i % 50 == 0:
                print('  %d/%d  ok=%d fail=%d' % (i, len(todo), ok, fail), flush=True)
    print('完成 ok=%d fail=%d' % (ok, fail))


if __name__ == '__main__':
    main()
