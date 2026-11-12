# -*- coding: utf-8 -*-
"""抓取 Pfam 全量家族目录（accession + 名称 + 类型），供选家族时按关键词筛选。

踩坑记录（2026-11-11）：
  1. InterPro 的 ?search= 端点在本机不可用——多词返回 204 空，单词返回 200 但结果为空。
  2. ?page=N 参数被忽略：无论 N 取何值都返回同一页（PF00001–PF00204），
     第一版脚本因此把 200 条抓了 121 遍。必须用响应里的 next（cursor 分页）。
输出 pfam_catalog.tsv：accession<TAB>type<TAB>name（按 accession 去重）
"""
import json
import os
import subprocess
import time

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, 'pfam_catalog.tsv')
START = 'https://www.ebi.ac.uk/interpro/api/entry/pfam/?page_size=200'


def get(url, retries=5):
    for k in range(retries):
        r = subprocess.run(['curl', '-sL', '-m', '120', url], capture_output=True)
        if r.returncode == 0 and r.stdout.strip():
            try:
                return json.loads(r.stdout)
            except ValueError:
                pass
        time.sleep(2 * (k + 1))
    return None


def main():
    seen, rows = set(), []
    url = START
    page = 0
    total = None
    while url:
        d = get(url)
        if d is None:
            print('page %d 失败，停止' % page, flush=True)
            break
        if total is None:
            total = d.get('count')
            print('count =', total, flush=True)
        for r in d.get('results', []):
            m = r['metadata']
            a = m['accession']
            if a in seen:
                continue
            seen.add(a)
            nm = m.get('name')
            if isinstance(nm, dict):
                nm = nm.get('name') or nm.get('short') or ''
            rows.append((a, m.get('type') or '', (nm or '').replace('\t', ' ')))
        page += 1
        if page % 25 == 0:
            print('  page %d  unique=%d' % (page, len(rows)), flush=True)
        url = d.get('next')
    with open(OUT, 'w', encoding='utf-8', newline='') as f:
        for a, t, n in rows:
            f.write('%s\t%s\t%s\n' % (a, t, n))
    print('saved %d unique rows -> %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
