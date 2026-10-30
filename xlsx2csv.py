# -*- coding: utf-8 -*-
"""
xlsx2csv.py —— 纯 stdlib 的 xlsx 工作表导出（指定 sheet 名）。

用法：
    python xlsx2csv.py <in.xlsx> <sheet_name> <out.csv>

说明：
  - 只读取单元格值与 sharedStrings，不做格式/公式求值；对本项目的补充表足够。
  - 列按绝对列号（A=0, B=1, ...）回填，空单元格写成空串，保证行列对齐。
  - 运行环境：Python 3.8+，仅标准库。
"""
import csv
import re
import sys
import zipfile
from xml.etree import ElementTree as ET

NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
ns = {'m': NS}


def col_idx(ref):
    m = re.match(r'([A-Z]+)', ref or 'A')
    c = m.group(1) if m else 'A'
    n = 0
    for ch in c:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def load(path, want):
    z = zipfile.ZipFile(path)
    ss = []
    if 'xl/sharedStrings.xml' in z.namelist():
        root = ET.fromstring(z.read('xl/sharedStrings.xml'))
        for si in root.findall('m:si', ns):
            ss.append(''.join(t.text or '' for t in si.iter('{%s}t' % NS)))

    wb = ET.fromstring(z.read('xl/workbook.xml'))
    rels = {}
    if 'xl/_rels/workbook.xml.rels' in z.namelist():
        r = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        for e in r:
            rels[e.get('Id')] = e.get('Target')
    rid_attr = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'

    target = None
    for s in wb.find('m:sheets', ns):
        if s.get('name') == want:
            t = rels.get(s.get(rid_attr), '')
            target = t if t.startswith('xl/') else 'xl/' + t.lstrip('/')
            break
    if target is None:
        raise SystemExit('sheet not found: %s' % want)

    root = ET.fromstring(z.read(target))
    rows = []
    for row in root.iter('{%s}row' % NS):
        vals = {}
        for c in row.findall('m:c', ns):
            v = c.find('m:v', ns)
            t = c.get('t')
            if t == 's' and v is not None:
                val = ss[int(v.text)]
            elif t == 'inlineStr':
                val = ''.join(x.text or '' for x in c.iter('{%s}t' % NS))
            else:
                val = v.text if v is not None else ''
            vals[col_idx(c.get('r'))] = val
        if vals:
            mx = max(vals)
            rows.append([vals.get(i, '') for i in range(mx + 1)])
    width = max(len(r) for r in rows) if rows else 0
    return [r + [''] * (width - len(r)) for r in rows]


def main():
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    rows = load(sys.argv[1], sys.argv[2])
    with open(sys.argv[3], 'w', newline='', encoding='utf-8-sig') as f:
        csv.writer(f).writerows(rows)
    print('rows=%d cols=%d -> %s' % (len(rows), len(rows[0]) if rows else 0, sys.argv[3]))


if __name__ == '__main__':
    main()
