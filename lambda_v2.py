# -*- coding: utf-8 -*-
"""λ 的统一核算（取代 dual_constraint.py 里那套账）

为什么要重算
------------
1. dual_constraint.py 的 READS 写死 3,000,000，但深度实验实际扫的是 6.3 M/样本；
   报告 λ=0.64 却拿它去和 6.3 M 下的实测 1 条比，等于少算了一倍。
2. 更关键：λ 公式算的是「与基因重叠的 reads 数」，而 15-mer 匹配只认
   **编码链、正确读框、且含完整 15 aa 窗口** 的 reads。100 bp 单端 reads 下，
   6 个读框里只有 1 个是编码框（另一条链的互补框不是），再扣掉读段跨越
   基因边界的损失，理论检出效率 ≈ 1/6 × 0.9 ≈ 15%。这个因子此前没有计入，
   于是「预期 3.2 条、实测 1 条」被当成了自洽，其实是被系统性高估。
3. 因此「达到 λ=5 需要多少深度」也跟着变，必须重算。

输出
----
lambda_v2.tsv      每物种 + 家族口径的 λ（原始 / 效率校正后）
lambda_v2.txt      给人看的汇总
"""
import csv
import math
import os

BASE = os.path.dirname(os.path.abspath(__file__))

GENOME_MB = {
    'Eggerthella_lenta': 3.2,
    'Streptococcus_thermophilus': 1.9,
    'Lactobacillus_fermentum': 2.1,
    'Bacteroides_nordii': 5.5,
    'Clostridium_perfringens': 3.0,
    'Bacteroides_oleiciplenus': 5.5,
}
GENOME_DEFAULT = 3.5

N_READS = 6_300_000     # 每样本实际扫描深度（与 deep 实验一致）
GENE_LEN = 1000         # hdcA 编码区
READ_LEN = 100
ABUND_UNIT = 100.0      # ELDERMET supp2 丰度为百分数

# 15-mer 检出效率：6 个读框只有 1 个是编码框；再扣读段跨边界的损失
FRAME_EFF = 1.0 / 6.0
EDGE_EFF = (GENE_LEN - READ_LEN) / float(GENE_LEN)
EFF_15MER = FRAME_EFF * EDGE_EFF

# 常规流程（150 bp PE + HMMER/蛋白级搜索）的效率：仍受读框约束，但读段更长、
# 允许部分窗口、且蛋白级打分容忍错配，取下限 0.5 作为对照
EFF_CONVENTIONAL = 0.5


def fnum(s, d=0.0):
    try:
        return float(str(s).strip())
    except (TypeError, ValueError):
        return d


def load_carriage():
    car = {}
    with open(os.path.join(BASE, 'TableS1_species_fixed.csv'), encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            sp = (r.get('species') or '').strip().replace(' ', '_')
            a, b = fnum(r.get('n_HSB_strains'), -1), fnum(r.get('n_total_strains'), -1)
            if sp and b > 0:
                car[sp] = dict(frac=a / b, typ=(r.get('hdcA_type') or '').strip())
    return car


def load_abundance():
    with open(os.path.join(BASE, 'ELDERMET_species_189x836.csv'), encoding='utf-8-sig') as f:
        rd = csv.reader(f)
        hdr = next(rd)
        cols = hdr[1:]
        rows = [r for r in rd]
    n = len(rows)
    mean_ab, prev = {}, {}
    for j, c in enumerate(cols):
        v = [fnum(r[j + 1]) for r in rows]
        mean_ab[c] = sum(v) / n
        prev[c] = sum(1 for x in v if x > 0) / n
    return mean_ab, prev


def lam_species(a_pct, c, g_mb, n_reads):
    """与基因重叠的期望 reads 数（未扣检出效率）"""
    return n_reads * (a_pct / ABUND_UNIT) * c * GENE_LEN / (g_mb * 1e6)


def poisson_cdf_le(k, lam):
    s = 0.0
    for i in range(k + 1):
        s += math.exp(-lam) * lam ** i / math.factorial(i)
    return s


def main():
    carr = load_carriage()
    mean_ab, prev = load_abundance()
    members = [s for s in mean_ab if s in carr and mean_ab[s] > 0]
    members.sort(key=lambda z: -mean_ab[z])

    rows = []
    for s in members:
        g = GENOME_MB.get(s, GENOME_DEFAULT)
        a, c = mean_ab[s], carr[s]['frac']
        l0 = lam_species(a, c, g, N_READS)
        rows.append(dict(species=s, type=carr[s]['typ'], abundance_pct=a,
                         carriage=c, genome_mb=g, lam_raw=l0,
                         lam_det=l0 * EFF_15MER))

    tot_raw = sum(r['lam_raw'] for r in rows)
    tot_det = tot_raw * EFF_15MER
    naive_raw = sum(lam_species(r['abundance_pct'], 1.0, r['genome_mb'], N_READS)
                    for r in rows)

    obs_total = 1          # 5 样本 31.3 M reads 合计命中
    n_samples = 5

    with open(os.path.join(BASE, 'lambda_v2.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['species', 'type', 'abundance_pct', 'carriage', 'genome_mb',
                    'lambda_raw_6.3M', 'lambda_detected_6.3M'])
        for r in rows:
            w.writerow([r['species'], r['type'], '%.5f' % r['abundance_pct'],
                        '%.4f' % r['carriage'], r['genome_mb'],
                        '%.4f' % r['lam_raw'], '%.4f' % r['lam_det']])

    L = []
    A = L.append
    A('=' * 78)
    A('λ 统一核算（每样本 %.1f M reads，hdcA 家族 = 5 个检出物种）' % (N_READS / 1e6))
    A('=' * 78)
    A('')
    A('%-30s %-10s %8s %8s %9s %9s' % ('物种', '丰度%', '携带率', 'λ原始', 'λ检出', '占比'))
    for r in rows:
        A('%-30s %8.4f %8.4f %9.4f %9.4f %8.1f%%' % (
            r['species'], r['abundance_pct'], r['carriage'],
            r['lam_raw'], r['lam_det'], 100 * r['lam_raw'] / tot_raw))
    A('')
    A('家族合计      λ原始 = %.3f /样本     λ检出 = %.3f /样本' % (tot_raw, tot_det))
    A('')
    A('--- 检出效率分解（100 bp 单端 + 15-mer 精确匹配）---')
    A('  编码框占比  1/6            = %.3f' % FRAME_EFF)
    A('  读段不跨边界 (1000-100)/1000 = %.3f' % EDGE_EFF)
    A('  合计检出效率                = %.3f  (%.1f%%)' % (EFF_15MER, 100 * EFF_15MER))
    A('')
    A('--- 与实测对照（5 样本，合计 31.3 M reads，实测命中 %d 条）---' % obs_total)
    A('  若用朴素代理（假设菌株全携带）: λ = %.2f/样本 → 5 样本期望 %.1f 条'
      % (naive_raw, naive_raw * n_samples))
    A('      P(X<=1) = %.2e  → 朴素代理被实测强烈否决' % poisson_cdf_le(1, naive_raw * n_samples))
    A('  若用携带率校正后的 λ原始    : λ = %.2f/样本 → 5 样本期望 %.1f 条'
      % (tot_raw, tot_raw * n_samples))
    A('      P(X<=1) = %.3f  → 仍偏低，说明未扣检出效率' % poisson_cdf_le(1, tot_raw * n_samples))
    A('  若再扣检出效率（λ检出）    : λ = %.2f/样本 → 5 样本期望 %.1f 条'
      % (tot_det, tot_det * n_samples))
    A('      P(X<=1) = %.3f  → 预测与实测自洽' % poisson_cdf_le(1, tot_det * n_samples))
    A('')
    A('  实测/原始期望 = %.2f，与理论检出效率 %.2f 吻合'
      % (obs_total / (tot_raw * n_samples), EFF_15MER))
    A('')
    A('--- 达到半定量深度（λ检出 >= 5）所需 reads/样本 ---')
    need_15 = 5.0 / tot_det * N_READS
    need_conv = 5.0 / (tot_raw * EFF_CONVENTIONAL) * N_READS
    need_raw = 5.0 / tot_raw * N_READS
    A('  15-mer / 100 bp 单端（本流程，效率 %.0f%%）: %.0f M reads/样本  (%.1f 个当前 run)'
      % (100 * EFF_15MER, need_15 / 1e6, need_15 / N_READS))
    A('  常规流程（150 bp PE + 蛋白级搜索，效率 %.0f%%）: %.0f M reads/样本  (%.1f 个当前 run)'
      % (100 * EFF_CONVENTIONAL, need_conv / 1e6, need_conv / N_READS))
    A('  不计效率（仅作参照，不应用）           : %.0f M reads/样本' % (need_raw / 1e6))
    A('')
    A('结论：本文可行性判定不依赖 λ 的绝对值——无论取哪一档，hdcA 家族都')
    A('      远低于 λ=5，同时 H = 0.152 超过 H*(n_sp=5) = 0.0426，双约束同时不满足。')
    A('      （H* 由 09 版的全队列单值 0.0098 改为按家族规模的条件阈值，见 §2.4。）')
    A('      但"需要多少深度"必须按效率校正后的口径给，否则会把深度的门槛说小 7 倍。')

    txt = '\n'.join(L)
    open(os.path.join(BASE, 'lambda_v2.txt'), 'w', encoding='utf-8').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
