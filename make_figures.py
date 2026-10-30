# -*- coding: utf-8 -*-
"""
make_figures.py —— 文章 ① 图 1–3（纯 SVG，期刊级矢量）
运行：python make_figures.py
输出：../图/fig1_dual_constraint.svg 等
"""
import os, csv, math, json
from svgplot import Fig, Ax, C_BLUE, C_ORANGE, C_GREEN, C_RED, C_PURPLE, C_SKY, C_GRAY, C_DARK

HERE = os.path.dirname(os.path.abspath(__file__))
# 输出目录可由环境变量 FIG_OUT 指定，用于在不同稿件版本间出图而不互相覆盖
OUT = os.path.normpath(os.environ.get("FIG_OUT") or os.path.join(HERE, "..", "图"))
os.makedirs(OUT, exist_ok=True)

# ---------------- 数据 ----------------
# 阈值（全队列，条件概率口径）：review_final.py §1
# 旧值 0.00031/0.00975/0.03610/0.06252 系"H 的 30% 分位数"，已废弃
HSTAR = {3: 0.03430, 5: 0.04261, 10: 0.07685, 20: 0.09230}
# 环境变量 HSTAR_JSON 指向 {"3":x,"5":x,"10":x,"20":x} 时改用该组阈值，
# 供"中位数版"稿件出图；不设该变量则用上面的单种子值，行为与旧版一致。
_hj = os.environ.get("HSTAR_JSON")
if _hj and os.path.exists(_hj):
    with open(_hj, encoding="utf-8") as _f:
        HSTAR = {int(k): float(v) for k, v in json.load(_f).items()}
    print("[make_figures] HSTAR <- %s : %s" % (_hj, HSTAR))
# 累积口径（仅作对照，见图 1a 虚线）
HSTAR_CUM = {3: 0.1162, 5: 0.1273, 10: 0.1185, 20: 0.1123}
# 双约束判决：dual_constraint.tsv
FAM = {
    "hdcA_all":  dict(label="hdcA family (n=20)", n=20, c_eff=0.278, H=0.15217, lam=1.083),
    "hdcA_contrib": dict(label="hdcA (contributing species, n=5)", n=5, c_eff=0.278,
                         H=0.15217, lam=1.083),
    "hdcA_hi":   dict(label="hdcA \u2265 40% (n=3)", n=3, c_eff=0.796, H=0.08078, lam=0.937),
    "elenta":    dict(label="E. lenta", n=1, c_eff=1.0, H=0.0, lam=0.102),
    "cperf":     dict(label="C. perfringens", n=1, c_eff=1.0, H=0.0, lam=0.739),
    "sthermo":   dict(label="S. thermophilus", n=1, c_eff=0.052, H=0.0, lam=0.117),
    "lferm":     dict(label="L. fermentum", n=1, c_eff=0.016, H=0.0, lam=0.029),
}
LAM_MIN = 5.0

rows = list(csv.DictReader(open(os.path.join(HERE, "TableS1_species_fixed.csv"), encoding="utf-8-sig")))
SPEC = []
for r in rows:
    k, n = int(r["n_HSB_strains"]), int(r["n_total_strains"])
    SPEC.append(dict(sp=r["species"], typ=r["hdcA_type"], cr=k / n, n=n))


# =========================================================
# 图 1 —— 双约束判决图
# =========================================================
def fig1():
    f = Fig(510, 330)
    ax = Ax(f, (62, 46, 232, 246), (2, 26), (0.025, 0.16))

    # --- Panel A：H* 随家族规模变化 ---
    f.text(40, 34, "a", size=11, weight="bold", anchor="middle")
    ax.grid(ys=[0.03, 0.05, 0.1])
    ax.frame()
    ax.xticks([3, 5, 10, 20], ["3", "5", "10", "20"])
    ax.yticks([0.03, 0.05, 0.1], ["0.03", "0.05", "0.10"])
    ax.xlabel("Species per family (n)", dy=16)
    ax.ylabel("Threshold  H*", dx=30)
    f.text(147, 26, "Conditional threshold for \u03c1 \u2265 0.90", size=6.5, fill="#555555")

    xs = [3, 5, 10, 20]
    ys = [HSTAR[n] for n in xs]

    # 累积口径（对照，虚线）——说明两种定义差别有多大
    yc = [HSTAR_CUM[n] for n in xs]
    f.polyline([(ax.X(x), ax.Y(y)) for x, y in zip(xs, yc)],
               stroke="#B9B9B9", sw=0.9, dash="3,2.5")
    f.text(ax.X(20) - 2, ax.Y(HSTAR_CUM[20]) + 10, "cumulative definition", size=5.6,
           fill="#8A8A8A", anchor="end", halo=True)

    f.polyline([(ax.X(x), ax.Y(y)) for x, y in zip(xs, ys)], stroke=C_BLUE, sw=1.4)
    # FIG_SPAN 决定标注报哪一段的倍数（默认 3-20，与单种子版一致）；
    # FIG_WEAK3=1 时把 n=3 点画成空心，提示该端点弱识别（中位数版）
    _span = os.environ.get("FIG_SPAN") or "3-20"
    _k0, _k1 = (int(t) for t in _span.split("-"))
    _weak3 = os.environ.get("FIG_WEAK3") == "1"
    for x, y in zip(xs, ys):
        if _weak3 and x == 3:
            f.circle(ax.X(x), ax.Y(y), 3.2, fill="#ffffff", stroke=C_BLUE, sw=1.2)
        else:
            f.circle(ax.X(x), ax.Y(y), 3.2, fill=C_BLUE, stroke="#ffffff", sw=1.0)
        f.text(ax.X(x), ax.Y(y) + 11, "%.3f" % y, size=5.8, anchor="middle",
               fill="#555555", halo=True)
    # 倍数注释（放在曲线左下方的空白处，避开累积线与曲线）
    fold = HSTAR[_k1] / HSTAR[_k0]
    f.text(72, 110, "H* spans %.1f-fold" % fold, size=6.5, fill=C_RED, anchor="start", halo=True)
    f.text(72, 119, "from n = %d to n = %d" % (_k0, _k1), size=6.5, fill=C_RED, anchor="start", halo=True)
    subset_fold = os.environ.get("FIG_SUBSET_FOLD", "1.4")
    f.text(72, 128, "within %s-fold across subsets" % subset_fold, size=6, fill="#666666", anchor="start", halo=True)

    # --- Panel B：归一化判决图 ---
    ax2 = Ax(f, (312, 46, 484, 246), (0.004, 200), (0.1, 1000))
    f.text(300, 34, "b", size=11, weight="bold", anchor="middle")

    # 四象限底色
    xm, ym = ax2.X(LAM_MIN), ax2.Y(1.0)
    ax2.f.rect(ax2.x0, ax2.y0, xm - ax2.x0, ym - ax2.y0, fill="#FDF0EC")
    ax2.f.rect(xm, ax2.y0, ax2.x1 - xm, ym - ax2.y0, fill="#FBF3E2")
    ax2.f.rect(ax2.x0, ym, xm - ax2.x0, ax2.y1 - ym, fill="#EAF3F8")
    ax2.f.rect(xm, ym, ax2.x1 - xm, ax2.y1 - ym, fill="#EAF6EF")

    ax2.frame()
    ax2.xticks([0.01, 0.1, 1, 5, 10, 100], ["0.01", "0.1", "1", "5", "10", "100"])
    ax2.yticks([0.1, 1, 10, 100, 1000], ["0.1", "1", "10", "100", "1000"])
    ax2.xlabel("Expected detectable hits per sample  \u03bbdet  (log)", dy=16)
    ax2.ylabel("Normalised dispersion  H / H*(n)", dx=36)

    # 判据线
    f.line(ax2.x0, ym, ax2.x1, ym, stroke=C_DARK, w=1.2, dash="4,2.5")
    f.line(xm, ax2.y0, xm, ax2.y1, stroke=C_DARK, w=1.2, dash="4,2.5")
    f.text(ax2.x0 + 5, ym + 13, "H = H*", size=6.5, anchor="start", fill=C_DARK, halo=True)
    f.text(xm + 3, ax2.y0 + 8, "\u03bb = 5", size=6.5, anchor="start", fill=C_DARK)

    # 象限标签
    f.text(ax2.x0 + 6, ax2.Y(25), "not feasible", size=6.5, fill="#9A5B45", weight="bold", anchor="start")
    f.text(ax2.x0 + 6, ax2.Y(12), "at current depth", size=6, fill="#9A5B45", anchor="start")
    f.text(ax2.x0 + 6, ax2.Y(0.22), "proxy unbiased,", size=6.5, fill="#3A6F8A", weight="bold", anchor="start")
    f.text(ax2.x0 + 6, ax2.Y(0.165), "but out of reach", size=6, fill="#3A6F8A", anchor="start")
    f.text(ax2.x1 - 5, ax2.Y(25), "gene-level", size=6.5, fill="#96702F", weight="bold", anchor="end")
    f.text(ax2.x1 - 5, ax2.Y(0.20), "proxy usable", size=6.5, fill="#2F7A55", weight="bold", anchor="end")

    # hdcA 的两种家族口径：n = 20（全部匹配种）与 n = 5（仅贡献非零丰度者，正文主口径）
    # λ 相同、H 相同，只有阈值随规模变，所以两点竖直排列——这正是要展示的：
    # 口径改变只移动 H/H*，不移动 λ，判定不变。
    for key in ("hdcA_all", "hdcA_contrib"):
        d = FAM[key]
        r = d["H"] / HSTAR[d["n"]]
        f.circle(ax2.X(d["lam"]), ax2.Y(r), 3.4, fill=C_RED, stroke="#ffffff", sw=1.0)
    x_a = ax2.X(FAM["hdcA_all"]["lam"])
    ya = ax2.Y(FAM["hdcA_all"]["H"] / HSTAR[20])
    yb = ax2.Y(FAM["hdcA_contrib"]["H"] / HSTAR[5])
    f.line(x_a, ya - 3.2, x_a, yb + 3.2, stroke=C_RED, w=1.0, dash="2,1.6", opacity=0.55)
    # 标签右对齐贴在 λ = 5 虚线左侧：原先左对齐从 ax2.x0+95 起排，文字会压在这条虚线上
    f.text(xm - 4, 156, "hdcA family", size=5.9, anchor="end", fill=C_DARK, halo=True)
    # 两个口径的比值由 HSTAR 现算，不手写——否则换阈值组（如中位数版）时标注会与点位不符
    _r20 = FAM["hdcA_all"]["H"] / HSTAR[20]
    _r5 = FAM["hdcA_contrib"]["H"] / HSTAR[5]
    f.text(xm - 4, 148, "n=20: %.1f   n=5: %.1f" % (_r20, _r5), size=5.2, anchor="end",
           fill="#555555", halo=True)

    # 加深测序的横向箭头
    y_a = ax2.Y(FAM["hdcA_all"]["H"] / HSTAR[20])
    f.arrow(ax2.X(0.72), y_a, xm - 1, y_a, stroke="#AFAFAF", dash="2,2", head=3.0)
    f.text(xm - 2, y_a + 13, "\u00d7 4.6 depth \u2192 29 M reads",
           size=5.4, fill="#8A8A8A", halo=True, anchor="end")

    # --- 底部条带：n = 1 家族（H 恒为 0，在对数轴上无位置）---
    sy = 286
    f.line(ax2.x0, sy, ax2.x1, sy, stroke="#C8C8C8", w=0.9)
    f.line(xm, 252, xm, sy + 7, stroke="#B0B0B0", w=0.8, dash="3,2")
    for key, dx, anc in (("elenta", -5, "end"), ("cperf", 5, "start")):
        d = FAM[key]
        f.circle(ax2.X(d["lam"]), sy, 2.8, fill=C_BLUE, stroke="#ffffff", sw=0.8)
        f.text(ax2.X(d["lam"]) + dx, sy - 9, d["label"], size=5.8, anchor=anc, fill=C_DARK)
        f.text(ax2.X(d["lam"]) + dx, sy + 10, "\u03bb = %.3f" % d["lam"], size=5.4,
               anchor=anc, fill="#777777")
    f.text(ax2.x0, 312, "n = 1 families: H \u2261 0, so only \u03bb can fail",
           size=5.8, anchor="start", fill="#555555")

    p = os.path.join(OUT, "fig1_dual_constraint.svg")
    f.save(p)
    return p


# =========================================================
# 图 2 —— 携带率双峰分布
# =========================================================
def fig2():
    f = Fig(430, 272)
    ax = Ax(f, (54, 52, 400, 230), (0, 1.0), (0, 74), xlog=False, ylog=False)

    NB = 20
    bins = [0.0] * NB
    byp = [0.0] * NB
    for s in SPEC:
        i = min(int(round(s["cr"] * NB)), NB - 1) if s["cr"] < 1 else NB - 1
        bins[i] += 1
        if s["typ"].startswith("pyruvoyl"):
            byp[i] += 1

    ax.grid(ys=[0, 20, 40, 60])
    ax.frame()
    ax.yticks([0, 20, 40, 60], ["0", "20", "40", "60"])
    ax.xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0", "0.2", "0.4", "0.6", "0.8", "1.0"])
    ax.xlabel("Within-species carriage rate of the PDD gene family", dy=17)
    ax.ylabel("Number of species", dx=28)
    f.text(227, 24, "Annotated catalogue: 96 species", size=7, fill="#555555")

    bw = ax.X(1.0 / NB) - ax.X(0)
    for i in range(NB):
        x = ax.X(i / NB)
        h = ax.Y(0) - ax.Y(byp[i]) if byp[i] else 0
        f.rect(x + 0.6, ax.Y(byp[i]) if byp[i] else ax.Y(0), bw - 1.2, h,
               fill=C_BLUE, stroke="none")
        oth = bins[i] - byp[i]
        if oth:
            f.rect(x + 0.6, ax.Y(bins[i]), bw - 1.2, ax.Y(byp[i]) - ax.Y(bins[i]),
                   fill=C_ORANGE, stroke="none")

    # 顶部 rug
    for typ, yy, col in [("pyruvoyl", 36, C_BLUE), ("PLP", 42, C_ORANGE)]:
        sel = [s for s in SPEC if s["typ"] == typ or (typ == "pyruvoyl" and s["typ"] == "pyruvoyl/PLP")]
        for j, s in enumerate(sel):
            f.line(ax.X(s["cr"]), yy, ax.X(s["cr"]), yy + 5, stroke=col, w=0.5, opacity=0.55)

    # 标注
    f.text(ax.X(0.94), 56, "61/96 (64%)", size=7, weight="bold", fill=C_BLUE, anchor="end")
    f.text(ax.X(0.94), 64, "fully carried", size=6, fill="#555555", anchor="end")
    f.line(ax.X(0.02), ax.Y(14), ax.X(0.18), ax.Y(14), stroke=C_RED, w=1.0)
    f.line(ax.X(0.02), ax.Y(14) - 3, ax.X(0.02), ax.Y(14) + 3, stroke=C_RED, w=1.0)
    f.line(ax.X(0.18), ax.Y(14) - 3, ax.X(0.18), ax.Y(14) + 3, stroke=C_RED, w=1.0)
    f.text(ax.X(0.19), ax.Y(14), "12/96 (12%) below 20%", size=6.2, anchor="start", fill=C_RED)

    # 图例
    f.rect(60, 60, 7, 7, fill=C_BLUE)
    f.text(70, 63.5, "pyruvoyl", size=6.2, anchor="start")
    f.rect(122, 60, 7, 7, fill=C_ORANGE)
    f.text(132, 63.5, "PLP-dependent", size=6.2, anchor="start")
    f.text(60, 76, "Median carriage = 1.00; range 0.2\u2013100%", size=5.8, anchor="start", fill="#666666")

    p = os.path.join(OUT, "fig2_carriage_bimodal.svg")
    f.save(p)
    return p


# =========================================================
# 图 3 —— 三层漏斗
# =========================================================
def fig3():
    f = Fig(510, 300)
    f.text(40, 30, "a", size=11, weight="bold", anchor="middle")
    f.text(300, 30, "b", size=11, weight="bold", anchor="middle")

    # ---------------- Panel a：水平漏斗 ----------------
    stages = [(96, "Annotated catalogue", C_GRAY, "100%"),
              (20, "Matched in ELDERMET", C_SKY, "21%"),
              (5, "Detected in cohort", C_BLUE, "25%"),
              (3, "Carriage \u2265 40%", C_GREEN, "60%")]
    X0, SC = 150, 1.35
    yc = [66, 106, 146, 186]
    hgt = 24
    f.text(288, 44, "retained", size=5.8, fill="#888888", anchor="start")

    # 漏斗底色（连接各条右端）
    poly = [(X0, yc[0] + hgt / 2)]
    for k, (cnt, _, _, _) in enumerate(stages):
        poly.append((X0 + cnt * SC, yc[k] + hgt / 2))
    for k in range(len(stages) - 1, -1, -1):
        poly.append((X0, yc[k] - hgt / 2))
    f.poly(poly, fill="#F4F4F4", stroke="none")

    for k, (cnt, lab, col, ret) in enumerate(stages):
        w = cnt * SC
        y = yc[k] - hgt / 2
        f.text(X0 - 8, yc[k], lab, size=6.6, anchor="end", fill="#444444")
        f.rect(X0, y, w, hgt, fill=col, stroke="none")
        if w >= 22:
            f.text(X0 + w / 2, yc[k], str(cnt), size=9, weight="bold", fill="#ffffff")
        else:
            f.text(X0 + w + 4, yc[k], str(cnt), size=9, weight="bold", anchor="start", fill=col)
        f.text(288, yc[k], ret, size=6.4, anchor="start", fill="#777777")
    f.text(X0, 214, "overall attrition 32\u00d7 (96 \u2192 3)", size=6.6, fill=C_RED, weight="bold", anchor="start")

    # ---------------- Panel b：按辅基类型分组柱 ----------------
    bx0, by0, bx1, by1 = 330, 118, 496, 210
    axb = Ax(f, (bx0, by0, bx1, by1), (0, 3), (0, 68), xlog=False, ylog=False)
    axb.grid(ys=[0, 20, 40, 60])
    axb.frame()
    axb.yticks([0, 20, 40, 60], ["0", "20", "40", "60"], size=6)
    axb.ylabel("Species", dx=16, size=7)

    groups = [("catalogue", 57, 38), ("matched", 15, 4), ("detected", 5, 0)]
    for gi, (name, vp, vl) in enumerate(groups):
        cx = axb.X(gi + 0.5)
        f.text(cx, by1 + 10, name, size=5.8, fill="#444444")
        for off, v, col in [(-9, vp, C_BLUE), (9, vl, C_ORANGE)]:
            if v > 0:
                f.rect(cx + off - 7, axb.Y(v), 14, axb.Y(0) - axb.Y(v), fill=col, stroke="none")
                f.text(cx + off, axb.Y(v) - 6, str(v), size=6.4, fill=col)
            else:
                f.text(cx + off, axb.Y(0) - 7, "0", size=7, fill=C_ORANGE, weight="bold")
                f.line(cx + off - 3, axb.Y(0) - 21, cx + off + 3, axb.Y(0) - 15, stroke=C_ORANGE, w=1.1)
                f.line(cx + off + 3, axb.Y(0) - 21, cx + off - 3, axb.Y(0) - 15, stroke=C_ORANGE, w=1.1)

    # 图例
    f.rect(332, 104, 6, 6, fill=C_BLUE)
    f.text(341, 107, "pyruvoyl", size=5.8, anchor="start")
    f.rect(382, 104, 6, 6, fill=C_ORANGE)
    f.text(391, 107, "PLP-dependent", size=5.8, anchor="start")

    f.text(300, 244, "PLP: 38 \u2192 4 \u2192 0", size=6.4, anchor="start", fill="#8A5A20", weight="bold")
    f.text(300, 254, "the 4 matched species have mean abundance = 0 across all 189 samples",
           size=5.8, anchor="start", fill="#8A5A20")
    f.text(300, 263, "\u2192 absence or near-absence below the detection floor (\u00a72.1)",
           size=5.8, anchor="start", fill="#8A5A20")

    p = os.path.join(OUT, "fig3_funnel.svg")
    f.save(p)
    return p


if __name__ == "__main__":
    for fn in (fig1, fig2, fig3):
        print("saved:", fn())
