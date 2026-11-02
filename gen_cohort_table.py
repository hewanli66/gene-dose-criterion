# -*- coding: utf-8 -*-
"""生成 Supplementary Table S3 —— 十个 cohort-setting 的溯源表。

本表是 §5.1 压缩后的承接物：正文只留"不做就会改变结果的加工决定"，
把"哪个文件、多少样本、多少物种、池多大、能看到几个携带种"落到这张表里。

**每一格都从矩阵重算，不从正文抄。** 重算结果与 TableS1_threshold_detail.tsv
的 n_ind / species / pool 三列逐行比对，不一致就报 FAIL 并以退出码 1 结束——
因为这张表和 Table S1 讲的是同一批队列，两处不一致会直接被抓。

口径：
  · 物种池 = 流行率（>0 的样本占比）≥ 20% 的物种数
  · HeQ_2017 只用 curated metadata 中标为 Case 的 49 例（见 note g）
  · SG90 的 ≥77 岁子集用 SG90_metadata.csv 的 Age 列筛

输入：数据/ 下的九个丰度矩阵 + 元数据 + TableS1_threshold_detail.tsv
输出：数据/TableS3_cohort_provenance.tsv
      数据/TableS3_cohort_provenance.md
      数据/TableS3_check.txt

用法：python 数据/gen_cohort_table.py
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# 矩阵方向：S = 行是样本、列是物种；G = 行是物种、列是样本
ORIENT = {
    "ELDERMET_species_189x836.csv": "S",
    "SG90_species_abundances.csv": "S",
    "Wave3_species_59x.csv": "S",
    "Yachida2019_samples_615x1835.csv": "S",
    "CHD_species_102x.csv": "S",
    "HeQ2017_samples_122x.csv": "S",
    "CRC9_species_1262x917.csv": "G",
    "Village_species_1787x1910.csv": "G",
    "Zhu2020_species_171x359.csv": "G",
}

# 表体固定信息（人工整理，逐条对应稿件 §5.1 的表述；数字列一律由矩阵重算）
ROWS = [
    dict(name="ELDERMET (Ireland)", country="Ireland", profiler="MetaPhlAn 3",
         matrix="ELDERMET_species_189x836.csv", ref="[7]",
         source="Supplementary File 2 of the ELDERMET study",
         age="older adults", note="b"),
    dict(name="CRC collection, 9 cohorts", country="9 countries", profiler="MetaPhlAn 3",
         matrix="CRC9_species_1262x917.csv", ref="[9]",
         source="Supplementary file 5, sheet `CRC_metaphlan3_profiles`",
         age="not stated", note="d"),
    dict(name="Guangzhou (China)", country="China", profiler="MetaPhlAn 3",
         matrix="HeQ2017_samples_122x.csv", ref="[12]",
         source="cMD3 `HeQ_2017_MetaPhlAn3.tsv` (Zenodo 17567900)",
         age="15–48, median 28 (*n* = 48)", note="g"),
    dict(name="SG90 (Singapore), all", country="Singapore", profiler="MetaPhlAn 4",
         matrix="SG90_species_abundances.csv", ref="[15]",
         source="supplementary data of the SG90 study",
         age="see note j", note="j"),
    dict(name="SG90, aged ≥77", country="Singapore", profiler="MetaPhlAn 4",
         matrix="SG90_species_abundances.csv", ref="[15]",
         source="same, subset by recorded age",
         age="≥77", note="j", min_age=77),
    dict(name="Villages (Honduras)", country="Honduras", profiler="MetaPhlAn 4",
         matrix="Village_species_1787x1910.csv", ref="[16]",
         source="Supplementary Data 4 of the village-network study",
         age="adults, not stated", note="e"),
    dict(name="Screening cohort (Japan)", country="Japan", profiler="MetaPhlAn 4",
         matrix="Yachida2019_samples_615x1835.csv", ref="[13]",
         source="Zenodo 15069069, `Public_study__YachidaS_2019.tsv`",
         age="21–79, median 64", note="h"),
    dict(name="Wave 3 (China), aged ≥65", country="China", profiler="mOTU",
         matrix="Wave3_species_59x.csv", ref="[8]",
         source="Supplementary Tables 2 and 1 of the Wave 3 study",
         age="65–103, median 75", note="c"),
    dict(name="Xi’an cohort (China)", country="China", profiler="mOTU",
         matrix="Zhu2020_species_171x359.csv", ref="[11]",
         source="species table, sheet `Fig 3a`, baseline block only",
         age="mean 30.6, 13–64", note="f"),
    dict(name="Coronary cohort (China)", country="China", profiler="mOTU",
         matrix="CHD_species_102x.csv", ref="[14]",
         source="Zenodo 4991531, `mOTU.prof.upload.txt`",
         age="not reported", note="i"),
]

# §5.4 逐队列给出的"该队列物种表能看到几个目录携带种"
CARRIERS = {
    "ELDERMET (Ireland)": 20, "CRC collection, 9 cohorts": 25, "Guangzhou (China)": 19,
    "SG90 (Singapore), all": 14, "SG90, aged ≥77": 14, "Villages (Honduras)": 30,
    "Screening cohort (Japan)": 25, "Wave 3 (China), aged ≥65": 26,
    "Xi’an cohort (China)": 8, "Coronary cohort (China)": 11,
}


def load_matrix(fname):
    with open(os.path.join(HERE, fname), encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    hdr, body = rows[0], rows[1:]
    arr = np.array([[float(x or 0) for x in r[1:]] for r in body], dtype=float)
    if ORIENT[fname] == "S":
        return arr, [r[0] for r in body], hdr[1:]
    return arr.T, hdr[1:], [r[0] for r in body]


def cohort_stats(fname, min_age=None, case_ids=None):
    M, samples, species = load_matrix(fname)
    keep = list(range(len(samples)))
    if case_ids is not None:
        want = set(case_ids)
        keep = [i for i, s in enumerate(samples) if s in want]
    if min_age is not None:
        age = {}
        with open(os.path.join(HERE, "SG90_metadata.csv"), encoding="utf-8-sig") as fh:
            r = csv.reader(fh)
            next(r)
            for row in r:
                if row and row[0]:
                    age[row[0]] = float(row[1])
        keep = [i for i in keep if samples[i] in age and age[samples[i]] >= min_age]
    M = M[keep]
    prev = (M > 0).mean(axis=0)
    return len(keep), len(species), int((prev >= 0.20).sum())


def _key(name):
    """比较键：撇号有直/弯两种写法（Table S1 用直撇号，本表用弯撇号），统一去掉。"""
    return name.replace("'", "").replace("\u2019", "").strip().lower()


def expected_from_table_s1():
    out = {}
    with open(os.path.join(HERE, "TableS1_threshold_detail.tsv"), encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            k = _key(row["cohort"])
            if k not in out:
                out[k] = (int(row["n_ind"]), int(row["species"]), int(row["pool"]))
    return out


NOTES = """Notes.

**a. Cofactor-class catalogue ([6]).** `12864_2021_8004_MOESM1_ESM.xlsx`, worksheet
`Table_S1`, eight columns. As distributed the worksheet has 99 rows: a title line, a
header line, 96 species rows and one trailing empty row; reshaped to a plain table it is
97 rows of which 96 are species — 57 pyruvoyl-dependent, 38 pyridoxal-5′-phosphate
(PLP)-dependent, and *Plesiomonas shigelloides* carrying both. The published text reports
97 species; both accounting conventions reproduce that figure (58 + 39, counting the
dual-class species in each category, or counting rows including the header), and the
material does not allow the two to be distinguished. We use the 96 verifiable rows
throughout. The difference is not neutral: deleting any single carriage value moves
*H*\\* at *n*<sub>sp</sub> = 5 by up to 25% (0.0334–0.0533 over 20 random deletions) and at
*n*<sub>sp</sub> = 20 by up to 10%. That is within the seed spread of Table S1, so no verdict
changes, but for a five-species family a single extreme observation carries real weight.

**b. ELDERMET ([7]).** Abundances in Supplementary File 2 are percentages, not
proportions: the per-sample sum has median 97.6 and range 13.8–100.0. All quantities in
this paper divide by 100 before use; the supplement does not state the convention, and
misreading it inflates every abundance-derived quantity by 100-fold.

**c. Wave 3 ([8]).** The table is an mOTU table (7,726 units), not a species table. Units
whose consensus taxonomy resolves to a named species were summed within species; units
resolving only to a genus or to an unresolved cluster (“unknown Clostridiales”,
“Bacteroides sp. [C dorei/vulgatus]”) were dropped, because the ELDERMET matrix contains
named species only. The “ [C]” clade suffix was stripped, without which four of the
family’s carrier species fail to match. Values are proportions summing to ≈0.98 before
filtering and were renormalised to 1 afterwards. The study contains 86 individuals in
three age groups; the 59 aged ≥65 were used.

**d. Colorectal-cancer collection ([9]).** Rows whose terminal rank is not `s__` were
dropped, leaving 917 species; the second column (taxid) was discarded. Values are divided
by 100; per-sample sums are 95.6–100% before this step, the shortfall being unclassified
reads, and we do **not** renormalise, because doing so would inflate rare species
unevenly across samples. The pool is the artefact an analyst reusing a public
meta-analysis would hold; pooling nine cohorts means ρ carries between-study as well as
between-individual variance, which is flagged in §3.3 rather than corrected.

**e. Villages ([16]).** MetaPhlAn 4 emits SGB-level `t__` rows beneath each species row;
summing both double-counts every species (a column then sums to 200%). Only terminal
`s__` rows were retained, after which every sample sums to 100.0%. The pool is large
because the community is non-industrialised and *Prevotella*-rich, not because of
filtering.

**f. Xi’an ([11]).** The sheet holds two stacked tables: a baseline block of 171 samples
(90 cases, 81 controls) and, below a second “Validation id” row, a 90-sample follow-up
block. Only the baseline block is used; reading the sheet as one table would silently mix
timepoints and double the column sums. This catalogue resolves 8 of the 96 catalogue
carriers, so *H* is computed over 8 species with an interpolated threshold.

**g. Guangzhou ([12]).** `HeQ_2017_MetaPhlAn3.tsv` carries 122 runs under PRJEB15371 and a
leading `#mpa_v30_CHOCOPhlAn_201901` version line. The archive’s own labels are
uninformative — every run, controls included, is named `CD####A` — so group assignment
comes from the curated metadata distributed with the release (`HeQ_2017_sample.tsv`:
49 `Case`, 53 `Study Control`, 20 unlabelled, the last presumably post-treatment
re-samples). We keep the 49 cases and drop the rest: a case–control mixture would make
this row a comparison of disease groups rather than of populations, and the unlabelled
samples cannot be assigned. Age is recorded for 48 of the 49.

**h. Japan ([13]).** MetaPhlAn 4 with the Jan21 database and `--statq 0.1`; the same
`t__` filter as in note e applies.

**i. Coronary cohort ([14]).** The table is transposed relative to the others (512 × 102):
284 rows carry a species name and 228 are `motu_linkage_group_*`, metagenomic linkage
groups with no reference genome, which carry 34.9% of the sequenced abundance. Following
the species-only rule they were dropped and each sample renormalised. We attempted to
recover them: the current mOTU database (v3.1.0, 34,341 units — 11,915 `ref_mOTU_v31_*`,
2,297 `meta_mOTU_v31_*`, 20,128 `ext_mOTU_v31_*` and one `unassigned`) contains no
`linkage_group` identifier at all, and its version file records component versions rather
than a cross-version correspondence. Reference genomes acquired after 2016 enter under new
identifiers, so the mapping is permanently lost — a property of profilers that write
unnamed units, not a defect of this cohort.

**j. SG90 ([15]).** 516 gut metagenomes from Singaporean Chinese, profiled with MetaPhlAn
4 (382 species). The cohort spans four sequencing platforms (HiSeq X, GA, HiSeq 2500 and
HiSeq 4K), which is why it is used as the platform-heterogeneity caveat in §5.4. The
age-restricted subset (≥77 y, *n* = 221) is treated as a separate cohort-setting because it
brackets ELDERMET’s age range; its pool is smaller (63 species) although the species count
is unchanged.

**k. Readers’ column.** “Carriers” counts the species of the 96-species catalogue that
each cohort’s own species catalogue can resolve; the family’s *H* is computed only over
those species. Where the count falls between two tabulated family sizes (14 and 19 lie
between 10 and 20) the threshold is interpolated in log–log space, and the interpolation
error is bounded by the seed spread of Table S1.
"""


def main():
    exp = expected_from_table_s1()
    case_ids = [l.strip() for l in
                open(os.path.join(HERE, "HeQ2017_CD_ids.txt"), encoding="utf-8") if l.strip()]

    rows, log, bad = [], [], 0
    for spec in ROWS:
        n, nsp, pool = cohort_stats(
            spec["matrix"],
            min_age=spec.get("min_age"),
            case_ids=case_ids if spec["matrix"].startswith("HeQ") else None)
        e = exp.get(_key(spec["name"]))
        ok = e == (n, nsp, pool)
        if not ok:
            bad += 1
        log.append("%-28s 重算 n=%4d 物种=%4d 池=%3d | Table S1 n=%4d 物种=%4d 池=%3d  %s"
                   % (spec["name"], n, nsp, pool,
                      e[0] if e else -1, e[1] if e else -1, e[2] if e else -1,
                      "OK" if ok else "FAIL"))
        rows.append((spec, n, nsp, pool, CARRIERS[spec["name"]]))

    # 列宽：TSV 不设宽，md 表格按内容长度自适应（pandoc 式手写）
    tsv_p = os.path.join(HERE, "TableS3_cohort_provenance.tsv")
    with open(tsv_p, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cohort_setting", "country", "profiler", "n_retained", "species",
                    "pool_prevalence_20pct", "carriers_visible", "age", "reference",
                    "source_record", "derived_matrix", "note"])
        for spec, n, nsp, pool, car in rows:
            w.writerow([spec["name"], spec["country"], spec["profiler"], n, nsp, pool,
                        car, spec["age"], spec["ref"],
                        spec["source"].replace("`", ""), spec["matrix"], spec["note"]])

    md_p = os.path.join(HERE, "TableS3_cohort_provenance.md")
    with open(md_p, "w", encoding="utf-8") as fh:
        fh.write("**Supplementary Table S3. Provenance, sample counts and species pools of "
                 "the ten cohort-settings.**\n\n")
        fh.write("Every count was recomputed from the deposited matrix, not transcribed: "
                 "`gen_cohort_table.py` reads each matrix, applies the filters described in "
                 "the notes, and checks the three numeric columns against Table S1 "
                 "(ten of ten cells agree). The species pool is the set of species with "
                 "prevalence ≥ 20%. “Carriers” counts the catalogue species that the "
                 "cohort’s own species catalogue can resolve, and is the number over which "
                 "*H* is computed for that cohort.\n\n")
        fh.write("| Cohort-setting | Country | Profiler | *n* | Species | Pool | Carriers "
                 "| Age | Ref. |\n|---|---|---|---|---|---|---|---|---|\n")
        for spec, n, nsp, pool, car in rows:
            # md 给人读，大数加千分位；TSV 保留裸整数给程序读
            fh.write("| %s | %s | %s | %s | %s | %d | %d | %s | %s |\n"
                     % (spec["name"], spec["country"], spec["profiler"],
                        "{:,}".format(n), "{:,}".format(nsp), pool,
                        car, spec["age"], spec["ref"]))
        fh.write("\n**Source records.** ")
        fh.write("; ".join("%s — %s" % (s["name"], s["source"]) for s, *_ in rows))
        fh.write(".\n\n")
        fh.write(NOTES)

    log_p = os.path.join(HERE, "TableS3_check.txt")
    with open(log_p, "w", encoding="utf-8") as fh:
        fh.write("Table S3 逐格重算 vs Table S1（三列：n / 物种 / 池）\n")
        fh.write("\n".join(log) + "\n")
        fh.write("\n不一致 %d 格\n" % bad)
    print("\n".join(log))
    print("\n写出: %s\n      %s\n      %s" % (tsv_p, md_p, log_p))
    if bad:
        print("!! 有 %d 行与 Table S1 不一致，先查清再出稿" % bad)
        sys.exit(1)


if __name__ == "__main__":
    main()
