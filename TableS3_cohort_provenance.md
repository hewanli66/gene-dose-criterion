**Supplementary Table S3. Provenance, sample counts and species pools of the ten cohort-settings.**

Every count was recomputed from the deposited matrix, not transcribed: `gen_cohort_table.py` reads each matrix, applies the filters described in the notes, and checks the three numeric columns against Table S1 (ten of ten cells agree). The species pool is the set of species with prevalence ≥ 20%. “Carriers” counts the catalogue species that the cohort’s own species catalogue can resolve, and is the number over which *H* is computed for that cohort.

| Cohort-setting | Country | Profiler | *n* | Species | Pool | Carriers | Age | Ref. |
|---|---|---|---|---|---|---|---|---|
| ELDERMET (Ireland) | Ireland | MetaPhlAn 3 | 189 | 836 | 149 | 20 | older adults | [7] |
| CRC collection, 9 cohorts | 9 countries | MetaPhlAn 3 | 1,262 | 917 | 164 | 25 | not stated | [9] |
| Guangzhou (China) | China | MetaPhlAn 3 | 49 | 575 | 121 | 19 | 15–48, median 28 (*n* = 48) | [12] |
| SG90 (Singapore), all | Singapore | MetaPhlAn 4 | 516 | 382 | 69 | 14 | see note j | [15] |
| SG90, aged ≥77 | Singapore | MetaPhlAn 4 | 221 | 382 | 63 | 14 | ≥77 | [15] |
| Villages (Honduras) | Honduras | MetaPhlAn 4 | 1,787 | 1,910 | 467 | 30 | adults, not stated | [16] |
| Screening cohort (Japan) | Japan | MetaPhlAn 4 | 615 | 1,835 | 300 | 25 | 21–79, median 64 | [13] |
| Wave 3 (China), aged ≥65 | China | mOTU | 59 | 3,694 | 128 | 26 | 65–103, median 75 | [8] |
| Xi’an cohort (China) | China | mOTU | 171 | 359 | 211 | 8 | mean 30.6, 13–64 | [11] |
| Coronary cohort (China) | China | mOTU | 102 | 284 | 116 | 11 | not reported | [14] |

**Source records.** ELDERMET (Ireland) — Supplementary File 2 of the ELDERMET study; CRC collection, 9 cohorts — Supplementary file 5, sheet `CRC_metaphlan3_profiles`; Guangzhou (China) — cMD3 `HeQ_2017_MetaPhlAn3.tsv` (Zenodo 17567900); SG90 (Singapore), all — supplementary data of the SG90 study; SG90, aged ≥77 — same, subset by recorded age; Villages (Honduras) — Supplementary Data 4 of the village-network study; Screening cohort (Japan) — Zenodo 15069069, `Public_study__YachidaS_2019.tsv`; Wave 3 (China), aged ≥65 — Supplementary Tables 2 and 1 of the Wave 3 study; Xi’an cohort (China) — species table, sheet `Fig 3a`, baseline block only; Coronary cohort (China) — Zenodo 4991531, `mOTU.prof.upload.txt`.

Notes.

**a. Cofactor-class catalogue ([6]).** `12864_2021_8004_MOESM1_ESM.xlsx`, worksheet
`Table_S1`, eight columns. As distributed the worksheet has 99 rows: a title line, a
header line, 96 species rows and one trailing empty row; reshaped to a plain table it is
97 rows of which 96 are species — 57 pyruvoyl-dependent, 38 pyridoxal-5′-phosphate
(PLP)-dependent, and *Plesiomonas shigelloides* carrying both. The published text reports
97 species; both accounting conventions reproduce that figure (58 + 39, counting the
dual-class species in each category, or counting rows including the header), and the
material does not allow the two to be distinguished. We use the 96 verifiable rows
throughout. The difference is not neutral: deleting any single carriage value moves
*H*\* at *n*<sub>sp</sub> = 5 by up to 25% (0.0334–0.0533 over 20 random deletions) and at
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
