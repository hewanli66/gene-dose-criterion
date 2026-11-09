# Reproducibility package

**Species abundance is not a gene dose: a two-constraint criterion for functional
inference in metagenomics**

This package contains the code and the derived tables behind every number, table and
figure in the manuscript. It is self-contained in one respect and deliberately not in
another:

- **Self-contained** for everything that is a *result*: the per-cell thresholds, the
  seed spreads, the sensitivity grid and the two supplementary tables can all be
  regenerated from the files deposited here, without re-running a single simulation.
- **Not self-contained** for the *primary abundance matrices*, which belong to the
  studies that produced them. They are not redistributed here. `data_sources.tsv`
  lists, file by file, which supplement or archive each one comes from.

No third-party Python packages are required. Everything is standard library; the SVG
plotting and the xlsx reader are the two small modules included at the bottom of the
list below.

---

## 1. Quick check

Two commands run against what is in this folder and take a few seconds:

```bash
python gen_supp_tables.py     # rebuilds Supplementary Tables S1 and S2,
                              # and cross-checks all 40 cells against Table 5
python gen_cohort_table.py    # rebuilds Supplementary Table S3 (cohort provenance)
                              # and checks its n / species / pool against Table S1
python diag_ratio_spread.py   # recomputes the seed spreads quoted in section 2.4
python gen_bsh_table.py       # rebuilds Supplementary Table S4 from domain_carriage.tsv
```

`gen_supp_tables.py` prints `Table S1 与 Table 5 逐格核对：40 格，0 处不一致`
(“40 cells, 0 mismatches”) as its first line. That line is the whole point of the
check: the supplementary table and the main table are generated from different
deposited files and must agree cell by cell.

---

## 2. What each file is

### Scripts

| File | Role | Reads | Writes |
|---|---|---|---|
| `review_final.py` | Shared library. Defines the simulated family ensemble (`families`), the conditional threshold `h_cond` — the largest *H* at which P(ρ ≥ 0.90 \| H ≈ h) ≥ 0.95 — and Spearman ρ. Its `__main__` reproduces section 2.7 and Table 4, taking the conditional thresholds from `hstar_median.json` so that its ratios match Table 4. | `TableS1_species_fixed.csv`, `ELDERMET_species_189x836.csv`, `hstar_median.json` | `review_M11_grid.tsv`, `final_numbers.txt` (local), console |
| `cross_cohort9.py` | Cohort-agnostic matrix loader, log–log interpolation of *H*\*, carriage matching for *hdcA*, the prevalence ≥ 20% pool rule. | species matrices | `cross_cohort9.tsv` |
| `cross_cohort9_seed.py` | **Table 5 and Table 6.** Eight seeds per cell across the ten cohort-settings, 40 cells in Table 5 and 36 in Table 6. | the ten matrices | `cross_cohort9_seed_t5.tsv`, `cross_cohort9_seed_t5_detail.txt`, `cross_cohort9_seed_t6.tsv` |
| `heq_cd_only.py` | Replaces the Guangzhou row with the 49 samples the curated metadata labels as cases (§2.1.6d). | `HeQ2017_samples_122x.csv`, `HeQ2017_CD_ids.txt` | `heq_cd_only.txt`, `cross_cohort9_seed_t5_heqcd.tsv` |
| `gen_supp_tables.py` | **Supplementary Tables S1 and S2**, plus the 40-cell cross-check against Table 5. | the four seed/grid result files | `TableS1_threshold_detail.{tsv,md}`, `TableS2_M11_grid.{tsv,md}`, `supp_tables_check.txt` |
| `diag_ratio_spread.py` | Per-size seed spreads and the spread of the 5→20 ratio, as quoted in section 2.4. | the two seed-detail files | console |
| `diag_carriage_frac.py` | Separates *c*<sub>eff</sub> (abundance-weighted mean carriage) from the fraction of individuals carrying at least one carrier species. The two are quoted side by side in section 2.4 and are not the same quantity. | three matrices | console |
| `diag_m11_grid.py` | The 3 × 3 sensitivity grid over fidelity cut-off and retention level. | ELDERMET matrix | `review_M11_grid_median.tsv`, `diag_m11_grid.txt` |
| `diag_window_width.py` | Width of the moving window (0.01 / 0.02 / 0.03) at *n*<sub>sp</sub> = 3 and 20. | ELDERMET matrix | `diag_window_width.txt` |
| `supplement_hstar.py` | The three ELDERMET residence subsets (all / community-dwelling / care). | matrices, `ELDERMET_clinical_189.csv` | `supplement_hstar.txt` |
| `make_figures.py` | Figures 1–6 as SVG. | `TableS1_species_fixed.csv`, `hstar_median.json` | SVG files in `$FIG_OUT` |
| `mpa2matrix.py`, `motu2matrix.py`, `motu_chd2matrix.py` | Convert the published supplementary tables (MetaPhlAn clade table, mOTU table, transposed CHD table) into the sample × species CSV the analysis reads. | third-party tables | the CSV inputs listed in `data_sources.tsv` |
| `cross_cohort.py`, `cross_cohort3.py`, `cross_cohort5.py`, `cross_cohort8.py`, `cross_cohort_seed.py` | The same computation at each earlier stage of the cohort set (one, three, five, eight cohorts, and the four-seed variant). Kept because Methods §2.4 traces the sequence of additions through them; `cross_cohort9_seed.py` supersedes all of them for the values reported in the paper. | species matrices | `cross_cohort*.tsv`, `cross_cohort*.txt` |
| `eff_calibration.py`, `lambda_v2.py`, `spikein_test.py` | Section 2.7: detector efficiency from an endogenous single-copy control (`rpoB`), the Poisson expectation for a 1-kb gene, and the in-silico spike-in. **ε is calibrated, not derived:** `lambda_v2.py` takes ε = 0.71 from `eff_calibration.py` (range 0.48–0.94) and no longer uses the retired 1/6 reading-frame argument, which applied only to a single fixed frame and is contradicted by the measured recall of 1.000. Genome sizes are those of the finished type strains where one exists — *C. perfringens* ATCC 13124 = 3,256,683 bp (CP000246). The script prints both sensitivity grids (ε, and *G* for *C. perfringens*). | reads, `rpob_ref.json`, `conserved_kmers.json` | `eff_calibration.txt`, `lambda_v2.tsv`, `lambda_v2.txt` |
| `fetch_rpob.py`, `extract_kmers.py`, `scan_rpob.py` | Fetch the `rpoB` references from the EBI Proteins API, build the conserved k-mer set and the full k-mer index, and scan reads against it. | EBI Proteins API, ENA reads | `rpob_ref.json`, `conserved_kmers.json`, `reference_kmers.json` |
| `domain_carriage.py` | **Supplementary Table S4**: strain-level carriage of bile salt hydrolase, measured genome by genome. Downloads proteomes from **UniParc** (not UniProtKB — see note below), scans them with a Pfam HMM through pyhmmer, and computes *c*<sub>*s*</sub>, *H* and λ. | UniParc/EBI and InterPro (network); local Pfam HMM files | `domain_family.tsv`, `domain_carriage.tsv`, `domain_hits.tsv` |
| `gen_bsh_table.py` | Lays out Supplementary Table S4 with the abundance weights and the per-species contribution to *H*, and re-derives *H* from those columns as a check. It also runs the single-species sensitivity test of section 5.10: setting a species to zero carriage removes a carrier, so *n*<sub>sp</sub> falls from 11 to 10 and the ratios are taken against *H*\*(10) = 0.0732 rather than *H*\*(11) = 0.0757. | `domain_carriage.tsv` | `TableS4_bsh_carriage.{tsv,md}` |
| `m6_discriminant.py` | Section 2.10 and Supplementary Table S5: whether *H* carries information beyond *c*<sub>eff</sub>. Enumerates sub-families of each cohort’s real *hdcA* carrier set (sizes 3, 4, 5, 6, 8, 10; up to 250 subsets per cohort and size), scores each against ρ on the real abundance profiles, and compares the two rules by AUC, confusion matrix and within-quartile ranking. Also builds the analytic families with fixed carriage, the synthetic sample from the calibration model, and the catalogue-band families. Seed 20261103. | species matrices, `cross_cohort9.tsv`, `TableS1_species_fixed.csv` | `m6_ensemble.tsv`, `m6_report.txt` (console) |
| `svgplot.py` | Minimal SVG plotting layer used by `make_figures.py`. No matplotlib. | — | — |
| `xlsx2csv.py` | Minimal xlsx reader used by `motu2matrix.py`. | — | — |

> **Why UniParc and not UniProtKB.** A proteome-indexed FASTA request against UniProtKB
> returns an empty record for every strain except the species reference, even though the
> proteome is listed with a full protein count. UniParc holds the sequences for all
> strains and indexes them by proteome, so `uniparc/stream?query=proteome:<UPID>` is the
> endpoint that works; the returned sequence count matches the listed protein count. This
> cost half an hour of empty downloads to find, and is recorded so that it does not cost
> the next reader the same. No functional annotation is needed for the HMM scan, so
> UniParc’s lack of annotation is not a limitation.

### Deposited results

| File | What it holds |
|---|---|
| `gen_cohort_table.py` | **Supplementary Table S3**: the ten cohort-settings with samples retained, species count, prevalence-filtered pool and catalogue carriers visible. | the nine abundance matrices in this folder | `TableS3_cohort_provenance.{tsv,md}`, `TableS3_check.txt` |
| `TableS3_cohort_provenance.tsv` | **Supplementary Table S3**, one row per cohort-setting, with the depository record, the file used and the filtering note. | | |
| `TableS1_threshold_detail.tsv` | **Supplementary Table S1**, long format: 40 rows (ten settings × four family sizes), with median, eight-seed min and max, and the range as a percentage of the median. |
| `TableS1_threshold_detail.md` | The same table laid out for reading, with footnotes. |
| `TableS2_M11_grid.tsv` | **Supplementary Table S2**, nine rows: fidelity cut-off × retention level, median threshold, seed range, and the *hdcA* margin on the median and on the least favourable seed. |
| `TableS2_M11_grid.md` | The same grid as a 3 × 3 panel. |
| `TableS4_bsh_carriage.tsv` | **Supplementary Table S4** (Additional file 4), one row per species: carrying and tested genomes, *c*<sub>*s*</sub>, mean relative abundance, renormalised weight, and contribution to *H*. The contribution column sums to *H*. |
| `TableS4_bsh_carriage.md` | The same table laid out for reading, with the *H* = 0.0364, *c*<sub>eff</sub> = 0.936 summary in the footnote. |
| `TableS4_domain_family.tsv` | The family-level line behind S4: *n*<sub>sp</sub>, *c*<sub>eff</sub>, *H*, *H*\*, the ratio and λ. |
| `domain_carriage.py` | The pipeline that produced the S4 numbers; kept here so the catalogue can be rebuilt rather than trusted. |
| `supp_tables_check.txt` | The cell-by-cell comparison of Table S1 against Table 5, and the pooled ranges. |
| `cross_cohort9_seed_t5_heqcd.tsv` | **Table 5**, one row per cohort-setting. |
| `cross_cohort9_seed_t6.tsv` | **Table 6**, ten rows matched to *n* = 59, twelve draws each. |
| `cross_cohort9_seed_t5_detail.txt` | All eight seed values behind every cell of Table 5, for nine of the ten settings. |
| `heq_cd_only.txt` | The same detail for the Guangzhou 49-case row, which replaces the 122-sample version in the file above. |
| `review_M11_grid_median.tsv` | The eight-seed medians of the 3 × 3 grid. |
| `diag_window_width.txt`, `diag_m11_grid.txt` | Console output of the two diagnostics, kept so that the quoted numbers have a traceable source. |
| `hstar_median.json` | The ELDERMET threshold at the four family sizes, as used by the figure script. |
| `TableS1_species_fixed.csv` | The carriage catalogue: 96 species, eight columns, strain-level *hdcA* carriage and cofactor class. This is a derived table, not primary data, and it is included so that the analysis can be re-run without re-parsing the publisher’s spreadsheet. |
| `HeQ2017_CD_ids.txt` | The 49 Guangzhou sample IDs labelled as cases. |
| `_ena_heq_runs.tsv` | The 122 run accessions of PRJEB15371, for checking against the 117 samples the source paper reports (§ reference 12). |
| `cross_cohort.tsv`, `cross_cohort3.tsv`, `cross_cohort5.tsv`, `cross_cohort8.tsv`, `cross_cohort_seed.tsv`, `cross_cohort9.tsv` | The threshold table at each stage of the cohort set (1, 3, 5, 8, 9 cohorts, plus the four-seed variant). Section 2.4 quotes the progression; `cross_cohort9_seed_t5_heqcd.tsv` is the final one. |
| `eff_calibration.txt`, `lambda_v2.tsv`, `lambda_v2.txt` | Console output behind section 2.7 — the efficiency estimate, the observed/expected ratios and the length-normalised comparison, with the ε and genome-size sensitivity grids. |
| `TableS5_m6_ensemble.tsv` | Supplementary Table S5: all 17,043 families behind section 2.10 — the 10,104 real sub-families, 4,320 synthetic, 2,430 catalogue-band and 189 analytic families, each with *c*<sub>eff</sub>, *H*, the measured ρ and the threshold applied. Re-score any subset by filtering on the `kind` column. |
| `supplement_hstar.txt` | The three ELDERMET residence subsets (all 189, community-dwelling 105, care 53). |
| `rpob_ref.json`, `conserved_kmers.json`, `reference_kmers.json` | `rpoB` reference sequences, the conserved k-mer set used for the read-level scan, and the full k-mer index (0.6 MB) so that the scan can be re-run without refetching. |

---

## 3. Reproducing each item from the primary data

Place the matrices named in `data_sources.tsv` in this folder, keeping the file names
as given; every script resolves its inputs relative to its own location. Missing
matrices are reported as `[缺文件]` and skipped rather than crashing, so partial runs
are possible.

| Manuscript item | Command | Cost |
|---|---|---|
| Supplementary Tables S1, S2 | `python gen_supp_tables.py` | seconds |
| Supplementary Table S3 | `python gen_cohort_table.py` | seconds |
| Section 2.4 seed spreads | `python diag_ratio_spread.py` | seconds |
| Supplementary Table S4 (rebuild) | `python gen_bsh_table.py` | seconds |
| Supplementary Table S4 (catalogue from scratch) | `python domain_carriage.py --domains PF02275,PF02286,PF05985 --species 12 --per 12` | ~25 min on a 2-core Haswell CPU, plus ~150 MB of downloads |
| Table 5, Table 6 | `python cross_cohort9_seed.py` then `python heq_cd_only.py` | ~26 min and ~12 s on a 2-core Haswell CPU |
| Section 2.4 sensitivity grid | `python diag_m11_grid.py` | ~33 s |
| Section 2.4 window width | `python diag_window_width.py` | not timed here |
| Section 2.4 residence subsets | `python supplement_hstar.py` | not timed here |
| Table 4, section 2.7 | `python review_final.py` | not timed here (4,000 bootstrap draws) |
| Section 2.10, Supplementary Table S5 | `python m6_discriminant.py` | ~4 min on a 2-core Haswell CPU |
| Figures 1–6 | `FIG_OUT=./figures python make_figures.py` | seconds |

Run times are from the machine the analysis was done on (Intel i3-4170, Windows 11).
They scale with the single core the simulation uses.

---

## 4. Determinism

- Every tabulated threshold is the **median of eight seeds**, and the eight seeds are
  fixed constants, listed at the top of `cross_cohort9_seed.py`. Re-running on the same
  inputs reproduces the same tables by construction.
- The seed spread is reported next to every value rather than hidden: Supplementary
  Table S1 gives the min and max of the eight seeds in all 40 cells.
- Interpolation between the four simulated family sizes is log–log and is a
  convenience, not a fitted form; where the family size is known, read the table.

---

## 5. What is not here

- **The primary abundance matrices.** Eleven published tables, each owned by the study
  that produced it. `data_sources.tsv` gives the supplement or archive for each.

  `ELDERMET_species_189x836.csv` is the one that cannot be skipped. Fourteen scripts
  name it, and eight of them open it without a missing-file guard, so a clean download
  stops with `FileNotFoundError` rather than reporting `[缺文件]`; it is also the matrix
  behind Table 4 and section 2.7. It is Supplementary
  file 2 of Ghosh et al. (2020, *eLife* 9:e50240, CC BY 4.0), available from
  <https://cdn.elifesciences.org/articles/50240/elife-50240-supp2-v1.xlsx>
  (635,165 bytes, md5 `c3c7f79f2c8ead99407f3f905abe9836`). Take **worksheet `B`** and
  write it out as comma-separated values with a leading `SampleID` column: the sheet
  carries the 189 sample identifiers in column A but leaves cell A1 empty, so the
  header row needs that name supplied. The file that results is 584,146 bytes, md5
  `16c2683a6c18c94867e1d8b34bf8592a`, and reproduces every number in Table 4.
  `review_final.py` prints these same instructions when the file is absent.

- **The 122-sample Guangzhou metadata table** published by curatedMetagenomicData
  Curation. The 49 case IDs derived from it are included; the full table should be
  taken from the repository.
- **Raw reads.** The five ELDERMET runs used to calibrate detection efficiency are ENA
  runs ERR3986546, ERR3986547, ERR3986549, ERR3986555 and ERR3986557 (PRJEB37017).

---

## 6. Licence

Code (`*.py`) is released under the MIT licence, see `LICENSE`.
The tables and result files (`*.tsv`, `*.csv`, `*.json`, `*.md`) are released under
CC BY 4.0.

`TableS1_species_fixed.csv` is derived from Additional file 1 of Mou et al. (2021,
*BMC Genomics* 22:695, CC BY 4.0) with GTDB r95 taxonomy appended; please cite that
paper alongside this repository if you reuse it.

## 7. Citation

Citation details will follow the article. Until it is published, please cite this
repository by its DOI.
