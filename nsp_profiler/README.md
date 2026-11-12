# Family size under two profiler versions, and the criterion on 538 families

Two pieces of work sit in this folder. The first is the evidence behind the paragraph in
§3.2 of the manuscript that asks whether a family size taken from a profiler holds still
across pipeline versions: it does not, and the numbers are in `humann3/`. The second is a
wider panel that applies the same criterion to 538 families with known carriers, and then
re-runs it under a second cohort's abundance weights and at the level of individuals; it is
in `multi_family/`. The panel is not quoted in the manuscript, and is deposited so that the
claim that the criterion is not specific to *hdcA* and bile salt hydrolase can be checked
rather than taken on trust.

Neither folder redistributes primary data. Both read the HMP2 / IBDMDB release, which is
public but large, and give the paths below for refetching it.

---

## 1. `humann3/` — HUMAnN 2 and HUMAnN 3 on the same five samples

### Data

The HMP2 (IBDMDB) release is served from a public Globus endpoint,
`g-227ca.190ebd.75bc.data.globus.org`, no authentication required. Three things are used:

| Path under `products/HMP2/MGX/` | What it is |
|---|---|
| `2018-05-04/taxonomic_profiles_3.tsv.gz` | MetaPhlAn 3 species table, 1,627 samples, relative abundance in per cent |
| `func_profile3_WGS/<SID>_humann3.tar.bz2` | HUMAnN 3 output, ~1.3 MB per sample |
| `2018-05-04/func_profiles/<SID>_humann2.tar.bz2` | HUMAnN 2 output, same sample set, older release |

The two HUMAnN directories sit at different depths — `func_profiles` is nested under
`2018-05-04/` while `func_profile3_WGS` is at the top level. Requesting the v2 packages at
the top level returns 404, which is worth knowing before concluding the endpoint is down.

Each package carries `*_genefamilies_cpm.tsv` (HUMAnN 3) or `*_genefamilies.tsv`
(HUMAnN 2), and both are **stratified**: a line for a family as a whole, then one line per
contributing species, written `UniRef90_XXXX|g__Genus.s__Species`. That stratified layer is
the `c_s · p_s` decomposition the criterion audits, so these files can be scored directly.
HUMAnN 2 reports RPK and HUMAnN 3 CPM; every quantity used here is normalised within the
family before it is compared, so the unit does not carry across.

### Scripts

| File | Role |
|---|---|
| `fetch_hmp2.py` | Downloads v3 packages. Sequential, one sample at a time. |
| `fetch_hmp2_v2.py` | Downloads the v2 packages from the deeper path. |
| `fetch_parallel.py` | Replaces the two above for large pulls: HEAD for length, K byte ranges in parallel against **one** file, per-chunk retry, merge, then `tar -tjf` before the package is given its real name. The endpoint throttles per connection (3–9 KB/s single); 8–16 ranges against one file reach 30–100 KB/s, while opening several connections to *different* files gives no speed-up at all. |
| `pair_v2v3.py` | Pairs the two versions. Sample names carry `_P` / `_TR` in v3 and not in v2, so pairing is on the base name. Carriers are chosen by a relative threshold — contribution ≥ 0.1% of the family total — rather than a CPM cut, because the two versions use different units and an absolute cut would select different species sets. Computes *H*, *H*\*, λ and the verdict per family-sample pair. Writes `out/paired_v2v3.tsv`. |
| `analyze_humann3.py` | The v3-only summary: family counts, the distribution of *H*/*H*\*, and the proxy error by *H*/*H*\* bin. Writes `out/humann3_families.tsv`, `out/humann3_by_sample.tsv`, `out/humann3_summary.txt`. |
| `decompose_v2v3.py` | The version comparison in the manuscript. Splits each change of verdict into the part caused by the family size moving and the part caused by *H* moving. Writes `out/decompose_v2v3.tsv`. |

### How the decomposition works

A verdict is *H*/*H*\*(*n*<sub>sp</sub>) > 1. Between v2 and v3 both terms move: the newer
pangenome resolves more species per family, and the intensities inside the carrier set are
estimated differently. Holding *H*\* at the v2 family size separates them:

```
ratio2   = H2 / H*(n2)      the v2 verdict
ratio3   = H3 / H*(n3)      the v3 verdict
ratio3n2 = H3 / H*(n2)      counterfactual: v3 intensities, threshold at the v2 family size
```

The flip between `ratio3` and `ratio2` is the total; `ratio3n2` against `ratio2` isolates the
change in *H*; `ratio3` against `ratio3n2` isolates the change in family size. The two do not
add up, because a family can move both ways at once, so the script reports all four
trajectory classes rather than two.

### Reproducing the paragraph in §3.2

```bash
python fetch_parallel.py --kind v3 --samples <SID> ... --blocks 16
python fetch_hmp2_v2.py                       # v2 packages
# unpack both into work/, keeping *_genefamilies_cpm.tsv (v3) and *_genefamilies.tsv (v2)
python pair_v2v3.py                           # -> out/paired_v2v3.tsv
python decompose_v2v3.py                      # defaults: --rel 0.001 --min-p 0.1 --min-sp 3
```

`humann3/decompose_default.log` is the console output of that last command at its default
arguments, which is the run the manuscript quotes. Every figure in the paragraph can be read
off it:

| Manuscript | Where it comes from |
|---|---|
| "2,187 family-sample pairs" | `家族-样本对 n=2187` |
| "agreed … for 42%" | `n_sp 在 v2/v3 间相同：923 / 2187 = 42.2%` |
| "the median rising from three to five" | `n_sp 中位：v2 = 3，v3 = 5` |
| "30% of verdicts changed, 464 towards admission against 198 towards rejection" | `整体翻转方向：fail->pass 464，pass->fail 198`, total 662 / 2,187 = 30.3% |
| "49% were admitted … against 4%" | The boundary-distance block, `>2` row: `>=+1` 49.0% against `0` 4.2% |
| "*H* was 6% higher where the carrier set was unchanged" | The `Δn_sp = 0` row: `ΔH/H2中位% = 5.9%` |

The 24.4% / 5.9% split between the *H* channel and the family-size channel is in the
`判决轨迹分类` line. Note that the largest single class is "neither" (67.9%), and that the
two channels are reported because they are nearly additive here, not because one dominates.

### What this folder cannot support

- **Five pairs.** The two directories overlap in 1,317 samples, but at 3–9 KB/s per
  connection the pull is a two-hour job and only five pairs were recovered. Per-sample family
  overlap ranges from 54 to 797; the 54-family sample (`HSMA33KQ`) has ρ(*H*<sub>2</sub>,
  *H*<sub>3</sub>) = −0.04 and supports no reading at all. The direction is reliable — the
  family size rises in 5 of 5 samples — the magnitude is not.
- **The six-sample summary moves.** `out/humann3_summary.txt` is written on whatever is
  unpacked in `work/`. A seventh sample could not be recovered intact; with six instead of
  seven, the share of families with *H*/*H*\* > 1 moves from 58.8% to 67.2% while the
  monotone rise of proxy error across bins holds in both. Quote a range, not a point, and
  say how many samples.
- **Copy number and carriage are not separated.** A per-species contribution is
  proportional to copies in that genome times that genome's abundance. The continuous
  intensity measured here is the equivalent of the strain-level carriage in the main text,
  not the same quantity, and nothing in this folder separates the two factors.

---

## 2. `multi_family/` — the criterion on 538 families with known carriers

Not quoted in the manuscript. Included because the manuscript's worked examples are one
gene (*hdcA*) and one domain family (bile salt hydrolase), and a reader is entitled to ask
whether the failure mode generalises.

### How the panel is built

`fetch_pfam_catalog.py` walks the InterPro API for the Pfam catalogue, 53,275 entries; the
deposited `pfam_catalog.tsv` is the result, so the walk does not have to be repeated.
`fetch_hmms.py` pulls the HMM for each candidate family from InterPro and checks the
decompressed file for an `HMMER3` header. A failed download is still a valid gzip of an
error page, so file size and gzip integrity together prove nothing. `run_panel_mp.py`
scans the twelve most abundant ELDERMET species' proteomes against each HMM with pyhmmer,
eight processes, and resolves carriage per species per family. `summarize_panel.py` applies
the two constraints and sorts the failures by cause.

Requires `pyhmmer` 0.12.3; everything else is standard library. The HMM files themselves
(25 MB) are not deposited, because InterPro serves them and the catalogue file records
which ones are needed.

### Scripts

| File | Role |
|---|---|
| `fetch_pfam_catalog.py` | InterPro API walk, cursor-paginated. `?page=N` is accepted and ignored by the endpoint, so a page-numbered loop fetches the first page once per iteration and returns 53,275 rows that are the same 440 rows repeated — use the `next` cursor. |
| `fetch_hmms.py` | Per-family HMM download and validation. |
| `run_multi_family.py` | Single-process scan, used for development. |
| `run_panel_mp.py` | The same scan with a process pool; this produced `out/carriage_panel.tsv`. |
| `summarize_panel.py` | Two-constraint verdicts, failure classes, and the quadrant counts. Writes `out/families_panel.tsv`, `out/panel_fail.tsv`, `out/panel_summary.txt`. |
| `run_crosscohort.py` | Re-runs the verdicts with HMP2 abundance weights in place of ELDERMET's, at three species-presence thresholds. Writes `out/crosscohort_hmp2.tsv` and the two sensitivity files. |
| `run_individual.py` | Per-individual verdicts for 1,581 HMP2 individuals across the panel, and the agreement between the cohort-mean verdict and the individual majority. Writes `out/individual_hmp2.tsv`, `out/individual_summary.txt`. |
| `check_r.R` | Cross-checks the Spearman correlations in `panel_summary.txt` against R's `cor.test`. |

### Headline numbers

- 1,410 families scanned, **538 with at least one carrier**. Under both constraints: 404
  pass, 112 fail on fidelity alone, 22 on detectability alone, **none on both**.
- Above the range where the verdict can be decided at all (*n*<sub>sp</sub> ≥ 5), the
  detectability constraint never fires — the smallest λ<sub>det</sub> is 48.45 against a
  requirement of 5 — so the quadrant plot collapses to two cells. 71 of the 404 passes are
  *n*<sub>sp</sub> = 1 degenerate cases; the honest pass rate is 282/349 = 80.8%.
- Failures concentrate in carbohydrate-active families: GH18, GH38, GH67, GH123, SusG
  (CBM58), hyaluronidase, and the choline/ethanolamine kinases with the betaine ABC
  transporter — families already known to be species-specific in the literature.
- Swapping the abundance weights for HMP2's flips 55/538 = 10.2% of verdicts, and 88 of the
  112 fidelity failures recur (79%). The same at three presence thresholds. **But**
  *n*<sub>sp</sub> is identical for every family across the two cohorts — these twelve are
  core gut species — so this tests the weight structure only, not a change in species
  composition.
- Cohort-mean verdicts agree with the individual majority for 456/538 = 84.8%. Part of the
  disagreement is detectability moving with each individual's depth, which is why the number
  is reported rather than used.

---

## 3. Files

```
humann3/     fetch*.py, pair_v2v3.py, analyze_humann3.py, decompose_v2v3.py
             out/paired_v2v3.tsv, out/decompose_v2v3.tsv,
             out/humann3_families.tsv, out/humann3_by_sample.tsv,
             out/humann3_summary.txt, decompose_default.log
multi_family/ fetch_pfam_catalog.py, fetch_hmms.py, run_*.py, summarize_panel.py, check_r.R
             pfam_catalog.tsv
             out/carriage_panel.tsv, out/families_panel.tsv, out/panel_fail.tsv,
             out/panel_summary.txt, out/crosscohort_hmp2.tsv,
             out/individual_hmp2.tsv, out/individual_summary.txt
```

`out/decompose_v2v3.tsv` is one row per family-sample pair with *n*<sub>sp</sub>, *H* and the
verdict under each of the three ratios, so any subset quoted above can be re-tallied by
filtering on the columns rather than by re-running the decomposition.
