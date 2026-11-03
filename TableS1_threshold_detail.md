**Supplementary Table S1. Conditional thresholds in ten cohort-settings, with the spread over eight seeds.**

Every entry is `median (min–max; range as % of the median)` over the eight simulation seeds, at the family size given by the column. `H*` is the conditional threshold of §2.4: the largest `H` at which P(ρ ≥ 0.90 | H ≈ h) ≥ 0.95, with a moving window of width 0.02, `min_n` = 40 and 1,200 simulated families per level. Cohorts are the ten settings of Table 5: nine cohorts plus one age-restricted subset of SG90.

| Cohort | *n* | Profiler | Species | Pool | *H*\*(3) | *H*\*(5) | *H*\*(10) | *H*\*(20) |
|---|---|---|---|---|---|---|---|---|
| ELDERMET (Ireland) | 189 | MetaPhlAn 3 | 836 | 149 | 0.0240 (0.0213–0.0453; 100%) | 0.0378 (0.0340–0.0506; 44%) | 0.0732 (0.0545–0.0810; 36%) | 0.0936 (0.0848–0.0978; 14%) |
| CRC collection, 9 cohorts | 1,262 | MetaPhlAn 3 | 917 | 164 | 0.0300 (0.0234–0.0579; 115%) | 0.0459 (0.0362–0.0589; 49%) | 0.0856 (0.0769–0.1025; 30%) | 0.1056 (0.0993–0.1163; 16%) |
| Guangzhou (China) | 49 | MetaPhlAn 3 | 575 | 121 | 0.0584 (0.0241–0.0831; 101%) | 0.0413 (0.0290–0.0509; 53%) | 0.0669 (0.0554–0.0897; 51%) | 0.0883 (0.0828–0.0977; 17%) |
| SG90 (Singapore), all | 516 | MetaPhlAn 4 | 382 | 69 | 0.0447 (0.0323–0.0872; 123%) | 0.0798 (0.0636–0.0965; 41%) | 0.1064 (0.1020–0.1169; 14%) | 0.0994 (0.0964–0.1024; 6%) |
| SG90, aged ≥77 | 221 | MetaPhlAn 4 | 382 | 63 | 0.0369 (0.0264–0.0499; 64%) | 0.0645 (0.0422–0.0967; 85%) | 0.1067 (0.1003–0.1126; 12%) | 0.0907 (0.0869–0.0920; 6%) |
| Villages (Honduras) | 1,787 | MetaPhlAn 4 | 1,910 | 467 | 0.0285 (0.0209–0.0491; 99%) | 0.0443 (0.0310–0.0540; 52%) | 0.0766 (0.0605–0.0860; 33%) | 0.0921 (0.0860–0.1003; 16%) |
| Screening cohort (Japan) | 615 | MetaPhlAn 4 | 1,835 | 300 | 0.0318 (0.0227–0.0399; 54%) | 0.0413 (0.0314–0.0837; 127%) | 0.0714 (0.0549–0.0922; 52%) | 0.1081 (0.0894–0.1185; 27%) |
| Wave 3 (China), aged ≥65 | 59 | mOTU | 3,694 | 128 | 0.0270 (0.0225–0.0321; 35%) | 0.0390 (0.0310–0.0619; 79%) | 0.0834 (0.0660–0.0938; 33%) | 0.0940 (0.0886–0.1015; 14%) |
| Xi'an cohort (China) | 171 | mOTU | 359 | 211 | 0.0338 (0.0200–0.0549; 103%) | 0.0525 (0.0354–0.0748; 75%) | 0.1045 (0.0689–0.1171; 46%) | 0.1125 (0.0934–0.1330; 35%) |
| Coronary cohort (China) | 102 | mOTU | 284 | 116 | 0.0241 (0.0218–0.0336; 49%) | 0.0367 (0.0288–0.0496; 57%) | 0.0806 (0.0534–0.0932; 49%) | 0.0936 (0.0879–0.0962; 9%) |

Pooled ranges across the ten settings, by family size: `n_sp` = 3, 35–123%; `n_sp` = 5, 41–127%; `n_sp` = 10, 12–52%; `n_sp` = 20, 6–35%. Over all 40 cells the range is 6–127% with a median of 47.5%. The `n_sp` = 3 column is weakly identified (it moves 1.8-fold with the width of the moving window, against 1.1-fold at `n_sp` = 20; §2.4) and is reported as an order-of-magnitude anchor only. No comparison in the paper rests on it; the single use is the robustness bound in §2.4, where the *hdcA* family still fails at the largest of the eight seed values at that size.

The Guangzhou row uses the 49 samples that the curated metadata labels as cases (Methods §2.1.6d); its *H*\*(3) exceeds its *H*\*(5), which is why the span of that cohort is quoted as max/min across the four sizes rather than as a ratio of endpoints (Table 5, note).

---

**Columns carried by the deposited TSV.** Beyond the four `H*` columns tabulated above, `TableS1_threshold_detail.tsv` reports, for each cohort-setting, the observed `c_eff`, `H`, `rho_all` and `rho_carriers` of the *hdcA* family, the ratio `H`/*H*\*, and three counts: `hdca_matched` (catalogue carriers present in the cohort's species table), `contributing_species` (those of them with a non-zero cohort-mean abundance), and `ratio_nsp` (the family size at which the ratio was evaluated).

The two new columns answer a reading that the older file invited. `ratio` is computed once per cohort, at `ratio_nsp`, and is therefore identical on all four `n_sp` rows of that cohort — not because the threshold is flat, but because the ratio does not vary with the row. The row that matters for the ratio is `n_sp = ratio_nsp`; the four `H*` columns vary with `n_sp` as they should. `contributing_species` makes the distinction visible: ELDERMET resolves 20 catalogue carriers of which 5 carry non-zero abundance, so its reported ratio of 1.63 is evaluated at `n_sp` = 20 while its §2.6 discussion quotes the same failure as 4.0-fold at `n_sp` = 5. Recomputing every cohort at its contributing count changes no verdict (six of nine fail under either convention; `recalc_contributing.py` in the archive).

In `SG90, aged ≥77` the two counts differ by one (14 matched, 13 contributing): one carrier has zero mean abundance within the age-restricted subset, so the ratio for that subset is evaluated at 14 and is, if anything, slightly conservative.
