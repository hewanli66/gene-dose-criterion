# 复现包（中文速览）

对应稿件：**Species abundance is not a gene dose: a two-constraint criterion for
functional inference in metagenomics**（中位数版 v3-M）。

英文正文见 `README.md`，本文件只做中文要点说明，便于读者核对与复用。

---

## 一、这个包能干什么

| 情况 | 能否直接跑 |
|---|---|
| 重建补充表 S1、S2，并与 Table 5 逐格核对 | **能**，几秒钟，包里数据齐全 |
| 重算正文 §2.4 引用的种子极差（35–123% 等） | **能**，几秒钟 |
| 重算 Table 5、Table 6 | 需要自备十一个丰度矩阵（见 `data_sources.tsv`），约 26 分钟（i3-4170 单核实测 1,557 秒） |
| 重出图 1–6 | 需要 ELDERMET 矩阵 |
| 从原始测序数据重跑 | 不能，本包不含 reads |

两条验证命令：

```bash
python gen_supp_tables.py
python diag_ratio_spread.py
```

第一条的第 1 行会打印 `Table S1 与 Table 5 逐格核对：40 格，0 处不一致`。
这句话就是这次补表的意义——补充表与主表来自不同的落盘文件，必须逐格对上。

## 二、为什么不把原始矩阵放进来

十一个丰度矩阵分属十项已发表研究，再分发涉及各自的许可与体量；正文
*Data availability* 也写明了"不重新分发"。所以包里只放：

- **代码**（27 个 .py，全部只用标准库）
- **派生产物**（阈值表、种子明细、敏感性网格等）
- **一张整理的携带率表** `TableS1_species_fixed.csv`（96 物种 × 8 列），
  它是分析的核心输入，性质上属于"派生表"而非原始丰度矩阵，故随包附带

`data_sources.tsv` 逐文件写明：这个文件属于哪个队列、出自哪份补充材料、
是否随包分发、被哪个脚本读取。

## 三、文件对应关系（速查）

| 稿件里的位置 | 包里的文件 |
|---|---|
| Table 5 | `cross_cohort9_seed_t5_heqcd.tsv` |
| Table 6 | `cross_cohort9_seed_t6.tsv` |
| Supplementary Table S1 | `TableS1_threshold_detail.tsv`（长表）+ `.md`（可读版） |
| Supplementary Table S2 | `TableS2_M11_grid.tsv` + `.md` |
| Table 1 的底表 | `TableS1_species_fixed.csv` |
| §2.4 的种子极差 | `cross_cohort9_seed_t5_detail.txt` + `heq_cd_only.txt` |

注意：Table 5 的广州行取自 `heq_cd_only.txt`（49 例），而
`cross_cohort9_seed_t5_detail.txt` 里广州仍是 122 份的旧口径。合并后的
`cross_cohort9_seed_t5_heqcd.tsv` 才是最终值，`gen_supp_tables.py` 已按此替换。

## 四、复现顺序

1. 按 `data_sources.tsv` 取回矩阵，文件名保持原样，放在本目录；
2. `python cross_cohort9_seed.py` → Table 5 与 Table 6 的底表；
3. `python heq_cd_only.py` → 用 49 例替换广州行，生成 `..._heqcd.tsv`；
4. `python gen_supp_tables.py` → S1、S2，并做 40 格核对；
5. 需要图时：`FIG_OUT=./figures python make_figures.py`。

缺矩阵不会崩，脚本会打印 `[缺文件]` 并跳过该项，可分批跑。

## 五、没放进去的东西

- 十一个丰度矩阵（见上）
- curatedMetagenomicData 发布的广州 122 份整理元数据；包里只放由它导出的
  49 个 case ID（`HeQ2017_CD_ids.txt`），`heq_cd_only.py` 在元数据缺失时自动回退到这个文件
- 原始 reads（五个 ELDERMET run 见 README §5）

## 六、许可

代码 MIT；表与结果文件 CC BY 4.0。`TableS1_species_fixed.csv` 派生自
Mou et al. 2021（BMC Genomics 22:695，CC BY 4.0），复用时请一并引用原文。
