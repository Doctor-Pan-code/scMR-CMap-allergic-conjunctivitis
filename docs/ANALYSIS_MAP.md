# Analysis map

Method section of the manuscript → script that implements it → table it writes → figure or table
in the paper that consumes it. Nothing in the manuscript is computed outside this table.

Path prefixes: `scripts/`, `tables/`, `figures/`.

---

## Day 1 — data sources and the druggable genome pool

| Manuscript section | Code | Output | Consumed by |
| --- | --- | --- | --- |
| 方法 · 数据来源与可成药基因组池 | `scripts/day1_druggable_pool.py` | `tables/day1_hgnc_protein_coding.csv`, `tables/day1_druggable_pool.csv`, `tables/day1_druggable_sources.csv`, `tables/day1_pool_summary.csv` | 表 1, `figures/fig1_study_design` |

`load_hgnc` keeps approved protein-coding records that carry an Ensembl gene ID (19,255).
`load_dgidb_sources` and `load_dgidb_interactions` classify each gene into tiers T1/T2/T3.
`_chrom` and `_median` are helpers used by the summary. The tier assignment is a pure join:
no statistical model is fitted here, which is why Day 1 has no p-values.

---

## Day 2 — cis-eQTL Mendelian randomization

| Manuscript section | Code | Output | Consumed by |
| --- | --- | --- | --- |
| 方法 · cis-eQTL 孟德尔随机化 | `scripts/day2_mr_screen.py` | `tables/day2_screen_lead_snp_mr.csv`, `tables/day2_instruments_lead.csv`, `tables/day2_screen_summary.csv`, `tables/day2_mr_confirmatory.csv`, `tables/day2_mr_heterogeneity.csv`, `tables/day2_mr_leave_one_out.csv`, `tables/day2_mr_steiger.csv` | 表 1, `figures/fig2_mr_screen` |

Two stages, both inside this one script.

**Stage A (`stage_a`)** — one lead cis-eQTL instrument per gene, Wald ratio against the
allergic-conjunctivitis endpoint. Instruments come from `collect_lead_snps`, are restricted to the
cis window by `cis_filter`, deduplicated by `positional_prune`, and allele-harmonised by `align`.
6,329 genes are screenable; BH FDR is computed across all of them.

**Stage B (`stage_b`)** — every gene at screening P < 0.05 is re-estimated with all independent cis
instruments under six estimators. `_outcome_meta` reads the endpoint metadata once so the stage
does not re-query per gene. Cochran's Q, the MR-Egger intercept, leave-one-out, Steiger and the
replication endpoint are all emitted as separate tables so each diagnostic can be inspected alone.

Estimators themselves live in `scripts/lib_stats.py`: `ivw`, `mr_egger`, `weighted_median`,
`simple_mode`, `weighted_mode`, `leave_one_out`, `f_statistic`, `steiger`.

`figures` (lower-case, inside the script) writes the intermediate screen figure data.

---

## Day 3 — colocalization and SMR/HEIDI

| Manuscript section | Code | Output | Consumed by |
| --- | --- | --- | --- |
| 方法 · 共定位与 SMR/HEIDI | `scripts/day3_coloc_smr.py` | `tables/day3_coloc_results.csv`, `tables/day3_coloc_variant_level.csv`, `tables/day3_eqtl_source.csv`, `tables/day3_smr_heidi.csv`, `tables/day3_candidate_evidence.csv` | 表 2, `figures/fig3_coloc_smr` |

`candidates` reads the Day 2 confirmatory table and keeps the genes that clear the screening
threshold. `gtex_gencode` maps Ensembl IDs to GTEx v8 identifiers; `gtex_cis_variants` pulls the
dense cis grid; `pick_tissue` chooses the GTEx tissue with the most usable variants for that gene.
`variant_position` and `variant_alleles` resolve the grid against Ensembl, and `to_effect` converts
GTEx's normalised effect scale to the log-odds scale that `coloc.abf` expects. `thin` and `align`
reduce the grid to the shared variant set with concordant alleles — this is what makes the
colocalization assumption checkable rather than assumed.

Posteriors come from `scripts/lib_stats.py:coloc_abf`. SMR and HEIDI come from
`lib_stats.py:smr_wald` and `lib_stats.py:heidi`. `_f`, `_i`, `_z`, `_z_from_p` and
`matplotlib_cmap` are parsing and colour helpers.

`day3_candidate_evidence.csv` is the merge point: MR, coloc and SMR/HEIDI side by side, with
`n_supporting_axes` and the `prioritised` flag that decides who enters Day 4.

---

## Day 4 — single-cell expression-state signature

| Manuscript section | Code | Output | Consumed by |
| --- | --- | --- | --- |
| 方法 · 单细胞谱系内分层与跨谱系重现 | `scripts/day4_scrna_signature.py` | `tables/day4_cell_qc.csv`, `tables/day4_lineage_scores.csv`, `tables/day4_lineage_feasibility.csv`, `tables/day4_stratification_design.csv`, `tables/day4_group_balance.csv`, `tables/day4_pseudobulk_groups.csv`, `tables/day4_signature_up.csv`, `tables/day4_signature_dn.csv`, `tables/day4_signature_summary.csv`, `tables/day4_lineage_filter_audit.csv`, `tables/day4_cross_lineage_signature.csv`, `tables/day4_cross_lineage_replication.csv`, `tables/day4_sensitivity_signature.csv`, `tables/day4_null_calibration.csv`, `tables/day4_target_gene_by_cluster.csv`, `tables/day4_cluster_composition.csv` | 表 3, `figures/fig4_singlecell` |
| 方法 · 单细胞谱系内分层与跨谱系重现 (candidate probe) | `scripts/day4_probe_genes.py` | `tables/day4_probe_genes.csv`, `tables/day4_probe_cells.csv`, `tables/day4_probe_gene_scores.csv` | `figures/fig4_singlecell` |

`day4_scrna_signature.py` is the largest script in the repository and deliberately avoids loading
the whole matrix: `mtx_dims` reads the MatrixMarket header, `iter_mtx` streams the file in chunks,
and `_parse_chunk` decodes each chunk in place. `load_features`, `load_barcodes`, `load_clusters`
and `load_locus_group` build the cell metadata; `candidate_table` assembles the 20 prioritised
genes; `pick_target` selects the stratifier; `residualise` regresses out library size and lineage
identity so the split is a state contrast, not a cell-type contrast.

Three controls sit in the same script and are the reason the signature is credible:

- `day4_null_calibration.csv` — the stratifier is permuted and the whole pipeline re-run; the
  resulting correlation (+0.010) shows the split procedure itself carries no direction.
- `day4_cross_lineage_replication.csv` — the same 600 genes re-tested inside six lineages; 73.5 %
  directional concordance.
- `day4_sensitivity_signature.csv` — half-split stability across two independent cell halves.

`day4_lineage_feasibility.csv` is the table that forces the design choice: detection of the target
runs 11.6–20.6 % in every lineage, so no lineage can serve as the contrast.

`day4_probe_genes.py` re-streams the matrix once more to emit per-cell CPM for all 20 candidates.
It is separate from the signature script because the signature run does not need per-cell output
and holding 27,600 × 20 in memory alongside the streaming state was the difference between a
4 GB and a 900 MB peak.

---

## Day 5 — connectivity-map reversal

| Manuscript section | Code | Output | Consumed by |
| --- | --- | --- | --- |
| 方法 · 连接图谱逆转筛选 | `scripts/day5_cmap_reversal.py` | `tables/day5_query_size_calibration.csv`, `tables/day5_null_by_size.csv`, `tables/day5_cmap_connectivity.csv`, `tables/day5_cmap_summary.csv`, `tables/day5_l1000cds2_raw.csv`, `tables/day5_positive_control_check.csv`, `tables/day5_enrichr_library_status.csv`, `tables/day5_enrichr_lincs_up.csv`, `tables/day5_enrichr_lincs_dn.csv` | 表 4, `figures/fig5_cmap` |

`load_landmarks` and `load_signature` intersect the 600-gene signature with the 978 L1000 landmark
genes (44 survive). `query` posts both orientations to L1000CDS2; `top_score` reduces each response
to a single connectivity score; `record` writes the row. `gluco_in` checks whether a glucocorticoid
appears in the top of the list — the positive control that has to pass before any other result is
read. `cls_of` assigns the perturbagen class used in the figure.

`day5_query_size_calibration.csv` is the known-answer experiment: the dexamethasone signature is
pushed through the identical pipeline at four query sizes, and the size at which glucocorticoids
reliably come back (15 genes per arm) is the size used for the real query. The permutation null
(`day5_null_by_size.csv`, 20 draws per size) is generated at the same sizes, so the observed score
and the null are compared at matched resolution.

Connectivity and reversal specificity are computed by `lib_stats.py:connectivity_score`;
`lib_stats.py:tau_score` provides the rank-agreement cross-check.

---

## Day 6 — phenome-wide safety sweep

| Manuscript section | Code | Output | Consumed by |
| --- | --- | --- | --- |
| 方法 · 全表型安全性扫描与五轴合成 | `scripts/day6_phewas_safety.py` | `tables/day6_phenotype_universe.csv`, `tables/day6_phenotypes_scanned.csv`, `tables/day6_phewas_all.csv`, `tables/day6_phewas_hits.csv`, `tables/day6_gene_safety_verdict.csv` | 表 5, `figures/fig6_safety_integration` |

`phenotype_universe` enumerates FinnGen R13 endpoints; `select_phenotypes` applies the case-count
filter, leaving 1,400 scanned endpoints. `fetch_finngen_meta` pulls the endpoint metadata once.
`domain_of` and `severity_of` annotate each endpoint from keyword rules in
`config/config.json` — the annotation is data-driven, not hand-assigned per hit.
`candidate_genes` reads the 14 genes that reached Day 6.

The verdict is a veto, not a score: one severe systemic domain association removes the candidate
regardless of genetic strength. 10 genes are safe, 3 monitor, 1 high-risk.

---

## Day 7 — five-axis integration

| Manuscript section | Code | Output | Consumed by |
| --- | --- | --- | --- |
| 方法 · 全表型安全性扫描与五轴合成 | `scripts/day7_integration.py` | `tables/day7_axis_scores.csv`, `tables/day7_composite_ranking.csv`, `tables/day7_candidate_drug_pairs.csv`, `tables/day7_project_manifest.csv` | 表 2, `figures/fig6_safety_integration` |

`load_axes` reads the four upstream tables. `scale_mr`, `scale_coloc` and `scale_safety` map each
axis to [0, 1]; the single-cell and CMap axes are already on that scale. Weights come from
`config/config.json:scoring.weights` (MR 0.30, coloc 0.20, single cell 0.20, CMap 0.20, safety
0.10) — they are read, never typed into the script. `_f` is a float parser.

`day7_composite_ranking.csv` is the ranked table; `candidate_tier` uses the Tier 1 cut-off of 0.55.
TLR1 reaches 0.7631 and is the only Tier 1 entry.

---

## Cross-cutting infrastructure

| File | Role |
| --- | --- |
| `scripts/lib_io.py` | `ROOT` resolution (relative to the file, so the repository runs anywhere), config loading, CSV/JSON read-write, logging, SVG/PNG saving, `token()` |
| `scripts/lib_api.py` | cached HTTP clients for OpenGWAS, Ensembl, L1000CDS2, Enrichr and DGIdb; every response is cached under `data/intermediate/cache/` |
| `scripts/lib_stats.py` | the whole statistics engine in pure Python: `ivw`, `mr_egger`, `weighted_median`, `simple_mode`, `weighted_mode`, `leave_one_out`, `f_statistic`, `steiger`, `coloc_abf`, `heidi`, `smr_wald`, `connectivity_score`, `tau_score`, `bh_fdr`, `chi2_sf`, `t_sf`, `norm_ppf` |
| `scripts/run_pipeline.py` | sequential runner for Days 3–7, one log per step, continues past a failing step |

---

## Report build

| Manuscript section | Code | Output |
| --- | --- | --- |
| 图注 / 图件 | `scripts/report/make_figs.py` | `figures/fig1..fig6 .png` + `.pdf`, `figures/fig_audit.json` |
| whole manuscript | `scripts/report/build_manuscript.py` | `manuscript_out/_build.docx` |

`make_figs.py` exposes `fig1` … `fig6` and a `main()` that writes all six. It reads `tables/*.csv`
at run time through `load(name)` / `t(name)`, asserts a 165 mm canvas and 7.0–8.0 pt type, and
prints a per-figure audit line reporting canvas width, the type-size range, the number of elements
outside the canvas and the number of text collisions. Every run in this repository reports zero for
the last two.

`build_manuscript.py` assembles the submission `.docx` from `manuscript_src/_norm.txt` (run-level
normalised text with style markers), `manuscript_src/_refs_nature.json` (formatted reference list)
and `manuscript_src/tpl.docx` (the Scientific Reports Word template). It writes to
`manuscript_out/_build.docx`. Figures are read from `figures/`.
