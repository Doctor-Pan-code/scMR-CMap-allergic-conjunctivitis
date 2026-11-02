# Data dictionary

Every result table the manuscript draws on. Row counts are from the archived run; column names are verbatim. All CSVs are UTF-8 with BOM so Excel opens them correctly, and all floats are written at six significant figures.

Path prefix: `tables/`.

---

## `day1_druggable_pool.csv`

**Rows:** 10,983

One row per druggable gene. `druggable_tier` is T1 (clinical-stage annotation plus a recorded drug interaction), T2 (clinical annotation only) or T3 (literature or experimental evidence only); `has_cis_eqtl` marks whether an OpenGWAS `eqtl-a-ENSG*` dataset exists and therefore whether the gene can be instrumented at all.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `hgnc_id` |
| `name` |
| `location` |
| `gene_group` |
| `alias_symbol` |
| `prev_symbol` |
| `entrez_id` |
| `uniprot_ids` |
| `dgidb_sources` |
| `clinical_sources` |
| `n_sources` |
| `n_drug_interactions` |
| `n_drugs` |
| `example_drugs` |
| `druggable_tier` |
| `eqtl_dataset` |
| `has_cis_eqtl` |
| `chromosome` |

## `day1_druggable_sources.csv`

**Rows:** 39

Per-source-database gene counts, with the evidence class each database contributes.

| Column |
| --- |
| `source_db` |
| `n_genes` |
| `class` |

## `day1_hgnc_protein_coding.csv`

**Rows:** 19,255

HGNC approved, protein-coding records carrying an Ensembl gene ID. This is the 19,255-gene starting universe.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `hgnc_id` |
| `name` |
| `location` |
| `gene_group` |
| `alias_symbol` |
| `prev_symbol` |
| `entrez_id` |
| `uniprot_ids` |

## `day1_pool_summary.csv`

**Rows:** 4

Tier-level counts and medians behind Table 1 of the manuscript.

| Column |
| --- |
| `tier` |
| `n_genes` |
| `n_with_cis_eqtl` |
| `median_n_drugs` |
| `median_n_interactions` |

## `day2_instruments_lead.csv`

**Rows:** 7,006

The lead instrument for every screened gene, with its exposure and outcome effect pair, kept so the screen can be audited SNP by SNP.

| Column |
| --- |
| `id` |
| `trait` |
| `chr` |
| `position` |
| `rsid` |
| `ea` |
| `nea` |
| `eaf` |
| `beta` |
| `se` |
| `p` |
| `n` |
| `symbol` |
| `ensembl_gene_id` |
| `druggable_tier` |
| `n_cis_snp` |

## `day2_mr_confirmatory.csv`

**Rows:** 494

Stage B: for every gene reaching screening P < 0.05, all independent cis instruments under six estimators, plus Cochran Q, the MR-Egger intercept, leave-one-out support, replication in an independent conjunctivitis endpoint, the Steiger directionality test and BH FDR.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `druggable_tier` |
| `eqtl_dataset` |
| `nsnp` |
| `snps` |
| `f_statistic` |
| `screen_wald_p` |
| `screen_wald_or` |
| `IVW_FE_b` |
| `IVW_FE_se` |
| `IVW_FE_p` |
| `IVW_FE_or` |
| `IVW_FE_or_lci` |
| `IVW_FE_or_uci` |
| `IVW_RE_b` |
| `IVW_RE_se` |
| `IVW_RE_p` |
| `IVW_RE_or` |
| `IVW_RE_or_lci` |
| `IVW_RE_or_uci` |
| `MR_Egger_b` |
| `MR_Egger_se` |
| `MR_Egger_p` |
| `MR_Egger_or` |
| `MR_Egger_or_lci` |
| `MR_Egger_or_uci` |
| `WeightedMedian_b` |
| `WeightedMedian_se` |
| `WeightedMedian_p` |
| `WeightedMedian_or` |
| `WeightedMedian_or_lci` |
| `WeightedMedian_or_uci` |
| `SimpleMode_b` |
| `SimpleMode_se` |
| `SimpleMode_p` |
| `SimpleMode_or` |
| `SimpleMode_or_lci` |
| `SimpleMode_or_uci` |
| `WeightedMode_b` |
| `WeightedMode_se` |
| `WeightedMode_p` |
| `WeightedMode_or` |
| `WeightedMode_or_lci` |
| `WeightedMode_or_uci` |
| `egger_intercept` |
| `egger_intercept_p` |
| `Q` |
| `Q_pval` |
| `replication_nsnp` |
| `replication_b` |
| `replication_p` |
| `replication_or` |
| `steiger_r2_exposure` |
| `steiger_r2_outcome` |
| `steiger_correct_direction` |
| `IVW_RE_fdr` |

## `day2_mr_heterogeneity.csv`

**Rows:** 527

Per-SNP Wald ratios for the confirmatory genes - the input to Cochran Q.

| Column |
| --- |
| `symbol` |
| `rsid` |
| `bx` |
| `by` |
| `sx` |
| `sy` |
| `ratio` |

## `day2_mr_leave_one_out.csv`

**Rows:** 9

Leave-one-out IVW estimates for the genes where a single instrument could have driven the result.

| Column |
| --- |
| `symbol` |
| `dropped_rsid` |
| `ivw_b` |
| `ivw_se` |
| `ivw_p` |

## `day2_mr_steiger.csv`

**Rows:** 494

Steiger r-squared for exposure and outcome, with the direction verdict.

| Column |
| --- |
| `symbol` |
| `r2_exposure` |
| `r2_outcome` |
| `correct_direction` |

## `day2_screen_lead_snp_mr.csv`

**Rows:** 6,329

Stage A screening: one lead cis-eQTL instrument per gene, Wald ratio against the FinnGen allergic-conjunctivitis endpoint. `wald_fdr` is Benjamini-Hochberg across the screened genes; `f_statistic` is the first-stage strength of the lead instrument.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `druggable_tier` |
| `eqtl_dataset` |
| `rsid` |
| `chr` |
| `position` |
| `ea` |
| `nea` |
| `eaf_exposure` |
| `beta_exposure` |
| `se_exposure` |
| `p_exposure` |
| `n_exposure` |
| `beta_outcome` |
| `se_outcome` |
| `p_outcome` |
| `eaf_outcome` |
| `wald_b` |
| `wald_se` |
| `wald_or` |
| `wald_or_lci` |
| `wald_or_uci` |
| `wald_p` |
| `f_statistic` |
| `wald_fdr` |

## `day2_screen_summary.csv`

**Rows:** 6

Headline screen counts: genes screened, counts at each significance threshold, median F.

| Column |
| --- |
| `metric` |
| `value` |

## `day3_candidate_evidence.csv`

**Rows:** 50

The merged candidate table: MR, coloc and SMR/HEIDI evidence side by side, with `n_supporting_axes` and the `prioritised` flag that decides who enters Day 4.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `druggable_tier` |
| `screen_p` |
| `ivw_p` |
| `ivw_or` |
| `nsnp_mr` |
| `coloc_PP_H4` |
| `coloc_grade` |
| `smr_p` |
| `heidi_p` |
| `smr_heidi_pass` |
| `n_supporting_axes` |
| `prioritised` |

## `day3_coloc_results.csv`

**Rows:** 50

coloc.abf posterior probabilities PP.H0-H4 per candidate over the shared cis variant grid. `status` records whether the gene had enough shared variants to be testable; `coloc_grade` is strong (PP.H4 >= 0.80), suggestive (0.50-0.80) or none.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `gtex_tissue` |
| `eqtl_dataset` |
| `chr` |
| `gene_start` |
| `gene_end` |
| `window_kb` |
| `nsnp` |
| `PP_H0` |
| `PP_H1` |
| `PP_H2` |
| `PP_H3` |
| `PP_H4` |
| `coloc_grade` |
| `lead_rsid` |
| `lead_p_eqtl` |
| `lead_p_gwas` |
| `status` |

## `day3_coloc_variant_level.csv`

**Rows:** 4,329

The variant-level effect grid that coloc.abf consumed, retained so the posterior can be recomputed without re-querying Ensembl or GTEx.

| Column |
| --- |
| `symbol` |
| `rsid` |
| `chr` |
| `position` |
| `beta_eqtl` |
| `se_eqtl` |
| `p_eqtl` |
| `beta_gwas` |
| `se_gwas` |
| `p_gwas` |

## `day3_eqtl_source.csv`

**Rows:** 50

Which GTEx tissue supplied the dense cis-eQTL grid for each gene, and how many variants that grid held.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `gtex_gencode_id` |
| `n_tissues_tested` |
| `n_variants_best_tissue` |
| `gtex_tissue` |
| `top_rsid` |
| `top_nes` |
| `top_p` |
| `status` |

## `day3_smr_heidi.csv`

**Rows:** 50

SMR Wald ratio of the eQTL effect on the GWAS effect, with the HEIDI heterogeneity test over the cis region. `smr_pass` and `heidi_pass` are the two independent filters.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `top_rsid` |
| `b_eqtl` |
| `se_eqtl` |
| `b_gwas` |
| `se_gwas` |
| `b_smr` |
| `se_smr` |
| `p_smr` |
| `heidi_nsnp` |
| `heidi_Q` |
| `heidi_p` |
| `smr_pass` |
| `heidi_pass` |
| `status` |

## `day4_cell_qc.csv`

**Rows:** 27,600

Per-barcode QC for all 27,600 cells of GSE242382 before filtering, including target-gene counts and CPM.

| Column |
| --- |
| `barcode` |
| `total_counts` |
| `n_genes` |
| `mito_pct` |
| `target_gene_counts` |
| `target_cpm` |
| `cluster` |
| `pass_qc` |

## `day4_cluster_composition.csv`

**Rows:** 5

Cluster sizes after QC.

| Column |
| --- |
| `cluster` |
| `n_cells` |
| `pct` |

## `day4_cross_lineage_replication.csv`

**Rows:** 6

Per-lineage replication summary: Pearson r against the primary contrast, sign agreement, and how much of the primary signature the lineage recovers.

| Column |
| --- |
| `lineage` |
| `cells_high` |
| `cells_low` |
| `genes_compared` |
| `pearson_r_vs_primary` |
| `sign_agreement_pct_all_genes` |
| `signature_genes_with_data` |
| `signature_genes_concordant` |
| `signature_concordance_pct` |
| `lineage_own_signature_n` |
| `primary_recovered_pct` |

## `day4_cross_lineage_signature.csv`

**Rows:** 600

The same 600 genes re-tested inside each of six lineages, with per-lineage log2FC and the number of lineages in which the direction repeats.

| Column |
| --- |
| `symbol` |
| `log2fc_primary` |
| `log2fc_myeloid` |
| `log2fc_fibroblast` |
| `log2fc_mast` |
| `log2fc_melanocyte` |
| `log2fc_tnk` |
| `log2fc_endothelial` |
| `lineages_available` |
| `lineages_concordant` |
| `concordance_pct` |

## `day4_group_balance.csv`

**Rows:** 3

The high vs. low group comparison: cell numbers, total UMI, target CPM, residual difference and lineage composition. The groups differ only in target-gene expression.

| Column |
| --- |
| `group` |
| `n_cells` |
| `mean_total_umi` |
| `median_total_umi` |
| `median_n_genes` |
| `mean_target_cpm` |
| `pct_detected` |
| `mean_mito_pct` |
| `mean_residual` |
| `lineage_mix` |

## `day4_lineage_feasibility.csv`

**Rows:** 9

Detection rate of the target gene in each lineage. This table is the reason the analysis is a state contrast and not a cell-type contrast: detection runs 11.6-20.6 percent everywhere, so the gene marks no lineage.

| Column |
| --- |
| `lineage` |
| `cells_total` |
| `cells_detecting_target` |
| `detection_pct` |
| `mean_target_cpm` |
| `eligible_for_primary_split` |
| `eligible_as_replication` |

## `day4_lineage_filter_audit.csv`

**Rows:** 117

Every gene removed from the signature and why: identity-dominated, depth-confounded or below the detection floor.

| Column |
| --- |
| `symbol` |
| `log2fc_cluster_balanced` |
| `lineage_depth_r2` |
| `action` |
| `reason` |

## `day4_lineage_scores.csv`

**Rows:** 27,282

Per-cell identity panel scores over the eight conjunctival lineages plus a proliferation panel, the dominant-lineage call, and the residualised target-gene expression that defines the state split.

| Column |
| --- |
| `barcode` |
| `cluster` |
| `dominant_lineage` |
| `total_counts` |
| `target_cpm` |
| `target_residual` |
| `target_residual_refit` |
| `score_bplasma` |
| `score_endothelial` |
| `score_epithelial` |
| `score_fibroblast` |
| `score_mast` |
| `score_melanocyte` |
| `score_myeloid` |
| `score_proliferating` |
| `score_tnk` |

## `day4_null_calibration.csv`

**Rows:** 8

The randomly permuted stratifier control: correlation and sign agreement of the null signature, showing that the split procedure itself introduces no direction.

| Column |
| --- |
| `metric` |
| `value` |

## `day4_probe_cells.csv`

**Rows:** 27,282

The full per-cell probe table behind `day4_probe_genes.csv`: QC columns, the eight lineage panel scores, and per-cell CPM for all 20 prioritised candidates. Long, but it is the only file that lets the detection-rate numbers be recomputed from raw counts.

| Column |
| --- |
| `barcode` |
| `cluster` |
| `total_counts` |
| `n_genes` |
| `mito_pct` |
| `score_bplasma` |
| `score_endothelial` |
| `score_epithelial` |
| `score_fibroblast` |
| `score_mast` |
| `score_melanocyte` |
| `score_myeloid` |
| `score_tnk` |
| `cpm_TXNDC12` |
| `cpm_PBXIP1` |
| `cpm_IGSF8` |
| `cpm_ADPRH` |
| `cpm_TLR10` |
| `cpm_TLR1` |
| `cpm_IL7R` |
| `cpm_NNT` |
| `cpm_LMAN2` |
| `cpm_HSD17B8` |
| `cpm_PHF1` |
| `cpm_HSF2` |
| `cpm_NFE2L3` |
| `cpm_OSGIN2` |
| `cpm_SFXN2` |
| `cpm_PRKAB1` |
| `cpm_ATP7B` |
| `cpm_CRK` |
| `cpm_SMAD4` |
| `cpm_UCKL1` |

## `day4_probe_gene_scores.csv`

**Rows:** 20

Candidate-level probe summary: detection rate, mean CPM, identity R-squared and the two identity panels each gene correlates with most strongly.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `coloc_PP_H4` |
| `screen_p` |
| `pct_detected` |
| `mean_cpm` |
| `median_cpm` |
| `lineage_r2` |
| `top_panel_1` |
| `corr_1` |
| `top_panel_2` |
| `corr_2` |

## `day4_pseudobulk_groups.csv`

**Rows:** 18,007

Pseudobulk CPM per gene per group with log2 fold change, cluster-balanced log2FC, cross-cluster agreement, identity R-squared, t statistic and BH FDR over 18,007 genes. The signature is selected from this table.

| Column |
| --- |
| `ensembl_gene_id` |
| `symbol` |
| `locus_group` |
| `cpm_high` |
| `cpm_low` |
| `log2fc` |
| `log2fc_cluster_balanced` |
| `cluster_agreement` |
| `lineage_depth_r2` |
| `t_stat` |
| `p_value` |
| `fdr` |
| `pct_high` |
| `pct_low` |

## `day4_sensitivity_signature.csv`

**Rows:** 590

Half-split stability of the signature across two independent cell halves, reporting which genes hold direction in both.

| Column |
| --- |
| `symbol` |
| `log2fc_primary` |
| `log2fc_half_a` |
| `log2fc_half_b` |
| `in_half_a` |
| `in_half_b` |
| `direction_concordant_both` |

## `day4_signature_dn.csv`

**Rows:** 300

The 300 down-regulated signature genes under the same filters.

| Column |
| --- |
| `ensembl_gene_id` |
| `symbol` |
| `locus_group` |
| `cpm_high` |
| `cpm_low` |
| `log2fc` |
| `log2fc_cluster_balanced` |
| `cluster_agreement` |
| `lineage_depth_r2` |
| `t_stat` |
| `p_value` |
| `fdr` |
| `pct_high` |
| `pct_low` |

## `day4_signature_summary.csv`

**Rows:** 28

Two-column (metric, value) summary of the signature build: cells total, genes tested, filters applied, genes surviving each filter, final signature size.

| Column |
| --- |
| `metric` |
| `value` |

## `day4_signature_up.csv`

**Rows:** 300

The 300 up-regulated signature genes with the filters that admitted them (absolute log2FC > 0.25, cluster agreement >= 0.6, identity R-squared < 0.1, detection >= 5 percent).

| Column |
| --- |
| `ensembl_gene_id` |
| `symbol` |
| `locus_group` |
| `cpm_high` |
| `cpm_low` |
| `log2fc` |
| `log2fc_cluster_balanced` |
| `cluster_agreement` |
| `lineage_depth_r2` |
| `t_stat` |
| `p_value` |
| `fdr` |
| `pct_high` |
| `pct_low` |

## `day4_stratification_design.csv`

**Rows:** 9

How many cells each lineage contributed, how many depth strata were used, and the library-size gap left after matching - the balance audit behind the split.

| Column |
| --- |
| `lineage` |
| `cells_in_universe` |
| `depth_strata` |
| `cells_high` |
| `cells_low` |
| `refit_r2_lineage_depth` |
| `library_size_gap_pct` |
| `note` |

## `day4_target_gene_by_cluster.csv`

**Rows:** 5

Target-gene expression summarised by author cluster.

| Column |
| --- |
| `cluster` |
| `n_cells` |
| `mean_target_cpm` |
| `median_target_cpm` |
| `pct_expressing` |

## `day5_cmap_connectivity.csv`

**Rows:** 46

Per-compound reverse and mimic scores, their difference (reversal specificity), the null z, the overlapping genes and the best matching perturbagen context. `direction` records whether the compound pushes the signature back (reverse) or forward (mimic).

| Column |
| --- |
| `compound` |
| `compound_class` |
| `n_instances` |
| `reverse_best_score` |
| `mimic_best_score` |
| `reversal_specificity` |
| `null_z` |
| `reversal_overlap_genes` |
| `mimic_overlap_genes` |
| `best_rank` |
| `best_cell_line` |
| `best_dose` |
| `best_time` |
| `direction` |

## `day5_cmap_summary.csv`

**Rows:** 17

Headline CMap metrics: landmark coverage of the signature, top score, null comparison, compound counts.

| Column |
| --- |
| `metric` |
| `value` |

## `day5_enrichr_library_status.csv`

**Rows:** 6

Which Enrichr LINCS libraries answered, how many terms they returned, and how long each took.

| Column |
| --- |
| `library` |
| `signature` |
| `status` |
| `n_terms` |
| `seconds` |

## `day5_enrichr_lincs_dn.csv`

**Rows:** 687

Enrichr LINCS results for the down arm.

| Column |
| --- |
| `library` |
| `rank` |
| `term` |
| `p_value` |
| `z_score` |
| `combined_score` |
| `overlap_genes` |
| `fdr` |
| `direction` |

## `day5_enrichr_lincs_up.csv`

**Rows:** 392

Enrichr LINCS results for the up arm - the independent reference-library cross-check.

| Column |
| --- |
| `library` |
| `rank` |
| `term` |
| `p_value` |
| `z_score` |
| `combined_score` |
| `overlap_genes` |
| `fdr` |
| `direction` |

## `day5_l1000cds2_raw.csv`

**Rows:** 100

The raw L1000CDS2 response for both query orientations, retained unmodified as the audit trail behind the connectivity table.

| Column |
| --- |
| `query_direction` |
| `rank` |
| `pert_desc` |
| `pert_id` |
| `cell_id` |
| `pert_dose` |
| `pert_time` |
| `score` |
| `DEGcount` |
| `pubchem_id` |
| `drugbank_id` |
| `overlap_up_in_dn` |
| `overlap_dn_in_up` |

## `day5_null_by_size.csv`

**Rows:** 4

The permutation null at each query size: 20 random landmark draws per size, observed score, null mean/sd/min/max, z and the permutation P.

| Column |
| --- |
| `query_n` |
| `observed_top_score` |
| `null_n` |
| `null_mean` |
| `null_sd` |
| `null_min` |
| `null_max` |
| `n_null_ge_observed` |
| `observed_z` |
| `permutation_p` |
| `null_draws_with_glucocorticoid` |

## `day5_positive_control_check.csv`

**Rows:** 5

Pass/fail record for each internal control: dexamethasone recovery and null separation.

| Column |
| --- |
| `control` |
| `expected` |
| `n_retrieved` |
| `best_rank` |
| `pass` |
| `detail` |
| `note` |

## `day5_query_size_calibration.csv`

**Rows:** 12

The known-answer control: a dexamethasone response signature pushed through the identical pipeline at four query sizes, recording whether glucocorticoids come back and at what rank. This is what fixes the primary query size at 15 genes per arm.

| Column |
| --- |
| `signature` |
| `query_n` |
| `orientation` |
| `compounds_returned` |
| `top_score` |
| `n_glucocorticoid` |
| `best_glucocorticoid_rank` |
| `glucocorticoids_retrieved` |
| `top5` |

## `day6_gene_safety_verdict.csv`

**Rows:** 14

The verdict per gene - safe, monitor or high_risk - with the domain breakdown and the stated reason.

| Column |
| --- |
| `gene` |
| `ensembl_gene_id` |
| `rsid` |
| `screen_p` |
| `coloc_PP_H4` |
| `n_phenotypes_tested` |
| `n_significant` |
| `n_ocular` |
| `n_severe` |
| `top_domains` |
| `safety_verdict` |
| `verdict_reason` |

## `day6_phenotype_universe.csv`

**Rows:** 2,803

The candidate FinnGen R13 endpoint universe before the case-count filter.

| Column |
| --- |
| `dataset` |
| `phenocode` |
| `trait` |
| `category` |
| `ncase` |
| `ncontrol` |
| `sample_size` |

## `day6_phenotypes_scanned.csv`

**Rows:** 1,400

The 1,400 endpoints that survived the filter and were actually scanned for every gene.

| Column |
| --- |
| `dataset` |
| `phenocode` |
| `trait` |
| `category` |
| `ncase` |
| `ncontrol` |
| `sample_size` |

## `day6_phewas_all.csv`

**Rows:** 19,600

Every gene-endpoint test (14 genes x 1,400 endpoints), with disease domain and severity annotation assigned from the endpoint metadata.

| Column |
| --- |
| `gene` |
| `rsid` |
| `dataset` |
| `phenocode` |
| `trait` |
| `category` |
| `domain` |
| `severity` |
| `ncase` |
| `beta` |
| `se` |
| `p` |
| `or` |
| `ea` |
| `nea` |
| `significant` |

## `day6_phewas_hits.csv`

**Rows:** 106

The 106 associations at P < 5e-5, annotated by domain and severity - the input to the three-tier verdict.

| Column |
| --- |
| `gene` |
| `rsid` |
| `dataset` |
| `phenocode` |
| `trait` |
| `category` |
| `domain` |
| `severity` |
| `ncase` |
| `beta` |
| `se` |
| `p` |
| `or` |
| `ea` |
| `nea` |
| `significant` |

## `day7_axis_scores.csv`

**Rows:** 494

The five axis scores mapped to [0, 1] for all 494 confirmatory genes, with the composite score.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `druggable_tier` |
| `score_mr` |
| `score_coloc` |
| `score_singlecell` |
| `score_cmap` |
| `score_safety` |
| `composite_score` |
| `ivw_or` |
| `ivw_p` |
| `coloc_PP_H4` |
| `safety_verdict` |
| `replication_p` |
| `steiger_ok` |

## `day7_candidate_drug_pairs.csv`

**Rows:** 48

Gene-compound pairs formed by crossing the strongest reversal compounds with the top-ranked genes, with the evidence basis for each pairing.

| Column |
| --- |
| `gene` |
| `compound` |
| `compound_class` |
| `reversal_specificity` |
| `reverse_best_score` |
| `mimic_best_score` |
| `null_z` |
| `best_cell_line` |
| `best_dose` |
| `best_time` |
| `evidence_basis` |

## `day7_composite_ranking.csv`

**Rows:** 494

The same table ranked, with `candidate_tier` (Tier 1 >= 0.55, Tier 2, Tier 3). Table 2 of the manuscript is the top of this table.

| Column |
| --- |
| `symbol` |
| `ensembl_gene_id` |
| `druggable_tier` |
| `score_mr` |
| `score_coloc` |
| `score_singlecell` |
| `score_cmap` |
| `score_safety` |
| `composite_score` |
| `ivw_or` |
| `ivw_p` |
| `coloc_PP_H4` |
| `safety_verdict` |
| `replication_p` |
| `steiger_ok` |
| `rank` |
| `candidate_tier` |

## `day7_project_manifest.csv`

**Rows:** 82

Inventory of every project artifact with size, used to verify the delivered package.

| Column |
| --- |
| `kind` |
| `file_name` |
| `relative_path` |
| `size_bytes` |

