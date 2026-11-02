# -*- coding: utf-8 -*-
"""Generate docs/DATA_DICTIONARY.md from the archived result tables."""
import csv
import glob
import io
import os

NOTES = {
    "day1_druggable_pool.csv":
        "One row per druggable gene. `druggable_tier` is T1 (clinical-stage annotation plus a recorded "
        "drug interaction), T2 (clinical annotation only) or T3 (literature or experimental evidence only); "
        "`has_cis_eqtl` marks whether an OpenGWAS `eqtl-a-ENSG*` dataset exists and therefore whether the "
        "gene can be instrumented at all.",
    "day1_hgnc_protein_coding.csv":
        "HGNC approved, protein-coding records carrying an Ensembl gene ID. This is the 19,255-gene "
        "starting universe.",
    "day1_pool_summary.csv":
        "Tier-level counts and medians behind Table 1 of the manuscript.",
    "day1_druggable_sources.csv":
        "Per-source-database gene counts, with the evidence class each database contributes.",
    "day2_screen_lead_snp_mr.csv":
        "Stage A screening: one lead cis-eQTL instrument per gene, Wald ratio against the FinnGen "
        "allergic-conjunctivitis endpoint. `wald_fdr` is Benjamini-Hochberg across the screened genes; "
        "`f_statistic` is the first-stage strength of the lead instrument.",
    "day2_instruments_lead.csv":
        "The lead instrument for every screened gene, with its exposure and outcome effect pair, kept so "
        "the screen can be audited SNP by SNP.",
    "day2_mr_confirmatory.csv":
        "Stage B: for every gene reaching screening P < 0.05, all independent cis instruments under six "
        "estimators, plus Cochran Q, the MR-Egger intercept, leave-one-out support, replication in an "
        "independent conjunctivitis endpoint, the Steiger directionality test and BH FDR.",
    "day2_mr_heterogeneity.csv":
        "Per-SNP Wald ratios for the confirmatory genes - the input to Cochran Q.",
    "day2_mr_leave_one_out.csv":
        "Leave-one-out IVW estimates for the genes where a single instrument could have driven the result.",
    "day2_mr_steiger.csv":
        "Steiger r-squared for exposure and outcome, with the direction verdict.",
    "day2_screen_summary.csv":
        "Headline screen counts: genes screened, counts at each significance threshold, median F.",
    "day3_coloc_results.csv":
        "coloc.abf posterior probabilities PP.H0-H4 per candidate over the shared cis variant grid. "
        "`status` records whether the gene had enough shared variants to be testable; `coloc_grade` is "
        "strong (PP.H4 >= 0.80), suggestive (0.50-0.80) or none.",
    "day3_coloc_variant_level.csv":
        "The variant-level effect grid that coloc.abf consumed, retained so the posterior can be "
        "recomputed without re-querying Ensembl or GTEx.",
    "day3_eqtl_source.csv":
        "Which GTEx tissue supplied the dense cis-eQTL grid for each gene, and how many variants that "
        "grid held.",
    "day3_smr_heidi.csv":
        "SMR Wald ratio of the eQTL effect on the GWAS effect, with the HEIDI heterogeneity test over the "
        "cis region. `smr_pass` and `heidi_pass` are the two independent filters.",
    "day3_candidate_evidence.csv":
        "The merged candidate table: MR, coloc and SMR/HEIDI evidence side by side, with "
        "`n_supporting_axes` and the `prioritised` flag that decides who enters Day 4.",
    "day4_cell_qc.csv":
        "Per-barcode QC for all 27,600 cells of GSE242382 before filtering, including target-gene counts "
        "and CPM.",
    "day4_lineage_scores.csv":
        "Per-cell identity panel scores over the eight conjunctival lineages plus a proliferation panel, "
        "the dominant-lineage call, and the residualised target-gene expression that defines the state "
        "split.",
    "day4_lineage_feasibility.csv":
        "Detection rate of the target gene in each lineage. This table is the reason the analysis is a "
        "state contrast and not a cell-type contrast: detection runs 11.6-20.6 percent everywhere, so the "
        "gene marks no lineage.",
    "day4_stratification_design.csv":
        "How many cells each lineage contributed, how many depth strata were used, and the library-size "
        "gap left after matching - the balance audit behind the split.",
    "day4_group_balance.csv":
        "The high vs. low group comparison: cell numbers, total UMI, target CPM, residual difference and "
        "lineage composition. The groups differ only in target-gene expression.",
    "day4_pseudobulk_groups.csv":
        "Pseudobulk CPM per gene per group with log2 fold change, cluster-balanced log2FC, cross-cluster "
        "agreement, identity R-squared, t statistic and BH FDR over 18,007 genes. The signature is "
        "selected from this table.",
    "day4_signature_up.csv":
        "The 300 up-regulated signature genes with the filters that admitted them (absolute log2FC > "
        "0.25, cluster agreement >= 0.6, identity R-squared < 0.1, detection >= 5 percent).",
    "day4_signature_dn.csv":
        "The 300 down-regulated signature genes under the same filters.",
    "day4_lineage_filter_audit.csv":
        "Every gene removed from the signature and why: identity-dominated, depth-confounded or below the "
        "detection floor.",
    "day4_cross_lineage_signature.csv":
        "The same 600 genes re-tested inside each of six lineages, with per-lineage log2FC and the number "
        "of lineages in which the direction repeats.",
    "day4_cross_lineage_replication.csv":
        "Per-lineage replication summary: Pearson r against the primary contrast, sign agreement, and how "
        "much of the primary signature the lineage recovers.",
    "day4_sensitivity_signature.csv":
        "Half-split stability of the signature across two independent cell halves, reporting which genes "
        "hold direction in both.",
    "day4_null_calibration.csv":
        "The randomly permuted stratifier control: correlation and sign agreement of the null signature, "
        "showing that the split procedure itself introduces no direction.",
    "day4_probe_genes.csv":
        "Per-cell CPM for all 20 prioritised candidate genes - the raw material for the detection-rate "
        "figure and for choosing the stratifier.",
    "day4_probe_gene_scores.csv":
        "Candidate-level probe summary: detection rate, mean CPM, identity R-squared and the two identity "
        "panels each gene correlates with most strongly.",
    "day4_probe_cells.csv":
        "The full per-cell probe table behind `day4_probe_genes.csv`: QC columns, the eight lineage "
        "panel scores, and per-cell CPM for all 20 prioritised candidates. Long, but it is the only file "
        "that lets the detection-rate numbers be recomputed from raw counts.",
    "day4_signature_summary.csv":
        "Two-column (metric, value) summary of the signature build: cells total, genes tested, filters "
        "applied, genes surviving each filter, final signature size.",
    "day4_target_gene_by_cluster.csv":
        "Target-gene expression summarised by author cluster.",
    "day4_cluster_composition.csv":
        "Cluster sizes after QC.",
    "day5_query_size_calibration.csv":
        "The known-answer control: a dexamethasone response signature pushed through the identical "
        "pipeline at four query sizes, recording whether glucocorticoids come back and at what rank. This "
        "is what fixes the primary query size at 15 genes per arm.",
    "day5_null_by_size.csv":
        "The permutation null at each query size: 20 random landmark draws per size, observed score, null "
        "mean/sd/min/max, z and the permutation P.",
    "day5_cmap_connectivity.csv":
        "Per-compound reverse and mimic scores, their difference (reversal specificity), the null z, the "
        "overlapping genes and the best matching perturbagen context. `direction` records whether the "
        "compound pushes the signature back (reverse) or forward (mimic).",
    "day5_cmap_summary.csv":
        "Headline CMap metrics: landmark coverage of the signature, top score, null comparison, compound "
        "counts.",
    "day5_l1000cds2_raw.csv":
        "The raw L1000CDS2 response for both query orientations, retained unmodified as the audit trail "
        "behind the connectivity table.",
    "day5_positive_control_check.csv":
        "Pass/fail record for each internal control: dexamethasone recovery and null separation.",
    "day5_enrichr_library_status.csv":
        "Which Enrichr LINCS libraries answered, how many terms they returned, and how long each took.",
    "day5_enrichr_lincs_up.csv":
        "Enrichr LINCS results for the up arm - the independent reference-library cross-check.",
    "day5_enrichr_lincs_dn.csv":
        "Enrichr LINCS results for the down arm.",
    "day6_phenotype_universe.csv":
        "The candidate FinnGen R13 endpoint universe before the case-count filter.",
    "day6_phenotypes_scanned.csv":
        "The 1,400 endpoints that survived the filter and were actually scanned for every gene.",
    "day6_phewas_all.csv":
        "Every gene-endpoint test (14 genes x 1,400 endpoints), with disease domain and severity "
        "annotation assigned from the endpoint metadata.",
    "day6_phewas_hits.csv":
        "The 106 associations at P < 5e-5, annotated by domain and severity - the input to the three-tier "
        "verdict.",
    "day6_gene_safety_verdict.csv":
        "The verdict per gene - safe, monitor or high_risk - with the domain breakdown and the stated "
        "reason.",
    "day7_axis_scores.csv":
        "The five axis scores mapped to [0, 1] for all 494 confirmatory genes, with the composite score.",
    "day7_composite_ranking.csv":
        "The same table ranked, with `candidate_tier` (Tier 1 >= 0.55, Tier 2, Tier 3). Table 2 of the "
        "manuscript is the top of this table.",
    "day7_candidate_drug_pairs.csv":
        "Gene-compound pairs formed by crossing the strongest reversal compounds with the top-ranked "
        "genes, with the evidence basis for each pairing.",
    "day7_project_manifest.csv":
        "Inventory of every project artifact with size, used to verify the delivered package.",
}


def main() -> None:
    rows = []
    for f in sorted(glob.glob("tables/*.csv")):
        name = os.path.basename(f)
        with open(f, encoding="utf-8-sig", newline="") as fh:
            r = csv.reader(fh)
            hdr = next(r, [])
            n = sum(1 for _ in r)
        rows.append((name, n, hdr))

    out = io.StringIO()
    out.write("# Data dictionary\n\n")
    out.write(
        "Every result table the manuscript draws on. Row counts are from the archived run; "
        "column names are verbatim. All CSVs are UTF-8 with BOM so Excel opens them correctly, "
        "and all floats are written at six significant figures.\n\n"
    )
    out.write("Path prefix: `tables/`.\n\n---\n\n")
    for name, n, hdr in rows:
        out.write("## `%s`\n\n" % name)
        out.write("**Rows:** %s\n\n" % format(n, ","))
        note = NOTES.get(name)
        if note:
            out.write(note + "\n\n")
        out.write("| Column |\n| --- |\n")
        for c in hdr:
            out.write("| `%s` |\n" % c)
        out.write("\n")

    text = out.getvalue()
    with open("docs/DATA_DICTIONARY.md", "w", encoding="utf-8") as fh:
        fh.write(text)
    print("wrote", len(text), "chars,", len(rows), "tables")


if __name__ == "__main__":
    main()
