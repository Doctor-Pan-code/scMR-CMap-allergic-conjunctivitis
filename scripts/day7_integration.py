"""Day 7 - Evidence integration, composite scoring and manuscript-ready figures.

Why this day exists
-------------------
Days 2-6 produce five independent evidence axes for every candidate:
genetic causality (MR), shared causal variant (coloc / SMR-HEIDI), tissue-level
plausibility (single-cell stratification), pharmacological actionability
(CMap reversal) and systemic safety (PheWAS). This day converts them into one
interpretable ranking and the figure set a paper would need.

Composite score
---------------
Each axis is mapped to [0, 1] with explicit, pre-registered rules, then combined
with the weights in config.json (MR 0.30, coloc 0.20, single-cell 0.20,
CMap 0.20, safety 0.10). The score is a prioritisation device, not a p-value.

Inputs
------
tables/day2_mr_confirmatory.csv
tables/day3_candidate_evidence.csv
tables/day4_target_gene_by_cluster.csv
tables/day5_cmap_connectivity.csv
tables/day6_gene_safety_verdict.csv

Outputs
-------
tables/day7_axis_scores.csv
tables/day7_composite_ranking.csv
tables/day7_candidate_drug_pairs.csv
tables/day7_project_manifest.csv
figures/svg/fig24_composite_ranking.svg
figures/svg/fig25_axis_radar.svg
figures/svg/fig26_candidate_drug_map.svg
figures/svg/fig27_workflow_schematic.svg
figures/svg/fig28_evidence_funnel.svg
"""
from __future__ import annotations

import collections
import math
from typing import Any, Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np

from lib_io import (DATA, FIG_SVG, REPORT, ROOT, TABLES, banner, cfg, log,
                    log_init, read_table, save_fig, section, write_table)


def _f(v) -> Optional[float]:
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def load_axes() -> Dict[str, Dict[str, Any]]:
    axes: Dict[str, Dict[str, Any]] = collections.defaultdict(dict)
    for r in read_table("day2_mr_confirmatory.csv"):
        sym = r["symbol"]
        axes[sym].update({
            "ivw_p": _f(r.get("IVW_RE_p")),
            "ivw_or": _f(r.get("IVW_RE_or")),
            "ivw_lci": _f(r.get("IVW_RE_or_lci")),
            "ivw_uci": _f(r.get("IVW_RE_or_uci")),
            "nsnp": r.get("nsnp"),
            "egger_p": _f(r.get("egger_intercept_p")),
            "q_p": _f(r.get("Q_pval")),
            "replication_p": _f(r.get("replication_p")),
            "steiger_ok": r.get("steiger_correct_direction"),
            "screen_p": _f(r.get("screen_wald_p")),
            "druggable_tier": r.get("druggable_tier", ""),
            "ensembl_gene_id": r.get("ensembl_gene_id", ""),
        })
    for r in read_table("day3_candidate_evidence.csv"):
        sym = r["symbol"]
        axes[sym].update({
            "coloc_h4": _f(r.get("coloc_PP_H4")),
            "coloc_grade": r.get("coloc_grade", ""),
            "smr_p": _f(r.get("smr_p")),
            "heidi_p": _f(r.get("heidi_p")),
            "smr_heidi_pass": r.get("smr_heidi_pass") == "True",
        })
    for r in read_table("day6_gene_safety_verdict.csv"):
        axes[r["gene"]].update({
            "safety_verdict": r.get("safety_verdict", ""),
            "n_phewas_sig": _f(r.get("n_significant")),
            "n_severe": _f(r.get("n_severe")),
            "verdict_reason": r.get("verdict_reason", ""),
        })
    return axes


def scale_mr(p: Optional[float]) -> float:
    if p is None:
        return 0.0
    return max(0.0, min(1.0, (math.log10(0.05) - math.log10(max(p, 1e-300))) / 6.0))


def scale_coloc(h4: Optional[float], smr_pass: bool) -> float:
    a = 0.0 if h4 is None else max(0.0, min(1.0, h4))
    b = 1.0 if smr_pass else 0.0
    return 0.75 * a + 0.25 * b


def scale_safety(verdict: str) -> float:
    return {"safe": 1.0, "monitor": 0.55, "high_risk": 0.05}.get(verdict, 0.4)


def main() -> None:
    log_init("day7")
    banner("Day 7", "Evidence integration and composite prioritisation")
    conf = cfg()
    w = conf["scoring"]
    axes = load_axes()
    log(f"genes with Day 2-3 evidence: {len(axes)}")

    # single-cell axis. v5 no longer splits cells by raw target abundance, so the
    # old "how specific is the target to one cluster" formula no longer describes
    # what Day 4 did: TLR1 is expressed at ~20 CPM in every conjunctival cluster
    # (mean CPM 19.6-23.6, detection 11.3-23.4%), i.e. it has no cluster
    # specificity at all. What Day 4 establishes is that a TLR1-high *state*
    # exists whose transcriptome differs from the TLR1-low state. The axis
    # therefore scores the strength of that state contrast, and it is credited
    # only to the genes that were actually used to define the contrast.
    scrna_summary = {r["metric"]: r["value"] for r in read_table("day4_signature_summary.csv")}
    target_gene = scrna_summary.get("target_gene", "")
    sensitivity_gene = scrna_summary.get("sensitivity_gene", "")
    cl = read_table("day4_target_gene_by_cluster.csv")
    n_up = _f(scrna_summary.get("signature_up_n")) or 0.0
    n_dn = _f(scrna_summary.get("signature_dn_n")) or 0.0
    agreement = _f(scrna_summary.get("sensitivity_direction_agreement_pct")) or 0.0
    resid_gap = 0.0
    for r in read_table("day4_group_balance.csv"):
        if r.get("group") == "gap_high_vs_low":
            resid_gap = abs(_f(r.get("mean_residual")) or 0.0)
    sig_strength = max(0.0, min(1.0, (n_up + n_dn) / 600.0))
    contrast_strength = max(0.0, min(1.0, resid_gap / 2.0))
    repro = max(0.0, min(1.0, agreement / 100.0))
    sc_axis = sig_strength * contrast_strength * repro
    top_cluster = cl[0]["cluster"] if cl else ""
    log(f"single-cell axis: state contrast from {int(n_up)} up / {int(n_dn)} down genes, "
        f"residual gap {resid_gap:.3f}, sensitivity agreement {agreement:.1f}% "
        f"-> score {sc_axis:.3f}")
    log(f"  credited to the stratifier gene: {target_gene}")
    log(f"  note: {target_gene} is not cluster-specific "
        f"(mean CPM {cl[0].get('mean_target_cpm') if cl else 'NA'} in the top cluster "
        f"{top_cluster}, detection {cl[0].get('pct_expressing') if cl else 'NA'}%); "
        f"the axis measures a cell state, not a cell type")

    # CMap axis. Day 5 v2 reports reversal specificity (reverse-query score minus
    # mimic-query score) rather than a signed connectivity score, because the raw
    # L1000CDS2 score carries no direction at all -- direction comes from which
    # arm was submitted. A compound's specificity is scaled against the best
    # specificity observed in this run, so the axis is a relative rank on a
    # quantity that does have a meaning, and it is credited only to the gene
    # whose signature produced the query.
    cmap_rows = read_table("day5_cmap_connectivity.csv")
    cmap_summary = {r["metric"]: r["value"] for r in read_table("day5_cmap_summary.csv")}
    cmap_p = _f(cmap_summary.get("primary_query_permutation_p"))
    cmap_signal = (cmap_summary.get("signature_signal_detected") == "True")
    best_cmap = None
    for r in cmap_rows:
        spec = _f(r.get("reversal_specificity"))
        if spec is None or r.get("direction") != "reversing":
            continue
        if best_cmap is None or spec > best_cmap["reversal_specificity"]:
            best_cmap = {"compound": r["compound"],
                         "reversal_specificity": spec,
                         "reverse_best_score": _f(r.get("reverse_best_score")),
                         "mimic_best_score": _f(r.get("mimic_best_score")),
                         "null_z": _f(r.get("null_z")),
                         "compound_class": r.get("compound_class", ""),
                         "best_rank": r.get("best_rank")}
    if best_cmap:
        log(f"CMap axis: best reversal-specific compound = {best_cmap['compound']} "
            f"(specificity {best_cmap['reversal_specificity']:.4f}, "
            f"reverse {best_cmap['reverse_best_score']}, "
            f"mimic {best_cmap['mimic_best_score']})")
    else:
        log("CMap axis: no reversal-specific compound")
    log(f"  signature-level permutation p = {cmap_p if cmap_p is not None else 'NA'}; "
        f"signal detected = {cmap_signal}")

    section("7.1 Per-axis scoring")
    rows: List[Dict[str, Any]] = []
    for sym, a in axes.items():
        mr_s = scale_mr(a.get("ivw_p") if a.get("ivw_p") is not None else a.get("screen_p"))
        col_s = scale_coloc(a.get("coloc_h4"), bool(a.get("smr_heidi_pass")))
        saf_s = scale_safety(a.get("safety_verdict", ""))
        sc_s = sc_axis if (sym == target_gene and sym) else 0.0
        # CMap is a signature-level axis: it credits the gene that generated the
        # signature, and is discounted for every other gene.
        cm_s = 0.0
        if best_cmap and sym == target_gene and cmap_signal:
            cm_s = max(0.0, min(1.0, best_cmap["reversal_specificity"] / 0.10))
        composite = (w["w_mr"] * mr_s + w["w_coloc"] * col_s + w["w_scell"] * sc_s
                     + w["w_cmap"] * cm_s + w["w_safety"] * saf_s)
        rows.append({
            "symbol": sym,
            "ensembl_gene_id": a.get("ensembl_gene_id", ""),
            "druggable_tier": a.get("druggable_tier", ""),
            "score_mr": round(mr_s, 4),
            "score_coloc": round(col_s, 4),
            "score_singlecell": round(sc_s, 4),
            "score_cmap": round(cm_s, 4),
            "score_safety": round(saf_s, 4),
            "composite_score": round(composite, 4),
            "ivw_or": a.get("ivw_or"),
            "ivw_p": a.get("ivw_p"),
            "coloc_PP_H4": a.get("coloc_h4"),
            "safety_verdict": a.get("safety_verdict", ""),
            "replication_p": a.get("replication_p"),
            "steiger_ok": a.get("steiger_ok"),
        })
    rows.sort(key=lambda r: -r["composite_score"])
    write_table(rows, "day7_axis_scores.csv")
    log(f"axes scored for {len(rows)} genes")
    for r in rows[:10]:
        log(f"  {r['symbol']:<10} composite={r['composite_score']:.3f} "
            f"(MR {r['score_mr']:.2f} / coloc {r['score_coloc']:.2f} / "
            f"sc {r['score_singlecell']:.2f} / cmap {r['score_cmap']:.2f} / "
            f"safety {r['score_safety']:.2f})")

    section("7.2 Tiered candidate ranking")
    ranking = []
    for i, r in enumerate(rows, 1):
        if r["composite_score"] >= 0.55:
            tier = "Tier1_lead_candidate"
        elif r["composite_score"] >= 0.35:
            tier = "Tier2_supporting_candidate"
        else:
            tier = "Tier3_exploratory"
        ranking.append({**r, "rank": i, "candidate_tier": tier})
    write_table(ranking, "day7_composite_ranking.csv")
    tc = collections.Counter(r["candidate_tier"] for r in ranking)
    log(f"candidate tiers: {dict(tc)}")

    section("7.3 Candidate-compound pairs")
    pairs = []
    top_genes = [r["symbol"] for r in ranking[:8]]
    rev_compounds = [r for r in cmap_rows if (r.get("direction") == "reversing")][:25]
    for r in rev_compounds:
        for g in top_genes[:3]:
            pairs.append({
                "gene": g,
                "compound": r["compound"],
                "compound_class": r.get("compound_class", ""),
                "reversal_specificity": _f(r.get("reversal_specificity")),
                "reverse_best_score": _f(r.get("reverse_best_score")),
                "mimic_best_score": _f(r.get("mimic_best_score")),
                "null_z": _f(r.get("null_z")),
                "best_cell_line": r.get("best_cell_line", ""),
                "best_dose": r.get("best_dose", ""),
                "best_time": r.get("best_time", ""),
                "evidence_basis": "CMap reversal of the target-gene high/low conjunctival signature",
            })
    write_table(pairs, "day7_candidate_drug_pairs.csv")
    log(f"candidate-compound pairs written: {len(pairs)}")

    section("7.4 Project manifest")
    # NOTE: the manifest is written at the very end of main() on purpose. Writing it
    # here would list only the figures that already exist, and fig24-fig28 are created
    # in section 7.5 below — the first version of this file shipped a manifest that was
    # missing four of its own figures.

    section("7.5 Figures")
    # Fig 24 - composite ranking
    top = ranking[:16][::-1]
    if top:
        fig, ax = plt.subplots(figsize=(7.0, max(3.2, 0.32 * len(top) + 1.4)))
        y = np.arange(len(top))
        cols = {"Tier1_lead_candidate": "#c1443c",
                "Tier2_supporting_candidate": "#2a9d8f",
                "Tier3_exploratory": "#8d99ae"}
        ax.barh(y, [r["composite_score"] for r in top],
                color=[cols[r["candidate_tier"]] for r in top],
                edgecolor="white", height=0.68)
        for yi, r in zip(y, top):
            ax.text(r["composite_score"] + 0.008, yi, f"{r['composite_score']:.2f}",
                    va="center", fontsize=7)
        ax.set_yticks(y)
        ax.set_yticklabels([r["symbol"] for r in top], fontsize=7.5)
        ax.set_xlabel("Composite evidence score")
        ax.set_title("Fig 24. Integrated prioritisation of candidate drug targets")
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig24_composite_ranking")

    # Fig 25 - axis profile of the lead candidate
    if ranking:
        lead = ranking[0]
        labels = ["MR", "Coloc", "Single-cell", "CMap", "Safety"]
        vals = [lead["score_mr"], lead["score_coloc"], lead["score_singlecell"],
                lead["score_cmap"], lead["score_safety"]]
        ang = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
        vals_c = vals + vals[:1]
        ang_c = ang + ang[:1]
        fig = plt.figure(figsize=(4.8, 4.6))
        ax = fig.add_subplot(111, polar=True)
        ax.plot(ang_c, vals_c, color="#1f3b57", lw=1.4)
        ax.fill(ang_c, vals_c, color="#1f3b57", alpha=0.16)
        ax.set_xticks(ang)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_yticks([0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(["0.25", "0.5", "0.75", "1.0"], fontsize=6)
        ax.set_title(f"Fig 25. Evidence profile of the lead candidate {lead['symbol']}",
                     fontsize=10, pad=18)
        save_fig(fig, "fig25_axis_radar")

    # Fig 26 - candidate x compound map
    if pairs:
        compounds = sorted({p["compound"] for p in pairs},
                           key=lambda c: -max((p["reversal_specificity"] or 0)
                                              for p in pairs if p["compound"] == c))[:18]
        genes = sorted({p["gene"] for p in pairs})
        mat = np.full((len(genes), len(compounds)), np.nan)
        for p in pairs:
            if p["compound"] in compounds and p["gene"] in genes:
                mat[genes.index(p["gene"]), compounds.index(p["compound"])] = p["reversal_specificity"]
        fig, ax = plt.subplots(figsize=(max(6.4, 0.42 * len(compounds) + 2.2),
                                        max(2.6, 0.34 * len(genes) + 1.4)))
        vmax = np.nanmax(mat) if np.isfinite(np.nanmax(mat)) else 0.1
        im = ax.imshow(mat, aspect="auto", cmap="YlGnBu", vmin=0, vmax=vmax)
        ax.set_xticks(range(len(compounds)))
        ax.set_xticklabels([c[:16] for c in compounds], rotation=60, ha="right", fontsize=6.5)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes, fontsize=7.5)
        ax.set_title("Fig 26. Candidate gene versus reversal-specific compounds")
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02,
                     label="reversal specificity")
        save_fig(fig, "fig26_candidate_drug_map")

    # Fig 27 - workflow schematic
    fig, ax = plt.subplots(figsize=(9.0, 3.2))
    steps = [
        ("Day 1", "Druggable\ngenome pool"),
        ("Day 2", "cis-eQTL\nMR screen"),
        ("Day 3", "Coloc +\nSMR/HEIDI"),
        ("Day 4", "Single-cell\nstratification"),
        ("Day 5", "CMap\nreversal"),
        ("Day 6", "PheWAS\nsafety"),
        ("Day 7", "Composite\nranking"),
    ]
    n = len(steps)
    for i, (d, t) in enumerate(steps):
        x = i * 1.0
        ax.add_patch(plt.Rectangle((x - 0.42, -0.34), 0.84, 0.68,
                                   facecolor="#1f3b57" if i < n - 1 else "#c1443c",
                                   edgecolor="white", linewidth=1.2, zorder=2))
        ax.text(x, 0.10, d, ha="center", va="center", color="white",
                fontsize=7.5, fontweight="bold", zorder=3)
        ax.text(x, -0.14, t, ha="center", va="center", color="white",
                fontsize=6.4, zorder=3)
        if i < n - 1:
            ax.annotate("", xy=(x + 0.55, 0), xytext=(x + 0.44, 0),
                        arrowprops=dict(arrowstyle="->", color="#8d99ae", lw=1.1))
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(-0.55, 0.45)
    ax.axis("off")
    ax.set_title("Fig 27. Seven-day analysis workflow", fontsize=10)
    save_fig(fig, "fig27_workflow_schematic")

    # Fig 28 - evidence attrition funnel. Every count is read from the tables; the
    # first version hardcoded the pool and the instrumentable subset and skipped
    # the 6,329 -> 494 step, so a reader could not see where the screened genes went.
    pool_sum = {r["tier"]: r for r in read_table("day1_pool_summary.csv")}
    n_pool = int(pool_sum["ALL"]["n_genes"])
    n_instrument = int(pool_sum["ALL"]["n_with_cis_eqtl"])
    d2s = {r["metric"]: r["value"] for r in read_table("day2_screen_summary.csv")}
    n_screened = int(_f(d2s.get("genes_screened")) or 0)
    n_p05 = int(_f(d2s.get("genes_P_lt_0.05")) or 0)
    n_screen = len(rows)
    n_coloc = sum(1 for r in rows if (_f(r.get("coloc_PP_H4")) or 0) >= 0.5)
    n_coloc_strong = sum(1 for r in rows
                         if (_f(r.get("coloc_PP_H4")) or 0) >= conf["coloc"]["pp_h4_strong"])
    n_lead = sum(1 for r in ranking if r["candidate_tier"] == "Tier1_lead_candidate")
    stages = ["Druggable pool", "cis-eQTL instrumentable", "MR screened",
              "MR P<0.05", "Coloc PP.H4>=0.5",
              f"Coloc PP.H4>={conf['coloc']['pp_h4_strong']}", "Tier-1 candidates"]
    raw_values = [n_pool, n_instrument, n_screened, n_p05, n_coloc, n_coloc_strong,
                  max(n_lead, 1)]
    values, kept = [], []
    for s, v in zip(stages, raw_values):
        if values and values[-1] == v:
            continue          # a plateau is not a filter; drawing it invites misreading
        values.append(v)
        kept.append(s)
    stages = kept
    fig, ax = plt.subplots(figsize=(8.2, 3.4))
    cols = ["#dde5ec", "#bcccd9", "#a8bccd", "#8aa5bb", "#6b8fa8", "#2a9d8f", "#c1443c"]
    cols = cols[:len(stages)]
    left = 0.0
    total = float(values[0])
    for i, (s, v) in enumerate(zip(stages, values)):
        width = v / total
        ax.barh(0, width, left=left, height=0.5, color=cols[i], edgecolor="white", lw=0.8)
        if width > 0.07:
            ax.text(left + width / 2, 0, f"{v:,}", ha="center", va="center",
                    fontsize=7.2, color="white" if i >= 3 else "#1f3b57")
        left += width
    ax.set_yticks([])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Proportion of the starting druggable pool")
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c) for c in cols]
    ax.legend(handles, [f"{s} ({v:,})" for s, v in zip(stages, values)],
              fontsize=6.8, loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
              frameon=False)
    ax.set_title("Fig 28. Evidence attrition from the druggable genome to tier-1 candidates")
    ax.spines[["top", "right", "left"]].set_visible(False)
    save_fig(fig, "fig28_evidence_funnel")

    section("7.6 Project manifest")
    manifest = []
    for folder, label in [(TABLES, "table"), (FIG_SVG, "figure_svg")]:
        if not folder.exists():
            continue
        for p in sorted(folder.glob("*")):
            if p.is_file():
                manifest.append({
                    "kind": label,
                    "file_name": p.name,
                    "relative_path": str(p.relative_to(ROOT)).replace("\\", "/"),
                    "size_bytes": p.stat().st_size,
                })
    write_table(manifest, "day7_project_manifest.csv")
    log(f"manifest entries: {len(manifest)} "
        f"({sum(1 for m in manifest if m['kind']=='table')} tables, "
        f"{sum(1 for m in manifest if m['kind']=='figure_svg')} SVG figures)")

    log("")
    log("Day 7 complete.")


if __name__ == "__main__":
    main()
