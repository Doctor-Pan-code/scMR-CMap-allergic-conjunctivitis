"""Day 3 - colocalization (coloc.abf) and SMR/HEIDI for the Day 2 candidate genes.

Dense cis-eQTL reference
------------------------
OpenGWAS eqtl-a datasets only store the SNPs that already reached significance,
so they cannot supply the variant-level effect grid that colocalization needs.
This script therefore rebuilds the cis-eQTL grid from the GTEx v10
single-tissue cis-eQTL release (all tested SNP-gene pairs, ~1e3 variants per
gene), pairs it with the OpenGWAS outcome GWAS on rsID, and runs:
  * coloc.abf on the shared cis variants (PP.H0-H4)
  * SMR Wald ratio at the top eQTL variant
  * HEIDI heterogeneity over the cis region, restricted by LD with the top SNP
"""
from __future__ import annotations

import math
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import lib_api as api
import lib_stats as st
from lib_io import (DATA, ROOT, banner, cache_get, cache_put, cfg, log, log_init,
                    read_table, save_fig, section, write_table)

MAX_REGION_VARIANTS = 120
ASSOC_BATCH = 60
MAX_CANDIDATES = 50
LD_PAIRWISE_MAX = 10
GTEX_BASE = "https://gtexportal.org/api/v2"
GTEX_TISSUE_PREF = [
    "Whole_Blood", "Cells_EBV-transformed_lymphocytes", "Skin_Not_Sun_Exposed_Suprapubic",
    "Skin_Sun_Exposed_Lower_leg", "Esophagus_Mucosa", "Lung", "Spleen",
    "Minor_Salivary_Gland", "Thyroid", "Adipose_Subcutaneous",
]

COLOC_COLUMNS = [
    "symbol", "ensembl_gene_id", "gtex_tissue", "eqtl_dataset", "chr", "gene_start",
    "gene_end", "window_kb", "nsnp", "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4",
    "coloc_grade", "lead_rsid", "lead_p_eqtl", "lead_p_gwas", "status",
]
VAR_COLUMNS = ["symbol", "rsid", "chr", "position", "beta_eqtl", "se_eqtl", "p_eqtl",
               "beta_gwas", "se_gwas", "p_gwas"]
SMR_COLUMNS = ["symbol", "ensembl_gene_id", "top_rsid", "b_eqtl", "se_eqtl", "b_gwas",
               "se_gwas", "b_smr", "se_smr", "p_smr", "heidi_nsnp", "heidi_Q",
               "heidi_p", "smr_pass", "heidi_pass", "status"]
EVIDENCE_COLUMNS = ["symbol", "ensembl_gene_id", "druggable_tier", "screen_p", "ivw_p",
                    "ivw_or", "nsnp_mr", "coloc_PP_H4", "coloc_grade", "smr_p",
                    "heidi_p", "smr_heidi_pass", "n_supporting_axes", "prioritised"]
SOURCE_COLUMNS = ["symbol", "ensembl_gene_id", "gtex_gencode_id", "n_tissues_tested",
                  "n_variants_best_tissue", "gtex_tissue", "top_rsid", "top_nes",
                  "top_p", "status"]


def _f(v) -> Optional[float]:
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _i(v) -> Optional[int]:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _z_from_p(p: float) -> Optional[float]:
    """Two-sided z score from a p value, robust for extremely small p."""
    if p is None or not math.isfinite(p) or p <= 0:
        return None
    if p >= 1.0:
        return 0.0
    half = p / 2.0
    if half < 1e-300:
        return None
    return abs(st.norm_ppf(half))


# ---------------------------------------------------------------------------
# candidate list from Day 2
# ---------------------------------------------------------------------------


def candidates() -> List[Dict[str, Any]]:
    rows = read_table("day2_mr_confirmatory.csv")
    out = []
    for r in rows:
        p = _f(r.get("screen_wald_p"))
        if p is None:
            continue
        out.append({
            "symbol": r["symbol"],
            "ensembl_gene_id": r["ensembl_gene_id"],
            "eqtl_dataset": r["eqtl_dataset"],
            "lead_rsid": (r.get("snps") or "").split(";")[0],
            "screen_p": p,
            "ivw_p": _f(r.get("IVW_RE_p")),
            "ivw_or": _f(r.get("IVW_RE_or")),
            "nsnp": _i(r.get("nsnp")),
            "druggable_tier": r.get("druggable_tier", ""),
        })
    out.sort(key=lambda r: r["screen_p"])
    return out[:MAX_CANDIDATES]


# ---------------------------------------------------------------------------
# GTEx cis-eQTL retrieval
# ---------------------------------------------------------------------------


def gtex_gencode(symbols: Sequence[str]) -> Dict[str, str]:
    """Batch-resolve gene symbols to versioned GENCODE ids under GENCODE v39."""
    key = "gtex_gencode_v39"
    cache = cache_get(key) or {}
    todo = [s for s in symbols if s not in cache]
    if todo:
        for i in range(0, len(todo), 20):
            chunk = todo[i:i + 20]
            r = api._get(f"{GTEX_BASE}/reference/gene",
                         params={"geneId": list(chunk), "gencodeVersion": "v39",
                                 "itemsPerPage": 100}, timeout=120)
            if r.status_code != 200:
                continue
            for rec in r.json().get("data", []):
                sym = rec.get("geneSymbol")
                if sym:
                    cache[sym] = rec.get("gencodeId") or ""
            for s in chunk:
                cache.setdefault(s, "")
            api._sleep()
        cache_put(key, cache)
    return {s: cache.get(s, "") for s in symbols}


def gtex_cis_variants(gencode_id: str) -> List[Dict[str, Any]]:
    """All single-tissue cis-eQTL records for one gene (every GTEx tissue)."""
    key = f"gtex_cis_{gencode_id}"
    hit = cache_get(key)
    if hit is not None:
        return hit
    r = api._get(f"{GTEX_BASE}/association/singleTissueEqtl",
                 params={"gencodeId": gencode_id, "itemsPerPage": 100000}, timeout=240)
    if r.status_code != 200:
        return []
    recs = r.json().get("data", []) or []
    cache_put(key, recs)
    api._sleep()
    return recs


def pick_tissue(recs: Sequence[Dict[str, Any]]) -> Tuple[str, List[Dict[str, Any]]]:
    """Choose the tissue carrying the strongest cis signal for this gene."""
    by_tis: Dict[str, List[Dict[str, Any]]] = {}
    for r in recs:
        by_tis.setdefault(r.get("tissueSiteDetailId") or "NA", []).append(r)
    if not by_tis:
        return "", []
    best, best_key = "", None
    for tis, rows in by_tis.items():
        top_p = min((x.get("pValue") for x in rows if x.get("pValue") is not None),
                    default=1.0)
        z = _z_from_p(top_p) or 0.0
        pref = 1 if tis in GTEX_TISSUE_PREF else 0
        rank = (round(z, 3), len(rows), pref)
        if best_key is None or rank > best_key:
            best_key, best = rank, tis
    return best, by_tis.get(best, [])


def variant_position(rec: Dict[str, Any]) -> Optional[int]:
    vid = rec.get("variantId") or ""
    parts = vid.split("_")
    if len(parts) >= 2:
        try:
            return int(parts[1])
        except ValueError:
            pass
    pos = rec.get("pos")
    return int(pos) if pos is not None else None


def variant_alleles(rec: Dict[str, Any]) -> Tuple[str, str]:
    parts = (rec.get("variantId") or "").split("_")
    if len(parts) >= 4:
        return parts[2].upper(), parts[3].upper()
    return "", ""


def to_effect(rec: Dict[str, Any]) -> Optional[Tuple[float, float, float]]:
    """GTEx gives nes and p; convert to (beta, se, p) with se = |nes|/|z|."""
    nes = _f(rec.get("nes"))
    p = _f(rec.get("pValue"))
    if nes is None or p is None:
        return None
    z = _z_from_p(p)
    if z is None or z <= 0:
        return None
    se = abs(nes) / z
    if se <= 0:
        return None
    return nes, se, p


def thin(items: Sequence[Dict[str, Any]], limit: int) -> List[Dict[str, Any]]:
    """Positional thinning that always keeps the strongest eQTL variant."""
    if len(items) <= limit:
        return list(items)
    ordered = sorted(items, key=lambda r: variant_position(r) or 0)
    step = max(1, len(ordered) // limit)
    keep = ordered[::step][:limit]
    top = max(items, key=lambda r: _z_from_p(_f(r.get("pValue")) or 1.0) or 0.0)
    if top not in keep:
        keep[-1] = top
    return keep


def align(b_e, ea_e, nea_e, b_o, ea_o, nea_o) -> Optional[Tuple[float, float]]:
    ea_e, nea_e = (ea_e or "").upper(), (nea_e or "").upper()
    ea_o, nea_o = (ea_o or "").upper(), (nea_o or "").upper()
    if not ea_e or not ea_o:
        return None
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if ea_e == ea_o and nea_e == nea_o:
        return b_e, b_o
    if ea_e == nea_o and nea_e == ea_o:
        return b_e, -b_o
    if comp.get(ea_e) == ea_o and comp.get(nea_e) == nea_o:
        return b_e, b_o
    if comp.get(ea_e) == nea_o and comp.get(nea_e) == ea_o:
        return b_e, -b_o
    return None


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> None:
    log_init("day3")
    banner("Day 3", "Colocalization, SMR and HEIDI")
    conf = cfg()
    primary = conf["outcomes"]["primary"]
    kb = int(conf["coloc"]["window_kb"])
    cand = candidates()
    log(f"candidate genes carried in from Day 2: {len(cand)}")

    gmap = gtex_gencode([c["symbol"] for c in cand])
    n_map = sum(1 for v in gmap.values() if v)
    log(f"GTEx GENCODE v39 mapping resolved for {n_map}/{len(cand)} genes")

    coloc_rows: List[Dict[str, Any]] = []
    var_rows: List[Dict[str, Any]] = []
    smr_rows: List[Dict[str, Any]] = []
    source_rows: List[Dict[str, Any]] = []
    regional_store: Dict[str, List[Dict[str, Any]]] = {}

    section("3.1 GTEx dense cis-eQTL assembly and coloc.abf")
    for idx, c in enumerate(cand, 1):
        sym = c["symbol"]
        gid = gmap.get(sym) or ""
        if not gid:
            log(f"  [{idx}] {sym}: no GTEx GENCODE record, skipped")
            coloc_rows.append({"symbol": sym, "ensembl_gene_id": c["ensembl_gene_id"],
                               "nsnp": 0, "status": "no_gtex_gene"})
            source_rows.append({"symbol": sym, "ensembl_gene_id": c["ensembl_gene_id"],
                                "status": "no_gtex_gene"})
            continue
        recs = gtex_cis_variants(gid)
        if not recs:
            log(f"  [{idx}] {sym}: no GTEx cis-eQTL in any tissue")
            coloc_rows.append({"symbol": sym, "ensembl_gene_id": c["ensembl_gene_id"],
                               "nsnp": 0, "status": "no_gtex_eqtl"})
            source_rows.append({"symbol": sym, "ensembl_gene_id": c["ensembl_gene_id"],
                                "gtex_gencode_id": gid, "n_tissues_tested": 0,
                                "n_variants_best_tissue": 0, "status": "no_gtex_eqtl"})
            continue

        tis, trows = pick_tissue(recs)
        n_tis = len({r.get("tissueSiteDetailId") for r in recs})
        use = thin(trows, MAX_REGION_VARIANTS)
        top = max(trows, key=lambda r: _z_from_p(_f(r.get("pValue")) or 1.0) or 0.0)
        source_rows.append({
            "symbol": sym, "ensembl_gene_id": c["ensembl_gene_id"],
            "gtex_gencode_id": gid, "n_tissues_tested": n_tis,
            "n_variants_best_tissue": len(trows), "gtex_tissue": tis,
            "top_rsid": top.get("snpId"), "top_nes": _f(top.get("nes")),
            "top_p": _f(top.get("pValue")), "status": "ok",
        })

        rsids = list(dict.fromkeys([r.get("snpId") for r in use if r.get("snpId")]))
        if c["lead_rsid"] and c["lead_rsid"] not in rsids:
            rsids.append(c["lead_rsid"])

        gw = api.og_associations_pair(rsids, primary, batch=ASSOC_BATCH)
        eq_by_rs = {}
        for r in use:
            eff = to_effect(r)
            if eff is None:
                continue
            eq_by_rs[r["snpId"]] = {"beta": eff[0], "se": eff[1], "p": eff[2],
                                    "ea": variant_alleles(r)[0], "nea": variant_alleles(r)[1],
                                    "pos": variant_position(r)}
        common = [s for s in rsids if s in eq_by_rs and s in gw]
        if len(common) < 3:
            log(f"  [{idx}] {sym}: only {len(common)} shared variants, coloc skipped")
            coloc_rows.append({"symbol": sym, "ensembl_gene_id": c["ensembl_gene_id"],
                               "gtex_tissue": tis, "nsnp": len(common),
                               "status": "insufficient_variants"})
            continue

        b1, s1, b2, s2, rows_local = [], [], [], [], []
        for snp in common:
            e, g = eq_by_rs[snp], gw[snp]
            al = align(e["beta"], e["ea"], e["nea"], g["beta"], g["ea"], g["nea"])
            if al is None:
                continue
            b1.append(al[0]); s1.append(e["se"])
            b2.append(al[1]); s2.append(g["se"])
            rows_local.append({"symbol": sym, "rsid": snp, "chr": g.get("chr", ""),
                               "position": e["pos"],
                               "beta_eqtl": al[0], "se_eqtl": e["se"], "p_eqtl": e["p"],
                               "beta_gwas": al[1], "se_gwas": g["se"], "p_gwas": g["p"]})
        if len(b1) < 3:
            log(f"  [{idx}] {sym}: allele harmonisation left <3 variants")
            coloc_rows.append({"symbol": sym, "ensembl_gene_id": c["ensembl_gene_id"],
                               "gtex_tissue": tis, "nsnp": len(b1),
                               "status": "harmonisation_failed"})
            continue

        res = st.coloc_abf(b1, s1, b2, s2)
        regional_store[sym] = rows_local
        for r in rows_local:
            var_rows.append(r)
        lead = c["lead_rsid"]
        lead_row = next((r for r in rows_local if r["rsid"] == lead), None)
        if lead_row is None and rows_local:
            lead_row = max(rows_local, key=lambda r: -(r["p_eqtl"] or 1.0))
        coloc_rows.append({
            "symbol": sym,
            "ensembl_gene_id": c["ensembl_gene_id"],
            "gtex_tissue": tis,
            "eqtl_dataset": c["eqtl_dataset"],
            "chr": rows_local[0]["chr"],
            "gene_start": min(r["position"] or 0 for r in rows_local),
            "gene_end": max(r["position"] or 0 for r in rows_local),
            "window_kb": kb,
            "nsnp": res["nsnp"],
            "PP_H0": res["PP.H0"], "PP_H1": res["PP.H1"], "PP_H2": res["PP.H2"],
            "PP_H3": res["PP.H3"], "PP_H4": res["PP.H4"],
            "coloc_grade": ("strong" if res["PP.H4"] >= conf["coloc"]["pp_h4_strong"]
                            else "suggestive" if res["PP.H4"] >= conf["coloc"]["pp_h4_suggestive"]
                            else "weak"),
            "lead_rsid": lead_row["rsid"] if lead_row else None,
            "lead_p_eqtl": lead_row["p_eqtl"] if lead_row else None,
            "lead_p_gwas": lead_row["p_gwas"] if lead_row else None,
            "status": "ok",
        })
        log(f"  [{idx}/{len(cand)}] {sym} [{tis}]: PP.H4={res['PP.H4']:.3f} "
            f"(n={res['nsnp']}, {len(recs)} GTEx records)")

    write_table(coloc_rows, "day3_coloc_results.csv", COLOC_COLUMNS)
    write_table(var_rows, "day3_coloc_variant_level.csv", VAR_COLUMNS)
    write_table(source_rows, "day3_eqtl_source.csv", SOURCE_COLUMNS)
    strong = sum(1 for r in coloc_rows if r.get("coloc_grade") == "strong")
    sugg = sum(1 for r in coloc_rows if r.get("coloc_grade") == "suggestive")
    log(f"coloc strong (PP.H4>=0.80): {strong}; suggestive (>=0.50): {sugg}")

    section("3.2 SMR and HEIDI")
    for idx, c in enumerate(cand, 1):
        sym = c["symbol"]
        recs = regional_store.get(sym)
        if not recs:
            smr_rows.append({"symbol": sym, "status": "no_variant_level_data"})
            continue
        best = max(recs, key=lambda r: abs(_z(r["beta_eqtl"], r["se_eqtl"]) or 0))
        wald = st.smr_wald(best["beta_eqtl"], best["se_eqtl"],
                           best["beta_gwas"], best["se_gwas"])
        ld_key = f"ld3_{sym}_{best['rsid']}"
        ld = cache_get(ld_key)
        if ld is None:
            ld = {}
            near = sorted(recs, key=lambda r: abs((r["position"] or 0)
                                                  - (best["position"] or 0)))[:LD_PAIRWISE_MAX]
            for r in near:
                if r["rsid"] == best["rsid"]:
                    ld[r["rsid"]] = 1.0
                    continue
                try:
                    pair = api.ensembl_ld_pairwise(best["rsid"], r["rsid"])
                except Exception:
                    pair = None
                if pair and pair.get("r2") is not None:
                    ld[r["rsid"]] = float(pair["r2"])
                time.sleep(0.05)
            cache_put(ld_key, ld)
        be = [r["beta_eqtl"] for r in recs]
        bg = [r["beta_gwas"] for r in recs]
        sg = [r["se_gwas"] for r in recs]
        r2 = [ld.get(r["rsid"]) for r in recs]
        heidi = st.heidi(be, bg, sg, ld_r2=r2)
        smr_rows.append({
            "symbol": sym,
            "ensembl_gene_id": c["ensembl_gene_id"],
            "top_rsid": best["rsid"],
            "b_eqtl": best["beta_eqtl"],
            "se_eqtl": best["se_eqtl"],
            "b_gwas": best["beta_gwas"],
            "se_gwas": best["se_gwas"],
            "b_smr": wald["b_smr"],
            "se_smr": wald["se_smr"],
            "p_smr": wald["p_smr"],
            "heidi_nsnp": heidi["nsnp"],
            "heidi_Q": heidi.get("Q"),
            "heidi_p": heidi.get("heidi_p"),
            "smr_pass": (wald["p_smr"] is not None and wald["p_smr"] < conf["smr"]["p_threshold"]),
            "heidi_pass": (heidi.get("heidi_p") is not None and heidi["heidi_p"] > conf["smr"]["heidi_p_threshold"]),
            "status": "ok",
        })
        if idx % 5 == 0 or idx == len(cand):
            log(f"  [{idx}/{len(cand)}] {sym}: p_SMR={wald['p_smr']}, "
                f"p_HEIDI={heidi.get('heidi_p')} (n_heidi={heidi['nsnp']})")

    write_table(smr_rows, "day3_smr_heidi.csv", SMR_COLUMNS)

    section("3.3 Integrated candidate evidence table")
    coloc_by = {r["symbol"]: r for r in coloc_rows}
    smr_by = {r["symbol"]: r for r in smr_rows}
    evidence = []
    for c in cand:
        sym = c["symbol"]
        co = coloc_by.get(sym, {})
        sm = smr_by.get(sym, {})
        n_ok = sum(1 for flag in [co.get("coloc_grade") in ("strong", "suggestive"),
                                  sm.get("smr_pass"), sm.get("heidi_pass")] if flag)
        evidence.append({
            "symbol": sym,
            "ensembl_gene_id": c["ensembl_gene_id"],
            "druggable_tier": c["druggable_tier"],
            "screen_p": c["screen_p"],
            "ivw_p": c["ivw_p"],
            "ivw_or": c["ivw_or"],
            "nsnp_mr": c["nsnp"],
            "coloc_PP_H4": co.get("PP_H4"),
            "coloc_grade": co.get("coloc_grade", "not_run"),
            "smr_p": sm.get("p_smr"),
            "heidi_p": sm.get("heidi_p"),
            "smr_heidi_pass": bool(sm.get("smr_pass") and sm.get("heidi_pass")),
            "n_supporting_axes": n_ok,
            "prioritised": bool(co.get("coloc_grade") in ("strong", "suggestive")
                                or (sm.get("smr_pass") and sm.get("heidi_pass"))),
        })
    evidence.sort(key=lambda r: (-(r["n_supporting_axes"] or 0), r["screen_p"]))
    write_table(evidence, "day3_candidate_evidence.csv", EVIDENCE_COLUMNS)
    n_pri = sum(1 for e in evidence if e["prioritised"])
    log(f"genes passing coloc or SMR+HEIDI: {n_pri}")

    section("3.4 Figures")
    plot_rows = [r for r in coloc_rows if r.get("PP_H4") is not None]
    plot_rows.sort(key=lambda r: r["PP_H4"], reverse=True)
    plot_rows = plot_rows[:24][::-1]
    if plot_rows:
        fig, ax = plt.subplots(figsize=(6.8, max(3.0, 0.30 * len(plot_rows) + 1.3)))
        y = np.arange(len(plot_rows))
        vals = [r["PP_H4"] for r in plot_rows]
        cols = ["#1f3b57" if v >= 0.8 else "#2a9d8f" if v >= 0.5 else "#c9d3dc" for v in vals]
        ax.barh(y, vals, color=cols, edgecolor="white", height=0.68)
        ax.axvline(0.8, color="#c1443c", ls="--", lw=0.9)
        ax.text(0.805, len(y) - 0.4, "PP.H4 = 0.80", color="#c1443c", fontsize=7)
        ax.set_yticks(y)
        ax.set_yticklabels([r["symbol"] for r in plot_rows], fontsize=7.5)
        ax.set_xlim(0, 1.02)
        ax.set_xlabel("Posterior probability of a shared causal variant (PP.H4)")
        ax.set_title("Fig 8. Colocalization between GTEx cis-eQTL and allergic conjunctivitis GWAS")
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig08_coloc_pp_h4")

    best_sym = None
    best_h4 = -1.0
    for r in coloc_rows:
        h4 = _f(r.get("PP_H4"))
        if h4 is not None and h4 > best_h4:
            best_h4, best_sym = h4, r["symbol"]
    if best_sym and regional_store.get(best_sym):
        recs = sorted(regional_store[best_sym], key=lambda r: r["position"] or 0)
        x = [(r["position"] or 0) / 1e6 for r in recs]
        fig, axes = plt.subplots(2, 1, figsize=(7.4, 5.0), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 1], "hspace": 0.12})
        ax = axes[0]
        ax.scatter(x, [-math.log10(max(r["p_eqtl"], 1e-300)) for r in recs], s=14,
                   c="#1f3b57", linewidths=0)
        ax.set_ylabel("-log10 P (cis-eQTL)")
        ax.set_title(f"Fig 9. Regional colocalization plot, {best_sym} (PP.H4 = {best_h4:.3f})")
        ax.spines[["top", "right"]].set_visible(False)
        ax2 = axes[1]
        ax2.scatter(x, [-math.log10(max(r["p_gwas"], 1e-300)) for r in recs], s=14,
                    c="#2a9d8f", linewidths=0)
        ax2.set_ylabel("-log10 P (GWAS)")
        ax2.set_xlabel("Position on chromosome (Mb, GRCh38)")
        ax2.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig09_coloc_regional_top_gene")

    srows = [r for r in smr_rows if r.get("p_smr") is not None and r.get("heidi_p") is not None]
    if srows:
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
        pass_ = [r for r in srows if r.get("smr_pass") and r.get("heidi_pass")]
        fail_ = [r for r in srows if not (r.get("smr_pass") and r.get("heidi_pass"))]
        ax.scatter([-math.log10(max(r["p_smr"], 1e-300)) for r in fail_],
                   [-math.log10(max(r["heidi_p"], 1e-300)) for r in fail_],
                   s=22, c="#c9d3dc", linewidths=0, label="fails SMR or HEIDI")
        ax.scatter([-math.log10(max(r["p_smr"], 1e-300)) for r in pass_],
                   [-math.log10(max(r["heidi_p"], 1e-300)) for r in pass_],
                   s=30, c="#c1443c", linewidths=0, label="passes SMR + HEIDI")
        ax.axhline(-math.log10(0.05), color="#555555", ls="--", lw=0.8)
        ax.axvline(-math.log10(0.05), color="#555555", ls="--", lw=0.8)
        for r in srows[:14]:
            ax.annotate(r["symbol"], (-math.log10(max(r["p_smr"], 1e-300)),
                                      -math.log10(max(r["heidi_p"], 1e-300))),
                        fontsize=6.5, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel("-log10 P (SMR)")
        ax.set_ylabel("-log10 P (HEIDI heterogeneity)")
        ax.set_title("Fig 10. SMR association versus HEIDI pleiotropy test")
        ax.legend(fontsize=7, loc="lower right")
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig10_smr_heidi_scatter")

    if evidence:
        axes_lbl = ["MR screen\nP<0.05", "IVW P<0.05", "coloc\nPP.H4>=0.5", "SMR+HEIDI\npass"]
        top = evidence[:28]
        mat = []
        for e in top:
            mat.append([
                1.0 if (e["screen_p"] or 1) < 0.05 else 0.0,
                1.0 if (e["ivw_p"] or 1) < 0.05 else 0.0,
                1.0 if (e["coloc_PP_H4"] or 0) >= 0.5 else 0.0,
                1.0 if e["smr_heidi_pass"] else 0.0,
            ])
        fig, ax = plt.subplots(figsize=(4.6, max(3.4, 0.26 * len(top) + 1.2)))
        ax.imshow(np.array(mat), cmap=matplotlib_cmap(), aspect="auto", vmin=0, vmax=1)
        ax.set_xticks(range(len(axes_lbl)))
        ax.set_xticklabels(axes_lbl, fontsize=7)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([e["symbol"] for e in top], fontsize=6.8)
        ax.set_title("Fig 11. Multi-axis evidence matrix for candidate genes")
        for i in range(len(top)):
            for j in range(len(axes_lbl)):
                ax.text(j, i, "Y" if mat[i][j] else "", ha="center", va="center",
                        fontsize=6, color="white")
        save_fig(fig, "fig11_evidence_matrix")

    log("")
    log("Day 3 complete.")


def matplotlib_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("ev", ["#f2f5f8", "#1f3b57"])


def _z(b, se):
    if b is None or se is None or se <= 0:
        return None
    return b / se


if __name__ == "__main__":
    main()
