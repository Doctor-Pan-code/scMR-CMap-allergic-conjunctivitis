"""Day 2 - Two-sample cis-eQTL Mendelian randomization screen.

Design (mirrors the framework of Ke et al. TVST 2025, PMID 41165397)
--------------------------------------------------------------------
Stage A  screening pass: one lead cis-eQTL instrument per druggable gene,
         Wald ratio against the allergic-conjunctivitis GWAS. Fast, genome-wide
         over the instrumentable pool.
Stage B  confirmatory pass: for every gene reaching the screening threshold,
         all independent cis-eQTL instruments are used with IVW (FE/RE),
         MR-Egger, weighted median, simple/weighted mode, Cochran Q,
         leave-one-out and the Steiger directionality test.

Inputs
------
tables/day1_druggable_pool.csv
OpenGWAS eqtl-a-ENSG* datasets + FinnGen R13 endpoints

Outputs
-------
tables/day2_screen_lead_snp_mr.csv
tables/day2_instruments_lead.csv
tables/day2_mr_confirmatory.csv
tables/day2_mr_heterogeneity.csv
tables/day2_mr_leave_one_out.csv
tables/day2_mr_steiger.csv
tables/day2_screen_summary.csv
figures/svg/fig04_mr_screen_manhattan.svg
figures/svg/fig05_mr_screen_volcano.svg
figures/svg/fig06_mr_forest_top_genes.svg
figures/svg/fig07_loo_representative.svg
"""
from __future__ import annotations

import collections
import csv
import math
import sys
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

import lib_api as api
import lib_stats as st
from lib_io import (DATA, ROOT, TABLES, banner, cfg, log, log_init, read_table,
                    save_fig, section, write_table)

TOPHITS_BATCH = 500
ASSOC_BATCH = 60


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def load_pool() -> List[Dict[str, str]]:
    rows = read_table("day1_druggable_pool.csv")
    out = []
    for r in rows:
        if r.get("has_cis_eqtl") != "True" and r.get("has_cis_eqtl") != "1":
            continue
        if not r.get("eqtl_dataset"):
            continue
        out.append(r)
    return out


def collect_lead_snps(pool: List[Dict[str, str]]) -> Dict[str, List[Dict[str, Any]]]:
    """Fetch genome-wide-significant cis-eQTL SNPs for every gene dataset."""
    ids = sorted({r["eqtl_dataset"] for r in pool})
    log(f"querying tophits for {len(ids)} eqtl-a datasets (batch={TOPHITS_BATCH})")
    by_gene: Dict[str, List[Dict[str, Any]]] = collections.defaultdict(list)
    t0 = time.time()
    for i in range(0, len(ids), TOPHITS_BATCH):
        chunk = ids[i:i + TOPHITS_BATCH]
        for attempt in range(3):
            try:
                rows = api.og_tophits(chunk, pval=float(cfg()["exposure"]["p_threshold"]))
                break
            except Exception as exc:
                log(f"  tophits batch {i} failed ({exc}); retry {attempt+1}")
                time.sleep(3 * (attempt + 1))
        else:
            rows = []
        for row in rows:
            by_gene[row["id"]].append(row)
        if (i // TOPHITS_BATCH) % 4 == 0:
            log(f"  {i + len(chunk)}/{len(ids)} datasets, {sum(len(v) for v in by_gene.values())} SNP rows, "
                f"{time.time()-t0:.0f}s")
    log(f"tophits done in {time.time()-t0:.0f}s; {len(by_gene)} genes with >=1 cis-eQTL SNP")
    return by_gene


def cis_filter(rows: List[Dict[str, Any]], window_bp: int = 1_000_000) -> List[Dict[str, Any]]:
    """Keep SNPs on the gene's modal chromosome within +/-window of the median position."""
    if not rows:
        return []
    chrom = collections.Counter(r["chr"] for r in rows).most_common(1)[0][0]
    sub = [r for r in rows if r["chr"] == chrom]
    if not sub:
        return []
    pos = sorted(r["position"] for r in sub)
    med = pos[len(pos) // 2]
    keep = [r for r in sub if abs(r["position"] - med) <= window_bp]
    return keep or sub


def positional_prune(rows: List[Dict[str, Any]], kb: int = 500) -> List[Dict[str, Any]]:
    """Approximate LD pruning by position (used when no LD panel is reachable)."""
    rows = sorted(rows, key=lambda r: r["p"])
    kept: List[Dict[str, Any]] = []
    for r in rows:
        if all(abs(r["position"] - k["position"]) > kb * 1000 for k in kept):
            kept.append(r)
    return kept


def align(b_e: float, ea_e: str, nea_e: str, eaf_e: Optional[float],
          b_o: float, ea_o: str, nea_o: str, eaf_o: Optional[float]) -> Optional[Tuple[float, float]]:
    """Harmonise exposure/outcome alleles. Returns (bx, by) or None if unusable."""
    ea_e, nea_e = (ea_e or "").upper(), (nea_e or "").upper()
    ea_o, nea_o = (ea_o or "").upper(), (nea_o or "").upper()
    if not ea_e or not ea_o:
        return None
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    palindromic = {frozenset("AT"), frozenset("CG")}
    if frozenset(ea_e + nea_e) in palindromic:
        # ambiguous strand: keep only when allele frequency is clearly off 0.5
        if eaf_e is None or eaf_o is None:
            return None
        if 0.35 < eaf_e < 0.65 or 0.35 < eaf_o < 0.65:
            return None
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
# Stage A
# ---------------------------------------------------------------------------


def stage_a(pool: List[Dict[str, str]], by_gene: Dict[str, List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
    section("2.1 Stage A - lead cis-eQTL selection and instrument table")
    primary = cfg()["outcomes"]["primary"]
    gene_of = {r["eqtl_dataset"]: r for r in pool}

    lead_rows: List[Dict[str, Any]] = []
    lead_snp_of: Dict[str, Dict[str, Any]] = {}
    for ds, rows in by_gene.items():
        meta = gene_of.get(ds)
        if meta is None:
            continue
        cis = cis_filter(rows)
        if not cis:
            continue
        best = min(cis, key=lambda r: r["p"])
        best = dict(best)
        best["symbol"] = meta["symbol"]
        best["ensembl_gene_id"] = meta["ensembl_gene_id"]
        best["druggable_tier"] = meta["druggable_tier"]
        best["n_cis_snp"] = len(cis)
        lead_snp_of[ds] = best
        lead_rows.append(best)
    log(f"genes with a usable lead cis-eQTL: {len(lead_rows)}")
    write_table(lead_rows, "day2_instruments_lead.csv")

    snps = sorted({r["rsid"] for r in lead_rows if r.get("rsid")})
    log(f"unique lead SNPs to look up in the outcome GWAS: {len(snps)}")

    section("2.2 Stage A - outcome lookup")
    t0 = time.time()
    out_primary = api.og_associations_pair(snps, primary, batch=ASSOC_BATCH)
    log(f"primary outcome rows: {len(out_primary)} / {len(snps)} in {time.time()-t0:.0f}s")

    section("2.3 Stage A - Wald ratio per gene")
    screen: List[Dict[str, Any]] = []
    for ds, inst in lead_snp_of.items():
        o = out_primary.get(inst["rsid"])
        if not o:
            continue
        aligned = align(inst["beta"], inst["ea"], inst["nea"], inst.get("eaf"),
                        o["beta"], o["ea"], o["nea"], o.get("eaf"))
        if aligned is None:
            continue
        bx, by = aligned
        sy = o["se"]
        if not sy or sy <= 0:
            continue
        b = by / bx if abs(bx) > 1e-12 else float("nan")
        se = abs(sy / bx) if abs(bx) > 1e-12 else float("nan")
        if not math.isfinite(b) or not math.isfinite(se) or se <= 0:
            continue
        z = b / se
        p = 2 * st.norm_sf(abs(z))
        f = st.f_statistic([bx], [inst["se"]], int(inst.get("n") or 31470))
        screen.append({
            "symbol": inst["symbol"],
            "ensembl_gene_id": inst["ensembl_gene_id"],
            "druggable_tier": inst["druggable_tier"],
            "eqtl_dataset": ds,
            "rsid": inst["rsid"],
            "chr": inst["chr"],
            "position": inst["position"],
            "ea": inst["ea"], "nea": inst["nea"],
            "eaf_exposure": inst.get("eaf"),
            "beta_exposure": bx,
            "se_exposure": inst["se"],
            "p_exposure": inst["p"],
            "n_exposure": inst.get("n"),
            "beta_outcome": by,
            "se_outcome": sy,
            "p_outcome": o["p"],
            "eaf_outcome": o.get("eaf"),
            "wald_b": b,
            "wald_se": se,
            "wald_or": math.exp(b) if abs(b) < 50 else float("nan"),
            "wald_or_lci": math.exp(b - 1.96 * se) if abs(b) < 50 else float("nan"),
            "wald_or_uci": math.exp(b + 1.96 * se) if abs(b) < 50 else float("nan"),
            "wald_p": p,
            "f_statistic": f,
        })
    log(f"genes with a valid Wald ratio: {len(screen)}")
    ps = [r["wald_p"] for r in screen]
    qs = st.bh_fdr(ps)
    for r, q in zip(screen, qs):
        r["wald_fdr"] = q
    write_table(screen, "day2_screen_lead_snp_mr.csv")

    n_p05 = sum(1 for r in screen if r["wald_p"] < 0.05)
    n_fdr = sum(1 for r in screen if r["wald_fdr"] < float(cfg()["mr"]["fdr_alpha"]))
    n_bonf = sum(1 for r in screen if r["wald_p"] < 0.05 / max(1, len(screen)))
    log(f"  P<0.05: {n_p05}   FDR<0.05: {n_fdr}   Bonferroni: {n_bonf}")
    summary = [
        {"metric": "genes_screened", "value": len(screen)},
        {"metric": "genes_P_lt_0.05", "value": n_p05},
        {"metric": "genes_FDR_lt_0.05", "value": n_fdr},
        {"metric": "genes_bonferroni", "value": n_bonf},
        {"metric": "median_F_statistic", "value": st.__dict__ and _median([r["f_statistic"] for r in screen])},
        {"metric": "outcome", "value": primary},
    ]
    write_table(summary, "day2_screen_summary.csv")
    return screen


# ---------------------------------------------------------------------------
# Stage B
# ---------------------------------------------------------------------------


def stage_b(pool: List[Dict[str, str]], by_gene: Dict[str, List[Dict[str, Any]]],
            screen: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    section("2.4 Stage B - multi-instrument MR for screening hits")
    conf = cfg()
    thr = 0.05
    hits = [r for r in screen if r["wald_p"] < thr]
    hits.sort(key=lambda r: r["wald_p"])
    log(f"genes carried forward (P<0.05): {len(hits)}")
    primary = conf["outcomes"]["primary"]
    replication = conf["outcomes"]["replication"]
    gene_of = {r["eqtl_dataset"]: r for r in pool}

    results: List[Dict[str, Any]] = []
    het_rows: List[Dict[str, Any]] = []
    loo_rows: List[Dict[str, Any]] = []
    steiger_rows: List[Dict[str, Any]] = []

    for idx, hit in enumerate(hits, 1):
        ds = hit["eqtl_dataset"]
        meta = gene_of.get(ds)
        if meta is None:
            continue
        cis = positional_prune(cis_filter(by_gene.get(ds, [])), kb=500)
        if not cis:
            cis = [dict(by_gene[ds][0])]
        snps = [r["rsid"] for r in cis if r.get("rsid")]
        if not snps:
            continue
        try:
            o_rows = api.og_associations_pair(snps, primary, batch=ASSOC_BATCH)
            r_rows = api.og_associations_pair(snps, replication, batch=ASSOC_BATCH)
        except Exception as exc:
            log(f"  [{idx}/{len(hits)}] {hit['symbol']} outcome lookup failed: {exc}")
            continue
        bx, by, sx, sy, used = [], [], [], [], []
        for inst in cis:
            o = o_rows.get(inst["rsid"])
            if not o:
                continue
            al = align(inst["beta"], inst["ea"], inst["nea"], inst.get("eaf"),
                       o["beta"], o["ea"], o["nea"], o.get("eaf"))
            if al is None:
                continue
            bx.append(al[0]); by.append(al[1])
            sx.append(inst["se"]); sy.append(o["se"])
            used.append(inst["rsid"])
        if not bx:
            continue
        n_exp = int(cis[0].get("n") or 31470)
        fstat = st.f_statistic(bx, sx, n_exp)
        row: Dict[str, Any] = {
            "symbol": hit["symbol"],
            "ensembl_gene_id": hit["ensembl_gene_id"],
            "druggable_tier": hit["druggable_tier"],
            "eqtl_dataset": ds,
            "nsnp": len(bx),
            "snps": ";".join(used),
            "f_statistic": fstat,
            "screen_wald_p": hit["wald_p"],
            "screen_wald_or": hit["wald_or"],
        }
        ivw_fe = st.ivw(bx, by, sx, sy, random_effects=False)
        ivw_re = st.ivw(bx, by, sx, sy, random_effects=True)
        egger = st.mr_egger(bx, by, sx, sy)
        wmed = st.weighted_median(bx, by, sx, sy)
        smode = st.simple_mode(bx, by, sy)
        wmode = st.weighted_mode(bx, by, sx, sy)
        for label, res in [("IVW_FE", ivw_fe), ("IVW_RE", ivw_re), ("MR_Egger", egger),
                           ("WeightedMedian", wmed), ("SimpleMode", smode), ("WeightedMode", wmode)]:
            row[f"{label}_b"] = res.get("b")
            row[f"{label}_se"] = res.get("se")
            row[f"{label}_p"] = res.get("p")
            row[f"{label}_or"] = math.exp(res["b"]) if res.get("b") is not None and abs(res["b"]) < 50 else None
            row[f"{label}_or_lci"] = (math.exp(res["b"] - 1.96 * res["se"])
                                      if res.get("b") is not None and res.get("se") and abs(res["b"]) < 50 else None)
            row[f"{label}_or_uci"] = (math.exp(res["b"] + 1.96 * res["se"])
                                      if res.get("b") is not None and res.get("se") and abs(res["b"]) < 50 else None)
        row["egger_intercept"] = egger.get("egger_intercept")
        row["egger_intercept_p"] = egger.get("egger_intercept_p")
        row["Q"] = ivw_re.get("Q")
        row["Q_pval"] = ivw_re.get("Q_pval")
        # replication (IVW on the same instruments)
        rbx, rby, rsx, rsy = [], [], [], []
        for inst in cis:
            o = r_rows.get(inst["rsid"])
            if not o:
                continue
            al = align(inst["beta"], inst["ea"], inst["nea"], inst.get("eaf"),
                       o["beta"], o["ea"], o["nea"], o.get("eaf"))
            if al is None:
                continue
            rbx.append(al[0]); rby.append(al[1]); rsx.append(inst["se"]); rsy.append(o["se"])
        if rbx:
            rep = st.ivw(rbx, rby, rsx, rsy, random_effects=True)
            row["replication_nsnp"] = rep.get("nsnp")
            row["replication_b"] = rep.get("b")
            row["replication_p"] = rep.get("p")
            row["replication_or"] = math.exp(rep["b"]) if rep.get("b") is not None and abs(rep["b"]) < 50 else None
        else:
            row["replication_nsnp"] = 0
        # Steiger
        meta_out = _outcome_meta(primary)
        sg = st.steiger(bx[0], by[0], n_exp, int(meta_out.get("sample_size") or 200000),
                        ncase=meta_out.get("ncase"), ncontrol=meta_out.get("ncontrol"))
        row["steiger_r2_exposure"] = sg["r2_exposure"]
        row["steiger_r2_outcome"] = sg["r2_outcome"]
        row["steiger_correct_direction"] = sg["correct_direction"]
        results.append(row)

        for snp, x, y, s_x, s_y in zip(used, bx, by, sx, sy):
            het_rows.append({"symbol": hit["symbol"], "rsid": snp,
                             "bx": x, "by": y, "sx": s_x, "sy": s_y,
                             "ratio": y / x if abs(x) > 1e-12 else None})
        for lo in st.leave_one_out(bx, by, sx, sy):
            loo_rows.append({"symbol": hit["symbol"],
                             "dropped_rsid": used[lo["dropped"]] if lo["dropped"] < len(used) else "",
                             "ivw_b": lo["b"], "ivw_se": lo["se"], "ivw_p": lo["p"]})
        steiger_rows.append({"symbol": hit["symbol"], **sg})

        if idx % 25 == 0 or idx == len(hits):
            log(f"  [{idx}/{len(hits)}] done, {len(results)} genes with MR results")

    ps = [r["IVW_RE_p"] for r in results if r.get("IVW_RE_p") is not None and math.isfinite(r["IVW_RE_p"])]
    if ps:
        qs = st.bh_fdr(ps)
        it = iter(qs)
        for r in results:
            p = r.get("IVW_RE_p")
            r["IVW_RE_fdr"] = next(it) if p is not None and math.isfinite(p) else None
    write_table(results, "day2_mr_confirmatory.csv")
    write_table(het_rows, "day2_mr_heterogeneity.csv")
    write_table(loo_rows, "day2_mr_leave_one_out.csv")
    write_table(steiger_rows, "day2_mr_steiger.csv")
    return results, het_rows, loo_rows, steiger_rows


_OUTCOME_META: Dict[str, Dict[str, Any]] = {}


def _outcome_meta(dataset: str) -> Dict[str, Any]:
    if dataset in _OUTCOME_META:
        return _OUTCOME_META[dataset]
    info = api.og_gwasinfo([dataset])
    rec = info[0] if info else {}
    _OUTCOME_META[dataset] = rec
    return rec


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------


def figures(screen: List[Dict[str, Any]], results: List[Dict[str, Any]],
            loo_rows: List[Dict[str, Any]]) -> None:
    section("2.5 Figures")

    # Fig 04 - Manhattan style by chromosome
    if screen:
        order = [str(i) for i in range(1, 23)] + ["X"]
        idx = {c: i for i, c in enumerate(order)}
        pts = [(idx.get(r["chr"], 23), r["position"], -math.log10(max(r["wald_p"], 1e-300)), r["symbol"])
               for r in screen]
        pts.sort()
        fig, ax = plt.subplots(figsize=(9.0, 3.8))
        xpos, ypos, cols = [], [], []
        offset = 0
        prev = None
        for c, pos, y, sym in pts:
            if c != prev:
                offset += 1
                prev = c
            xpos.append(offset * 3e7 + pos)
            ypos.append(y)
            cols.append("#1f3b57" if c % 2 == 0 else "#8d99ae")
        ax.scatter(xpos, ypos, s=9, c=cols, linewidths=0)
        ax.axhline(-math.log10(0.05), color="#c1443c", ls="--", lw=0.8)
        ax.text(ax.get_xlim()[0], -math.log10(0.05) + 0.08, "P = 0.05", color="#c1443c", fontsize=7)
        ax.set_xlabel("Chromosome (genes ordered by lead cis-eQTL position)")
        ax.set_ylabel("-log10(P) of cis-eQTL MR (Wald ratio)")
        ax.set_title("Fig 4. Genome-wide cis-eQTL MR screen of the druggable pool in allergic conjunctivitis")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xticks([])
        save_fig(fig, "fig04_mr_screen_manhattan")

    # Fig 05 - volcano
    if screen:
        fig, ax = plt.subplots(figsize=(6.6, 4.6))
        xs = [r["wald_b"] for r in screen]
        ys = [-math.log10(max(r["wald_p"], 1e-300)) for r in screen]
        ax.scatter(xs, ys, s=10, c="#8d99ae", linewidths=0, alpha=0.65, label="P >= 0.05")
        sig = [r for r in screen if r["wald_p"] < 0.05]
        ax.scatter([r["wald_b"] for r in sig], [-math.log10(max(r["wald_p"], 1e-300)) for r in sig],
                   s=18, c="#c1443c", linewidths=0, label="P < 0.05")
        ax.axhline(-math.log10(0.05), color="#c1443c", ls="--", lw=0.8)
        ax.axvline(0, color="#555555", lw=0.6)
        for r in sorted(screen, key=lambda r: r["wald_p"])[:12]:
            ax.annotate(r["symbol"], (r["wald_b"], -math.log10(max(r["wald_p"], 1e-300))),
                        fontsize=6.5, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel("Wald ratio (log OR per 1-SD genetically predicted expression)")
        ax.set_ylabel("-log10(P)")
        ax.set_title("Fig 5. Effect-size versus significance of the cis-eQTL MR screen")
        ax.legend(fontsize=7, loc="upper left")
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig05_mr_screen_volcano")

    # Fig 06 - forest of top genes (IVW)
    ok = [r for r in results if r.get("IVW_RE_or") and r.get("IVW_RE_or_lci") and r.get("IVW_RE_or_uci")]
    ok.sort(key=lambda r: r.get("IVW_RE_p") or 1)
    top = ok[:20][::-1]
    if top:
        fig, ax = plt.subplots(figsize=(6.8, max(3.2, 0.32 * len(top) + 1.4)))
        y = np.arange(len(top))
        ors = [r["IVW_RE_or"] for r in top]
        lo = [r["IVW_RE_or_lci"] for r in top]
        hi = [r["IVW_RE_or_uci"] for r in top]
        ax.errorbar(ors, y, xerr=[np.array(ors) - np.array(lo), np.array(hi) - np.array(ors)],
                    fmt="o", color="#1f3b57", ecolor="#8d99ae", elinewidth=1.2,
                    capsize=2.5, markersize=4)
        ax.axvline(1.0, color="#c1443c", ls="--", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r['symbol']} ({r['nsnp']} SNP)" for r in top], fontsize=7.5)
        ax.set_xscale("log")
        ax.set_xlabel("OR for allergic conjunctivitis per 1-SD higher genetically predicted expression")
        ax.set_title("Fig 6. Multi-instrument IVW estimates for the leading candidate genes")
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig06_mr_forest_top_genes")

    # Fig 07 - leave-one-out for the strongest gene
    if results and loo_rows:
        best = sorted([r for r in results if r.get("IVW_RE_p") is not None],
                      key=lambda r: r["IVW_RE_p"])[0]
        sub = [r for r in loo_rows if r["symbol"] == best["symbol"]]
        if sub:
            fig, ax = plt.subplots(figsize=(6.4, max(2.8, 0.3 * (len(sub) + 1) + 1.2)))
            labels = ["All SNPs"] + [f"drop {r['dropped_rsid']}" for r in sub]
            vals = [best["IVW_RE_b"]] + [r["ivw_b"] for r in sub]
            ses = [best["IVW_RE_se"]] + [r["ivw_se"] for r in sub]
            y = np.arange(len(labels))[::-1]
            ax.errorbar(vals, y, xerr=[1.96 * s for s in ses], fmt="o", color="#1f3b57",
                        ecolor="#8d99ae", capsize=2.5, markersize=4)
            ax.axvline(0, color="#c1443c", ls="--", lw=0.8)
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=7)
            ax.set_xlabel("IVW beta (log OR)")
            ax.set_title(f"Fig 7. Leave-one-out sensitivity analysis, {best['symbol']}")
            ax.spines[["top", "right"]].set_visible(False)
            save_fig(fig, "fig07_loo_representative")


def _median(vals: Sequence[float]) -> float:
    v = sorted(x for x in vals if x is not None and math.isfinite(x))
    if not v:
        return float("nan")
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


# ---------------------------------------------------------------------------


def main() -> None:
    log_init("day2")
    banner("Day 2", "Two-sample cis-eQTL Mendelian randomization screen")
    pool = load_pool()
    log(f"instrumentable druggable genes from Day 1: {len(pool)}")
    by_gene = collect_lead_snps(pool)
    screen = stage_a(pool, by_gene)
    results, het_rows, loo_rows, steiger_rows = stage_b(pool, by_gene, screen)
    figures(screen, results, loo_rows)
    log("")
    log(f"Day 2 complete: {len(screen)} screened, {len(results)} confirmatory MR")


if __name__ == "__main__":
    main()
