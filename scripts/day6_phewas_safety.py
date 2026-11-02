"""Day 6 - Phenome-wide safety sweep (PheWAS) as the fallback filter.

Why this day exists
-------------------
A gene that survives MR + coloc + single-cell + CMap is still not a drug target
until its lifelong modulation is shown not to move unrelated disease risk. This
day takes each prioritised gene's lead cis-eQTL and asks what else it does across
the FinnGen R13 phenome (2,700+ endpoints spanning every ICD chapter), then
converts the answer into a three-tier clinical verdict following the scheme used
by the POAG PheWAS paper (PMID 40827287).

Tiers
-----
safe        no phenotype reaches the sweep threshold
monitor     significant only outside the ocular / atopic domain
high_risk   significant in a severe systemic domain (neoplasm, cardiovascular
            death, infection, renal failure, suicide/self-harm)

Inputs
------
tables/day3_candidate_evidence.csv
tables/day2_instruments_lead.csv
data/opengwas_catalog.json, FinnGen R13 phenotype metadata (cached)

Outputs
-------
tables/day6_phenotype_universe.csv
tables/day6_phewas_hits.csv
tables/day6_phewas_all.csv
tables/day6_gene_safety_verdict.csv
figures/svg/fig21_phewas_manhattan.svg
figures/svg/fig22_phewas_domain_burden.svg
figures/svg/fig23_safety_verdict_matrix.svg
"""
from __future__ import annotations

import collections
import json
import math
import time
from typing import Any, Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np

import lib_api as api
from lib_io import (DATA, ROOT, TABLES, banner, cache_get, cache_put, cfg, log,
                    log_init, read_table, save_fig, section, write_table)

PHEWAS_BATCH = 56
MIN_CASES = 50

SEVERE_KEYWORDS = [
    "neoplasm", "cancer", "carcinoma", "malignant", "lymphoma", "leukemia",
    "myocardial", "cardiac arrest", "sudden", "stroke", "intracerebral",
    "pulmonary embolism", "sepsis", "septic", "renal failure", "liver failure",
    "suicide", "self-harm", "intentional", "anaphyla", "aplastic",
    "agranulocytosis", "neutropenia", "tuberculosis", "hepatotoxicity",
]
OCULAR_KEYWORDS = [
    "eye", "ocular", "ophthal", "conjunctiv", "cornea", "corneal", "retina",
    "retinopathy", "glaucoma", "cataract", "uveitis", "iritis", "iridocyclitis",
    "lacrimal", "keratitis", "keratoconus", "macula", "myopia", "blepharitis",
    "chalazion", "pterygium", "vitreous", "optic",
]
ATOPIC_KEYWORDS = [
    "atopic", "allerg", "asthma", "eczema", "dermatitis", "rhinitis",
    "urticaria", "eosinophil", "ige", "hay fever", "conjunctivitis",
]


def phenotype_universe() -> List[Dict[str, Any]]:
    """FinnGen R13 endpoints available in OpenGWAS plus their categories."""
    cached = cache_get("phewas_universe_v2")
    if cached:
        return cached
    catalog = json.loads((DATA / "opengwas_catalog.json").read_text(encoding="utf-8"))
    meta = fetch_finngen_meta()
    rows = []
    for ds, rec in catalog.items():
        if not ds.startswith("finn-b-"):
            continue
        pheno = ds[len("finn-b-"):]
        m = meta.get(pheno, {})
        trait = (rec.get("trait") or "").strip()
        rows.append({
            "dataset": ds,
            "phenocode": pheno,
            "trait": trait,
            "category": m.get("category", "unmapped"),
            "ncase": rec.get("ncase"),
            "ncontrol": rec.get("ncontrol"),
            "sample_size": rec.get("sample_size"),
        })
    cache_put("phewas_universe_v2", rows)
    return rows


def select_phenotypes(universe: List[Dict[str, Any]], max_scan: int) -> List[Dict[str, Any]]:
    """Keep endpoints with enough cases, then cap the sweep deterministically.

    Every endpoint whose label touches the eye, the atopic/immune axis or the
    severe-systemic keyword list is kept unconditionally, because those are the
    ones the safety verdict is actually read from. The remaining budget is filled
    with the largest endpoints, so the sweep stays reproducible across runs.
    """
    keep = []
    for r in universe:
        ncase = r.get("ncase")
        try:
            ncase = int(ncase) if ncase is not None else 0
        except (TypeError, ValueError):
            ncase = 0
        if ncase < MIN_CASES:
            continue
        text = f"{r.get('trait','')} {r.get('category','')}".lower()
        r["_priority"] = bool(
            any(k in text for k in OCULAR_KEYWORDS)
            or any(k in text for k in ATOPIC_KEYWORDS)
            or any(k in text for k in SEVERE_KEYWORDS)
        )
        r["_ncase"] = ncase
        keep.append(r)
    forced = [r for r in keep if r["_priority"]]
    rest = [r for r in keep if not r["_priority"]]
    rest.sort(key=lambda r: -r["_ncase"])
    budget = max(0, max_scan - len(forced))
    chosen = forced + rest[:budget]
    chosen.sort(key=lambda r: (r.get("category") or "", -r["_ncase"]))
    return chosen


def fetch_finngen_meta() -> Dict[str, Dict[str, Any]]:
    cached = cache_get("finngen_phenos_meta")
    if cached:
        return cached
    out: Dict[str, Dict[str, Any]] = {}
    try:
        import requests
        r = requests.get("https://r13.finngen.fi/api/phenos", timeout=180)
        if r.status_code == 200:
            for p in r.json():
                out[p.get("phenocode", "")] = {
                    "category": p.get("category", "unmapped"),
                    "phenostring": p.get("phenostring", ""),
                    "num_cases": p.get("num_cases"),
                    "num_controls": p.get("num_controls"),
                }
    except Exception as exc:
        log(f"  FinnGen metadata fetch failed ({exc}); falling back to trait prefixes")
    cache_put("finngen_phenos_meta", out)
    return out


def domain_of(row: Dict[str, Any]) -> str:
    text = f"{row.get('trait','')} {row.get('category','')}".lower()
    if any(k in text for k in OCULAR_KEYWORDS):
        return "ocular"
    if any(k in text for k in ATOPIC_KEYWORDS):
        return "atopic_immune"
    cat = (row.get("category") or "").lower()
    if cat:
        return cat
    return "other"


def severity_of(row: Dict[str, Any]) -> str:
    text = f"{row.get('trait','')} {row.get('category','')}".lower()
    if any(k in text for k in SEVERE_KEYWORDS):
        return "severe"
    return "non_severe"


def candidate_genes(limit: int = 14) -> List[Dict[str, Any]]:
    ev = read_table("day3_candidate_evidence.csv")
    inst = {r["symbol"]: r for r in read_table("day2_instruments_lead.csv")}
    out = []
    for r in ev:
        sym = r["symbol"]
        if sym not in inst:
            continue
        try:
            p = float(r.get("screen_p") or "1")
        except ValueError:
            p = 1.0
        out.append({
            "symbol": sym,
            "ensembl_gene_id": r.get("ensembl_gene_id", ""),
            "rsid": inst[sym].get("rsid", ""),
            "beta": inst[sym].get("beta", ""),
            "se": inst[sym].get("se", ""),
            "ea": inst[sym].get("ea", ""),
            "nea": inst[sym].get("nea", ""),
            "screen_p": p,
            "coloc_PP_H4": r.get("coloc_PP_H4"),
            "prioritised": r.get("prioritised"),
        })
    pri = [x for x in out if x["prioritised"] == "True"]
    pool = pri or out
    pool.sort(key=lambda x: x["screen_p"])
    return pool[:limit]


def main() -> None:
    log_init("day6")
    banner("Day 6", "Phenome-wide safety sweep")
    conf = cfg()
    thr = float(conf["phewas"]["p"])
    genes = candidate_genes(int(conf["phewas"].get("n_genes", 14)))
    log(f"genes entering PheWAS: {len(genes)} -> {', '.join(g['symbol'] for g in genes)}")

    section("6.1 Phenotype universe")
    universe_all = phenotype_universe()
    log(f"FinnGen R13 endpoints available through OpenGWAS: {len(universe_all)}")
    cat_counts = collections.Counter(r["category"] for r in universe_all)
    log(f"  categories: {len(cat_counts)} (largest: {cat_counts.most_common(3)})")
    write_table(universe_all, "day6_phenotype_universe.csv")
    max_scan = int(conf["phewas"].get("max_phenotypes_scan",
                                     conf["phewas"].get("max_scan", 1400)))
    universe = select_phenotypes(universe_all, max_scan)
    n_forced = sum(1 for r in universe if r.get("_priority"))
    log(f"  swept endpoints: {len(universe)} (>= {MIN_CASES} cases; "
        f"{n_forced} kept unconditionally on eye/atopic/severe labels)")
    write_table([{k: v for k, v in r.items() if not k.startswith("_")} for r in universe],
                "day6_phenotypes_scanned.csv")
    ds_ids = [r["dataset"] for r in universe]

    section("6.2 Per-gene phenome sweep")
    # Resume support. The first full run lost 64 endpoints for one gene because an
    # oversized batch returned a bare 401; without this block the only way to fill the
    # gap would be to re-sweep all 14 genes (~90 min). Rows already on disk are kept,
    # and each gene is queried only for the datasets it is missing.
    all_rows: List[Dict[str, Any]] = []
    if (TABLES / "day6_phewas_all.csv").exists():
        all_rows = read_table("day6_phewas_all.csv")
        # read_table returns everything as str (write_table formats floats with %.6g),
        # so the numeric fields the verdict code compares must be coerced back.
        for r in all_rows:
            r["significant"] = str(r.get("significant")).lower() == "true"
            for k in ("beta", "se", "p", "or", "ncase"):
                v = r.get(k)
                if v in (None, "", "None", "nan"):
                    r[k] = None
                else:
                    try:
                        r[k] = float(v)
                    except (TypeError, ValueError):
                        r[k] = None
        done = collections.defaultdict(set)
        for r in all_rows:
            done[r["gene"]].add(r["dataset"])
        log(f"  resuming: {len(all_rows):,} gene-phenotype rows already on disk")
    else:
        done = collections.defaultdict(set)

    for gi, g in enumerate(genes, 1):
        snp = g["rsid"]
        if not snp:
            log(f"  [{gi}] {g['symbol']}: no lead SNP, skipped")
            continue
        missing = [d for d in ds_ids if d not in done[g["symbol"]]]
        if not missing:
            log(f"  [{gi}/{len(genes)}] {g['symbol']} ({snp}): already complete "
                f"({len(done[g['symbol']])}/{len(ds_ids)}), skipped")
            continue
        t0 = time.time()
        hits = api.og_associations_batched(snp, missing, batch=PHEWAS_BATCH)
        log(f"  [{gi}/{len(genes)}] {g['symbol']} ({snp}): {len(hits)}/{len(missing)} "
            f"missing phenotypes returned, {time.time()-t0:.0f}s "
            f"(total now {len(done[g['symbol']]) + len(hits)}/{len(ds_ids)})")
        for r in universe:
            o = hits.get(r["dataset"])
            if not o:
                continue
            b, se, p = o.get("beta"), o.get("se"), o.get("p")
            if p is None:
                continue
            all_rows.append({
                "gene": g["symbol"],
                "rsid": snp,
                "dataset": r["dataset"],
                "phenocode": r["phenocode"],
                "trait": r["trait"],
                "category": r["category"],
                "domain": domain_of(r),
                "severity": severity_of(r),
                "ncase": r.get("ncase"),
                "beta": b, "se": se, "p": p,
                "or": math.exp(b) if b is not None and abs(b) < 50 else None,
                "ea": o.get("ea"), "nea": o.get("nea"),
                "significant": bool(p is not None and p < thr),
            })
    write_table(all_rows, "day6_phewas_all.csv")
    hits_rows = [r for r in all_rows if r["significant"]]
    hits_rows.sort(key=lambda r: (r["gene"], r["p"]))
    write_table(hits_rows, "day6_phewas_hits.csv")
    log(f"total gene-phenotype tests: {len(all_rows):,}")
    log(f"significant (P<{thr}): {len(hits_rows):,}")

    section("6.3 Three-tier safety verdict")
    verdicts = []
    for g in genes:
        sub = [r for r in all_rows if r["gene"] == g["symbol"]]
        sh = [r for r in sub if r["significant"]]
        dom = collections.Counter(r["domain"] for r in sh)
        sev = [r for r in sh if r["severity"] == "severe"]
        ocular = [r for r in sh if r["domain"] == "ocular"]
        if not sh:
            verdict, reason = "safe", "no phenotype reaches the sweep threshold"
        elif sev:
            verdict = "high_risk"
            reason = ("significant in severe systemic domain: "
                      + "; ".join(sorted({r["trait"] for r in sev})[:4]))
        elif all(r["domain"] in ("ocular", "atopic_immune") for r in sh):
            verdict, reason = "safe", "all significant phenotypes are on-target (ocular/atopic)"
        else:
            verdict = "monitor"
            reason = ("off-domain signals: "
                      + "; ".join(sorted({r["domain"] for r in sh if r["domain"] not in
                                          ("ocular", "atopic_immune")})[:5]))
        verdicts.append({
            "gene": g["symbol"],
            "ensembl_gene_id": g["ensembl_gene_id"],
            "rsid": g["rsid"],
            "screen_p": g["screen_p"],
            "coloc_PP_H4": g["coloc_PP_H4"],
            "n_phenotypes_tested": len(sub),
            "n_significant": len(sh),
            "n_ocular": len(ocular),
            "n_severe": len(sev),
            "top_domains": ";".join(f"{k}:{v}" for k, v in dom.most_common(5)),
            "safety_verdict": verdict,
            "verdict_reason": reason,
        })
    verdicts.sort(key=lambda r: ({"safe": 0, "monitor": 1, "high_risk": 2}[r["safety_verdict"]],
                                 r["screen_p"]))
    write_table(verdicts, "day6_gene_safety_verdict.csv")
    vc = collections.Counter(v["safety_verdict"] for v in verdicts)
    log(f"verdicts: {dict(vc)}")

    section("6.4 Figures")
    # Fig 21 - PheWAS Manhattan for the primary candidate.
    # Sorting verdicts puts a `safe` gene first, and a gene with one ocular hit makes an
    # uninformative Manhattan plot. The figure is meant to show what the lead candidate
    # does across the phenome, so it is pinned to the gene that produced the Day 4
    # signature; if that gene is absent, fall back to the one with the most signals.
    lead = None
    try:
        summ = {r["metric"]: r["value"]
                for r in read_table("day4_signature_summary.csv")}
        cand = summ.get("target_gene")
        if cand and any(v["gene"] == cand for v in verdicts):
            lead = cand
    except FileNotFoundError:
        lead = None
    if lead is None:
        lead = max(verdicts, key=lambda v: v["n_significant"])["gene"] if verdicts else None
    log(f"  Fig 21 plots the phenome sweep of: {lead}")
    if lead:
        sub = [r for r in all_rows if r["gene"] == lead]
        sub.sort(key=lambda r: (r["domain"], r["p"]))
        if sub:
            fig, ax = plt.subplots(figsize=(9.0, 3.8))
            doms = sorted({r["domain"] for r in sub})
            cmap = plt.get_cmap("tab20")
            for i, d in enumerate(doms):
                rows = [r for r in sub if r["domain"] == d]
                ax.scatter([j for j in range(len(rows))],
                           [-math.log10(max(r["p"], 1e-300)) for r in rows],
                           s=10, color=cmap(i % 20), linewidths=0, label=d[:18])
            ax.axhline(-math.log10(thr), color="#c1443c", ls="--", lw=0.8)
            ax.text(0, -math.log10(thr) + 0.15, f"P = {thr:g}", color="#c1443c", fontsize=7)
            ax.set_xlabel("FinnGen R13 endpoints grouped by domain")
            ax.set_ylabel("-log10 P")
            ax.set_title(f"Fig 21. Phenome-wide sweep of the lead cis-eQTL of {lead}")
            ax.set_xticks([])
            ax.legend(fontsize=5.5, ncol=3, loc="upper right")
            ax.spines[["top", "right"]].set_visible(False)
            save_fig(fig, "fig21_phewas_manhattan")

    # Fig 22 - domain burden across genes
    if verdicts:
        doms = collections.Counter()
        for r in hits_rows:
            doms[r["domain"]] += 1
        top = doms.most_common(14)
        if top:
            fig, ax = plt.subplots(figsize=(7.4, 3.6))
            ax.bar([k for k, _ in top], [v for _, v in top], color="#1f3b57",
                   edgecolor="white")
            ax.set_ylabel("Significant gene-phenotype pairs")
            ax.set_xlabel("Phenotype domain")
            ax.tick_params(axis="x", rotation=45, labelsize=7)
            ax.set_title("Fig 22. Distribution of PheWAS signals across phenotype domains")
            ax.spines[["top", "right"]].set_visible(False)
            save_fig(fig, "fig22_phewas_domain_burden")

    # Fig 23 - safety verdict matrix
    if verdicts:
        genes_l = [v["gene"] for v in verdicts]
        metrics = ["n_significant", "n_ocular", "n_severe"]
        mat = np.array([[v[m] for m in metrics] for v in verdicts], dtype=float)
        mat = np.log10(mat + 1)
        fig, ax = plt.subplots(figsize=(4.8, max(3.0, 0.30 * len(genes_l) + 1.3)))
        im = ax.imshow(mat, aspect="auto", cmap="YlOrRd")
        ax.set_xticks(range(len(metrics)))
        ax.set_xticklabels(["all significant", "ocular", "severe"], fontsize=7.5)
        ax.set_yticks(range(len(genes_l)))
        ax.set_yticklabels([f"{v['gene']} [{v['safety_verdict']}]" for v in verdicts], fontsize=7)
        ax.set_title("Fig 23. Safety verdict matrix (log10 count + 1)")
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
        save_fig(fig, "fig23_safety_verdict_matrix")

    log("")
    log("Day 6 complete.")


if __name__ == "__main__":
    main()
