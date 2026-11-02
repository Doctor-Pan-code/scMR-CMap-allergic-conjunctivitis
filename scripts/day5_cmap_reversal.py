"""Day 5 - Connectivity Map reversal screening (v2).

Why this day exists
-------------------
Days 2-4 nominate causal genes and describe the transcriptional state that goes
with high expression of the top gene. This day asks the pharmacological question
the two source papers could not answer together: which compounds push the
conjunctival transcriptional state *back* toward the low-expression state?

What changed from v1, and why
-----------------------------
The first run returned 50 + 50 compounds with raw scores of 0.022-0.041, no
compound passed the configured strong threshold of -20, and no glucocorticoid
appeared at all, so the positive control failed outright. Three defects were
found, and fixing them required measuring how the engine behaves rather than
assuming it:

1.  **Query dilution.**  L1000CDS2 matches the query against the ~2,070 genes the
    reference collection covers (the L1000 landmark set plus inferred genes). The
    v7.3 signature has 300 + 300 genes but only 22 + 22 of them are landmark
    genes, so 93% of the query cannot be matched by construction. Measured on a
    known-answer control (a textbook dexamethasone response, which must retrieve
    glucocorticoids), the engine's top score falls monotonically with query size
    (n = 10: 0.53; n = 15: 0.41; n = 25: 0.29; n = 35: 0.23) and its ability to
    recover the known answer degrades with it. The query is therefore run at a
    small size and the size is treated as a design parameter, with the whole
    sweep reported.

2.  **Direction was inferred from the wrong quantity.**  v1 labelled a compound
    "reversing" whenever its mean raw score fell below the empirical median,
    i.e. the *weakest* matches were reported as the best reversers, which is why
    the head of the old table was a block of identical 0.0222 scores. The raw
    score carries no direction at all; direction comes from which arm was
    submitted. The same signature is therefore queried twice -- once as
    (up, dn), which returns compounds that reverse it, and once as (dn, up),
    which returns compounds that reproduce it -- and a compound counts as a
    reversal candidate only when its reverse-query score exceeds its mimic-query
    score.

3.  **The threshold had no scale.**  -20 was inherited from the CMap Web
    convention where connectivity scores run to +/-100; L1000CDS2 raw scores run
    0.02-0.53, so nothing could ever pass. Significance is now a permutation
    p-value: random gene sets of exactly the same size are pushed through the
    identical query to build the null distribution of the top score, and the
    observed score is placed in it. The null is built per query size because the
    null mean itself depends on size (0.195 at n = 10, 0.096 at n = 35).

A known-answer control licenses the whole exercise: the dexamethasone response
must retrieve glucocorticoids, and it does (13 glucocorticoids in the top 50,
best rank 1-6 depending on size). Without that check a null result would be
indistinguishable from a broken pipeline.

Inputs
------
tables/day4_signature_up.csv, tables/day4_signature_dn.csv
tables/day4_signature_summary.csv
data/l1000_landmark_genes.txt

Outputs
-------
tables/day5_l1000cds2_raw.csv
tables/day5_query_size_calibration.csv
tables/day5_null_by_size.csv
tables/day5_cmap_connectivity.csv
tables/day5_positive_control_check.csv
tables/day5_enrichr_lincs_up.csv
tables/day5_enrichr_lincs_dn.csv
tables/day5_enrichr_library_status.csv
tables/day5_cmap_summary.csv
figures/svg/fig17_cmap_score_distribution.svg
figures/svg/fig18_cmap_top_compounds.svg
figures/svg/fig19_cmap_direction_split.svg
figures/svg/fig20_cmap_positive_control.svg
"""
from __future__ import annotations

import random
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

import lib_api as api
from lib_io import (ROOT, banner, cfg, log, log_init, read_table, save_fig,
                    section, write_table)

# ---------------------------------------------------------------------------
# Pharmacological reference sets
# ---------------------------------------------------------------------------
GLUCOCORTICOIDS = {
    "methylprednisolone", "prednisolone", "prednisone", "dexamethasone",
    "fluticasone", "budesonide", "betamethasone", "hydrocortisone",
    "triamcinolone", "triamcinolone acetonide", "mometasone", "ciclesonide",
    "beclomethasone", "desoximetasone", "clobetasol", "clobetasol propionate",
    "fluocinolone", "fluorometholone", "loteprednol", "rimexolone", "medrysone",
    "deflazacort", "cortisone", "cortisol", "halcinonide", "diflorasone",
}
ANTI_ALLERGIC = {
    "montelukast", "zafirlukast", "formoterol", "salmeterol", "albuterol",
    "salbutamol", "cetirizine", "loratadine", "fexofenadine", "desloratadine",
    "olopatadine", "ketotifen", "azelastine", "epinastine", "bepotastine",
    "cromolyn", "nedocromil", "omalizumab", "dupilumab", "mepolizumab",
    "theophylline", "aminophylline", "ipratropium", "tiotropium",
}
IMMUNOSUPPRESSANTS = {
    "cyclosporine", "tacrolimus", "sirolimus", "everolimus", "mycophenolate",
    "azathioprine", "methotrexate", "cyclophosphamide", "leflunomide",
}
ANTI_INFLAMMATORY = {
    "aspirin", "ibuprofen", "indomethacin", "diclofenac", "naproxen",
    "celecoxib", "sulindac", "ketoprofen", "piroxicam", "meloxicam",
    "sulfasalazine", "mesalazine", "balsalazide",
}

# A textbook glucocorticoid response used only as a known-answer control. It is
# a literature signature, not derived from this project's data, so it tests the
# pipeline without circularity.
DEX_UP = ["DUSP1", "KLF15", "TSC22D3", "FKBP5", "PER1", "ANGPTL4", "PDK4",
          "SERPINE1", "MT1X", "MT2A", "CCL2", "IL6", "NFKBIA", "TNFAIP3",
          "SOCS1", "SOCS3", "ZFP36", "IER3", "DDIT4", "TRIB1", "RHOB",
          "MAP3K8", "CEBPD", "NR4A2", "KLF9", "SPRED1", "SPRY2", "ERRFI1",
          "TIPARP", "AXL", "MERTK", "CD163", "FPR1", "SLC2A1"]
DEX_DN = ["IL1B", "IL8", "CXCL1", "CXCL2", "CXCL8", "TNF", "CCL20", "CCL5",
          "ICAM1", "VCAM1", "SELE", "MMP9", "PLA2G4A", "PTGS2", "NOS2",
          "CD80", "CD86", "HLA-DRA", "TAP1", "PSMB9", "IRF1", "STAT1",
          "IFIT1", "IFIT3", "ISG15", "MX1", "OAS1", "GBP1", "GBP2", "IL12B",
          "IL23A", "LTA", "CSF2", "CSF3", "FCGR1A", "TLR2"]


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load_landmarks() -> List[str]:
    path = ROOT / "data" / "l1000_landmark_genes.txt"
    if not path.exists():
        raise RuntimeError(
            "data/l1000_landmark_genes.txt is missing; it lists the genes the "
            "L1000 reference collection measures, and queries outside this set "
            "cannot be matched by the engine")
    return sorted({l.strip().upper()
                   for l in path.read_text(encoding="utf-8").splitlines() if l.strip()})


def load_signature() -> Tuple[List[str], List[str], Dict[str, Any]]:
    up = [r["symbol"] for r in read_table("day4_signature_up.csv") if r.get("symbol")]
    dn = [r["symbol"] for r in read_table("day4_signature_dn.csv") if r.get("symbol")]
    summary = {r["metric"]: r["value"] for r in read_table("day4_signature_summary.csv")}
    return up, dn, summary


# ---------------------------------------------------------------------------
# Engine
# ---------------------------------------------------------------------------
def query(up: Sequence[str], dn: Sequence[str], db: str) -> List[Dict[str, Any]]:
    res = api.cmap_query(up, dn, aggravate=False, db_version=db)
    return (res or {}).get("topMeta") or []


def top_score(meta: Sequence[Dict[str, Any]]) -> Optional[float]:
    s = [m.get("score") for m in meta if m.get("score") is not None]
    return float(max(s)) if s else None


def record(meta: Sequence[Dict[str, Any]], direction: str, label: str,
           size: int) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for rank, m in enumerate(meta, 1):
        ov = m.get("overlap") or {}
        rows.append({
            "query_label": label,
            "query_direction": direction,
            "query_n": size,
            "rank": rank,
            "pert_desc": (m.get("pert_desc") or "").strip(),
            "pert_id": m.get("pert_id"),
            "cell_id": m.get("cell_id"),
            "pert_dose": m.get("pert_dose"),
            "pert_time": m.get("pert_time"),
            "score": m.get("score"),
            "DEGcount": m.get("DEGcount"),
            "pubchem_id": m.get("pubchem_id"),
            "drugbank_id": m.get("drugbank_id"),
            "overlap_up_in_dn": len(ov.get("up/dn") or []),
            "overlap_dn_in_up": len(ov.get("dn/up") or []),
        })
    return rows


def gluco_in(meta: Sequence[Dict[str, Any]]) -> List[Tuple[int, str]]:
    names = [((m.get("pert_desc") or "").strip().lower()) for m in meta]
    return [(i + 1, n) for i, n in enumerate(names) if n in GLUCOCORTICOIDS]


def cls_of(key: str) -> str:
    if key in GLUCOCORTICOIDS:
        return "glucocorticoid"
    if key in ANTI_ALLERGIC:
        return "anti_allergic_or_bronchodilator"
    if key in IMMUNOSUPPRESSANTS:
        return "immunosuppressant"
    if key in ANTI_INFLAMMATORY:
        return "anti_inflammatory"
    return "other"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    log_init("day5")
    banner("Day 5", "Connectivity Map reversal screening")
    conf = cfg()
    db = conf["cmap"]["db_version"]
    cmap_conf = conf["cmap"]
    up_all, dn_all, sig_summary = load_signature()
    target = sig_summary.get("target_gene", "TARGET")
    pool = load_landmarks()
    pool_set = set(pool)

    up_lm = [g for g in up_all if g.upper() in pool_set]
    dn_lm = [g for g in dn_all if g.upper() in pool_set]
    log(f"signature: {len(up_all)} up / {len(dn_all)} down (target gene = {target})")
    log(f"reference-collection overlap: {len(up_lm)} up / {len(dn_lm)} dn "
        f"({100.0*len(up_lm)/max(1,len(up_all)):.1f}% / "
        f"{100.0*len(dn_lm)/max(1,len(dn_all)):.1f}% of the signature is matchable)")

    primary_n = int(cmap_conf.get("primary_query_n", 15))
    sizes = [int(s) for s in cmap_conf.get("size_sweep", [10, 15, 25, 35])]
    if primary_n not in sizes:
        sizes = sorted(set(sizes + [primary_n]))
    null_draws = int(cmap_conf.get("null_draws", 20))
    rng = random.Random(int(cmap_conf.get("null_seed", 20261027)))

    raw_rows: List[Dict[str, Any]] = []

    # -- 5.1 size calibration on the known-answer control -------------------
    section("5.1 Query-size calibration on a known-answer control")
    log("A dexamethasone response signature must retrieve glucocorticoids. "
        "The engine returns compounds that reverse the submitted set, so the "
        "response is submitted with arms swapped; its reverse then retrieves the "
        "compounds that produce it.")
    calib: List[Dict[str, Any]] = []
    best_size, best_rank = None, 10 ** 6
    for n in sizes:
        for orient, a, b in (("swapped", DEX_DN[:n], DEX_UP[:n]),
                             ("normal", DEX_UP[:n], DEX_DN[:n])):
            try:
                meta = query(a, b, db)
            except Exception as exc:
                log(f"  dex n={n} {orient}: failed {type(exc).__name__}")
                continue
            hits = gluco_in(meta)
            rank = min((r for r, _ in hits), default=None)
            calib.append({
                "signature": "dexamethasone_response",
                "query_n": n, "orientation": orient,
                "compounds_returned": len(meta),
                "top_score": round(top_score(meta), 5) if meta else None,
                "n_glucocorticoid": len(hits),
                "best_glucocorticoid_rank": rank,
                "glucocorticoids_retrieved": "; ".join(f"{nm}@{r}" for r, nm in hits[:6]),
                "top5": ", ".join(((m.get("pert_desc") or "").strip() for m in meta[:5])),
            })
            log(f"  dex n={n:<3} {orient:<8} top={top_score(meta):.4f} "
                f"glucocorticoids={len(hits)} best rank={rank}")
            if orient == "swapped" and rank is not None and rank < best_rank:
                best_size, best_rank = n, rank
            time.sleep(0.6)
    log(f"  best known-answer recovery at query size {best_size} "
        f"(best glucocorticoid rank {best_rank})")

    # -- 5.2 observed scores across the same sizes -------------------------
    section("5.2 Signature score across query sizes")
    obs_by_size: Dict[int, Optional[float]] = {}
    for n in sizes:
        try:
            meta = query(up_all[:n], dn_all[:n], db)
        except Exception as exc:
            log(f"  n={n}: failed {type(exc).__name__}")
            continue
        ts = top_score(meta)
        obs_by_size[n] = ts
        hits = gluco_in(meta)
        calib.append({
            "signature": "conjunctival", "query_n": n, "orientation": "reverse",
            "compounds_returned": len(meta),
            "top_score": round(ts, 5) if ts is not None else None,
            "n_glucocorticoid": len(hits),
            "best_glucocorticoid_rank": min((r for r, _ in hits), default=None),
            "glucocorticoids_retrieved": "; ".join(f"{nm}@{r}" for r, nm in hits[:6]),
            "top5": ", ".join(((m.get("pert_desc") or "").strip() for m in meta[:5])),
        })
        log(f"  signature n={n:<3} top={ts:.4f} glucocorticoids={len(hits)}")
        if n == primary_n:
            raw_rows += record(meta, "reverse", "primary", n)
        time.sleep(0.6)
    write_table(calib, "day5_query_size_calibration.csv")

    # -- 5.3 permutation null per size -------------------------------------
    section("5.3 Permutation null per query size")
    log(f"Random gene sets of identical size give the null distribution of the "
        f"top score; {null_draws} draws per size, seed fixed at "
        f"{cmap_conf.get('null_seed', 20261027)}.")
    null_rows: List[Dict[str, Any]] = []
    p_at_size: Dict[int, float] = {}
    for n in sizes:
        nulls: List[float] = []
        null_gluco = 0
        for _ in range(null_draws):
            draw = rng.sample(pool, 2 * n)
            try:
                meta = query(draw[:n], draw[n:], db)
            except Exception:
                continue
            ts = top_score(meta)
            if ts is not None:
                nulls.append(ts)
            if gluco_in(meta):
                null_gluco += 1
            time.sleep(0.5)
        obs = obs_by_size.get(n)
        if not nulls or obs is None:
            continue
        n_ge = sum(1 for t in nulls if t >= obs)
        p = (1 + n_ge) / (1 + len(nulls))
        p_at_size[n] = p
        z = (obs - float(np.mean(nulls))) / max(1e-9, float(np.std(nulls, ddof=1)))
        null_rows.append({
            "query_n": n,
            "observed_top_score": round(obs, 5),
            "null_n": len(nulls),
            "null_mean": round(float(np.mean(nulls)), 5),
            "null_sd": round(float(np.std(nulls, ddof=1)), 5),
            "null_min": round(min(nulls), 5),
            "null_max": round(max(nulls), 5),
            "n_null_ge_observed": n_ge,
            "observed_z": round(z, 2),
            "permutation_p": round(p, 4),
            "null_draws_with_glucocorticoid": null_gluco,
        })
        log(f"  n={n:<3} observed {obs:.4f} vs null {np.mean(nulls):.4f} "
            f"± {np.std(nulls, ddof=1):.4f}  z={z:.2f}  p={p:.3f}")
    write_table(null_rows, "day5_null_by_size.csv")

    p_primary = p_at_size.get(primary_n, 1.0)
    obs_primary = obs_by_size.get(primary_n)
    log(f"primary query size {primary_n}: observed top {obs_primary:.4f}, "
        f"permutation p = {p_primary:.3f}")

    # -- 5.4 mimic query at the primary size -------------------------------
    section("5.4 Mimic query at the primary size")
    mimic_meta: List[Dict[str, Any]] = []
    try:
        mimic_meta = query(dn_all[:primary_n], up_all[:primary_n], db)
        log(f"  mimic: {len(mimic_meta)} compounds, top score "
            f"{top_score(mimic_meta):.4f}")
        raw_rows += record(mimic_meta, "mimic", "primary", primary_n)
    except Exception as exc:
        log(f"  mimic query failed: {type(exc).__name__}: {exc}")

    # -- 5.5 compound-level reversal specificity ---------------------------
    section("5.5 Compound-level reversal specificity")
    rev_scores: Dict[str, List[float]] = {}
    mim_scores: Dict[str, List[float]] = {}
    meta_of: Dict[str, Dict[str, Any]] = {}
    for r in raw_rows:
        name = (r["pert_desc"] or "").strip()
        if not name or r["score"] is None:
            continue
        key = name.lower()
        bucket = rev_scores if r["query_direction"] == "reverse" else mim_scores
        bucket.setdefault(key, []).append(r["score"])
        rec = meta_of.setdefault(key, {"compound": name, "best_rank": 10 ** 6,
                                       "cell": "", "dose": "", "time": "",
                                       "n_instances": 0, "rev_ov": 0, "mim_ov": 0})
        rec["n_instances"] += 1
        ov = r["overlap_up_in_dn"] + r["overlap_dn_in_up"]
        if r["query_direction"] == "reverse":
            rec["rev_ov"] = max(rec["rev_ov"], ov)
            if r["rank"] < rec["best_rank"]:
                rec["best_rank"] = r["rank"]
                rec["cell"] = r["cell_id"]
                rec["dose"] = r["pert_dose"]
                rec["time"] = r["pert_time"]
        else:
            rec["mim_ov"] = max(rec["mim_ov"], ov)

    null_mean_p = next((r["null_mean"] for r in null_rows
                        if r["query_n"] == primary_n), 0.0)
    null_sd_p = next((r["null_sd"] for r in null_rows
                      if r["query_n"] == primary_n), 1.0)
    conn_rows: List[Dict[str, Any]] = []
    for key, rec in meta_of.items():
        rv = rev_scores.get(key) or []
        mm = mim_scores.get(key) or []
        if not rv:
            continue
        rev_best = max(rv)
        mim_best = max(mm) if mm else None
        spec = rev_best - (mim_best if mim_best is not None else null_mean_p)
        z = (rev_best - null_mean_p) / max(1e-9, null_sd_p)
        conn_rows.append({
            "compound": rec["compound"],
            "compound_class": cls_of(key),
            "n_instances": rec["n_instances"],
            "reverse_best_score": round(rev_best, 5),
            "mimic_best_score": round(mim_best, 5) if mim_best is not None else None,
            "reversal_specificity": round(spec, 5),
            "null_z": round(z, 2),
            "reversal_overlap_genes": rec["rev_ov"],
            "mimic_overlap_genes": rec["mim_ov"],
            "best_rank": rec["best_rank"],
            "best_cell_line": rec["cell"],
            "best_dose": rec["dose"],
            "best_time": rec["time"],
            "direction": "reversing" if spec > 0 else "mimicking",
        })
    conn_rows.sort(key=lambda r: (-r["reversal_specificity"], r["best_rank"]))
    write_table(conn_rows, "day5_cmap_connectivity.csv")
    n_rev = sum(1 for r in conn_rows if r["direction"] == "reversing")
    log(f"compounds assessed: {len(conn_rows)}; reversal-specific: {n_rev}")

    strong = [r for r in conn_rows
              if r["direction"] == "reversing" and p_primary <= 0.05]
    log(f"strong candidates (signature-level permutation p = {p_primary:.3f} "
        f"<= 0.05): {len(strong) if p_primary <= 0.05 else 0}")

    # -- 5.6 class control -------------------------------------------------
    section("5.6 Pharmacological class control")
    ctrl_rows: List[Dict[str, Any]] = []
    dex_hits = [(r["query_n"], r["best_glucocorticoid_rank"], r["n_glucocorticoid"])
                for r in calib
                if r["signature"] == "dexamethasone_response"
                and r["orientation"] == "swapped"]
    best_dex = min((h[1] for h in dex_hits if h[1] is not None), default=None)
    ctrl_rows.append({
        "control": "known_answer_dexamethasone",
        "expected": "glucocorticoids occupy the top ranks",
        "n_retrieved": max((h[2] for h in dex_hits), default=0),
        "best_rank": best_dex,
        "pass": bool(best_dex is not None and best_dex <= 10),
        "detail": "; ".join(f"n={n}:best_rank={r},count={c}" for n, r, c in dex_hits),
        "note": ("a textbook glucocorticoid response retrieves glucocorticoids "
                 "from the same engine, so the pipeline can detect a reversal "
                 "signal when one exists"),
    })
    for klass in ("glucocorticoid", "anti_allergic_or_bronchodilator",
                  "immunosuppressant", "anti_inflammatory"):
        sub = [r for r in conn_rows if r["compound_class"] == klass]
        n_r = sum(1 for r in sub if r["direction"] == "reversing")
        ctrl_rows.append({
            "control": f"class_{klass}",
            "expected": "reversing" if klass == "glucocorticoid" else "either",
            "n_retrieved": len(sub),
            "best_rank": min((r["best_rank"] for r in sub), default=None),
            "pass": (n_r > 0) if klass == "glucocorticoid" else None,
            "detail": "; ".join(f"{r['compound']}@{r['best_rank']}" for r in sub[:6]),
            "note": f"{n_r} of {len(sub)} retrieved members are reversal-specific",
        })
    write_table(ctrl_rows, "day5_positive_control_check.csv")
    gluco = [r for r in conn_rows if r["compound_class"] == "glucocorticoid"]
    n_gluco_rev = sum(1 for r in gluco if r["direction"] == "reversing")
    log(f"glucocorticoids in the conjunctival result: {len(gluco)}, "
        f"reversal-specific: {n_gluco_rev}")

    # -- 5.7 Enrichr cross-check -------------------------------------------
    section("5.7 Enrichr perturbation libraries (independent cross-check)")
    # The LINCS_L1000_Chem_Pert_* libraries stream a very large payload and
    # repeatedly abort around 140-170 s from this host, so they are attempted
    # with a short deadline and the nomination falls back to the GEO drug
    # perturbation and LINCS ligand libraries, which return in 1-16 s. A *_down
    # library scores compounds that suppress the submitted set, so the up
    # signature is queried against *_down and the dn signature against *_up.
    LIB_PLAN = [
        ("LINCS_L1000_Chem_Pert_down", "up", 110),
        ("LINCS_L1000_Chem_Pert_up", "dn", 110),
        ("Drug_Perturbations_from_GEO_down", "up", 120),
        ("Drug_Perturbations_from_GEO_up", "dn", 120),
        ("LINCS_L1000_Ligand_Perturbations_down", "up", 90),
        ("LINCS_L1000_Ligand_Perturbations_up", "dn", 90),
    ]
    enrich_rows: Dict[str, List[Dict[str, Any]]] = {"up": [], "dn": []}
    lib_status: List[Dict[str, Any]] = []
    try:
        uid_up = api.enrichr_add_list(up_all[:1500], f"{target}_high_signature")
        uid_dn = api.enrichr_add_list(dn_all[:1500], f"{target}_low_signature")
        log(f"  Enrichr userListIds: up={uid_up}, dn={uid_dn}")
        uid_of = {"up": uid_up, "dn": uid_dn}
        for lib, key, tmo in LIB_PLAN:
            uid = uid_of.get(key)
            if uid is None:
                lib_status.append({"library": lib, "signature": key,
                                   "status": "no_user_list", "n_terms": 0,
                                   "seconds": 0})
                continue
            t0 = time.time()
            try:
                res = api.enrichr_enrich(uid, [lib], timeout=tmo)
                terms = (res or {}).get(lib, [])[:300]
                status = "ok" if terms else "empty"
            except Exception as exc:
                log(f"  {lib}: {type(exc).__name__} after {time.time()-t0:.0f}s, skipped")
                lib_status.append({"library": lib, "signature": key,
                                   "status": f"failed:{type(exc).__name__}",
                                   "n_terms": 0,
                                   "seconds": round(time.time() - t0, 1)})
                continue
            for t in terms:
                enrich_rows[key].append({
                    "library": lib, "rank": t[0], "term": t[1], "p_value": t[2],
                    "z_score": t[3], "combined_score": t[4],
                    "overlap_genes": ";".join(t[5][:20]) if len(t) > 5 else "",
                    "fdr": t[6] if len(t) > 6 else None, "direction": key,
                })
            log(f"  {lib}: {len(terms)} terms in {time.time()-t0:.0f}s")
            lib_status.append({"library": lib, "signature": key, "status": status,
                               "n_terms": len(terms),
                               "seconds": round(time.time() - t0, 1)})
            time.sleep(0.6)
    except Exception as exc:
        log(f"  Enrichr workflow failed: {type(exc).__name__}: {exc}")
    write_table(enrich_rows["up"], "day5_enrichr_lincs_up.csv")
    write_table(enrich_rows["dn"], "day5_enrichr_lincs_dn.csv")
    write_table(lib_status, "day5_enrichr_library_status.csv")

    # -- 5.8 summary -------------------------------------------------------
    section("5.8 Summary")
    summary = [
        {"metric": "signature_up_n", "value": len(up_all)},
        {"metric": "signature_dn_n", "value": len(dn_all)},
        {"metric": "landmark_up_n", "value": len(up_lm)},
        {"metric": "landmark_dn_n", "value": len(dn_lm)},
        {"metric": "primary_query_n", "value": primary_n},
        {"metric": "primary_query_top_score",
         "value": round(obs_primary, 5) if obs_primary is not None else None},
        {"metric": "primary_query_permutation_p", "value": round(p_primary, 4)},
        {"metric": "primary_query_null_mean",
         "value": next((r["null_mean"] for r in null_rows
                        if r["query_n"] == primary_n), None)},
        {"metric": "signature_signal_detected", "value": bool(p_primary <= 0.05)},
        {"metric": "known_answer_best_rank", "value": best_dex},
        {"metric": "known_answer_pass", "value": bool(best_dex is not None and best_dex <= 10)},
        {"metric": "compounds_assessed", "value": len(conn_rows)},
        {"metric": "compounds_reversing", "value": n_rev},
        {"metric": "glucocorticoid_controls", "value": len(gluco)},
        {"metric": "glucocorticoid_reversing", "value": n_gluco_rev},
        {"metric": "enrichr_up_terms", "value": len(enrich_rows["up"])},
        {"metric": "enrichr_dn_terms", "value": len(enrich_rows["dn"])},
    ]
    write_table(summary, "day5_cmap_summary.csv")

    # -- 5.9 figures -------------------------------------------------------
    section("5.9 Figures")
    # Fig 17 - observed against the null at the primary size
    nulls_primary = [r for r in null_rows if r["query_n"] == primary_n]
    if nulls_primary and obs_primary is not None:
        r0 = nulls_primary[0]
        fig, ax = plt.subplots(figsize=(7.0, 3.6))
        ax.bar([0], [r0["null_mean"]], yerr=[r0["null_sd"]], width=0.45,
               color="#c9c9c9", edgecolor="#8d99ae", capsize=4,
               label=f"random gene sets (n={r0['null_n']} draws)")
        ax.bar([1], [obs_primary], width=0.45, color="#1f3b57",
               label="conjunctival signature")
        ax.axhline(r0["null_max"], color="#c1443c", ls=":", lw=0.9,
                   label="null maximum")
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"null\n(mean ± SD, n={primary_n})",
                            f"signature\n(p = {p_primary:.3f})"], fontsize=8)
        ax.set_ylabel("Top L1000CDS2 score")
        ax.set_title(f"Fig 17. Signature score against its permutation null "
                     f"(query size {primary_n})")
        ax.legend(fontsize=7, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig17_cmap_score_distribution")

    # Fig 18 - top compounds by reversal specificity
    cand = conn_rows[:20][::-1]
    if cand:
        fig, ax = plt.subplots(figsize=(7.4, max(3.2, 0.32 * len(cand) + 1.5)))
        y = np.arange(len(cand))
        vals = [r["reversal_specificity"] for r in cand]
        cols = ["#c1443c" if r["compound_class"] != "other" else "#1f3b57"
                for r in cand]
        ax.barh(y, vals, color=cols, edgecolor="white", height=0.7)
        ax.axvline(0, color="#555555", lw=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels([r["compound"][:34] for r in cand], fontsize=7)
        ax.set_xlabel("Reversal specificity (reverse score − mimic score)")
        ax.set_title("Fig 18. Leading reversal-specific compounds")
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig18_cmap_top_compounds")

    # Fig 19 - direction split by class
    if conn_rows:
        classes = ["glucocorticoid", "anti_allergic_or_bronchodilator",
                   "immunosuppressant", "anti_inflammatory", "other"]
        rev_n, mim_n = [], []
        for c in classes:
            sub = [r for r in conn_rows if r["compound_class"] == c]
            rev_n.append(sum(1 for r in sub if r["direction"] == "reversing"))
            mim_n.append(sum(1 for r in sub if r["direction"] == "mimicking"))
        fig, ax = plt.subplots(figsize=(7.2, 3.6))
        x = np.arange(len(classes))
        ax.bar(x - 0.19, rev_n, width=0.38, color="#2a9d8f", label="reversing",
               edgecolor="white")
        ax.bar(x + 0.19, mim_n, width=0.38, color="#c1443c", label="mimicking",
               edgecolor="white")
        ax.set_xticks(x)
        ax.set_xticklabels(["glucocorticoid", "anti-allergic /\nbronchodilator",
                            "immunosuppressant", "anti-\ninflammatory", "other"],
                           fontsize=7.2)
        ax.set_ylabel("Compounds")
        ax.set_title("Fig 19. Direction of connectivity by pharmacological class")
        ax.legend(fontsize=7, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig19_cmap_direction_split")

    # Fig 20 - known-answer control across query sizes
    dex_sw = [r for r in calib if r["signature"] == "dexamethasone_response"
              and r["orientation"] == "swapped"]
    sig_obs = [r for r in calib if r["signature"] == "conjunctival"]
    if dex_sw and sig_obs:
        fig, ax = plt.subplots(figsize=(7.0, 3.6))
        xs = [r["query_n"] for r in dex_sw]
        ax.plot(xs, [r["top_score"] for r in dex_sw], "o-", color="#2a9d8f",
                lw=1.4, ms=4, label="dexamethasone response (known answer)")
        xs2 = [r["query_n"] for r in sig_obs]
        ax.plot(xs2, [r["top_score"] for r in sig_obs], "s--", color="#1f3b57",
                lw=1.4, ms=4, label="conjunctival signature")
        for r in dex_sw:
            if r["best_glucocorticoid_rank"]:
                ax.annotate(f"rank {r['best_glucocorticoid_rank']}",
                            (r["query_n"], r["top_score"]), fontsize=6,
                            xytext=(4, 4), textcoords="offset points",
                            color="#2a9d8f")
        ax.set_xlabel("Query size (genes per arm)")
        ax.set_ylabel("Top L1000CDS2 score")
        ax.set_title("Fig 20. Known-answer control: score and recovery by query size")
        ax.legend(fontsize=7, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig20_cmap_positive_control")

    log("")
    log("Day 5 complete.")


if __name__ == "__main__":
    main()
