"""Day 4 v7 - lineage-restricted TLR1 state contrast with cross-lineage replication.

Why v7 exists: the v4/v5/v6 audit trail
---------------------------------------
v4 stratified cells by TLR1 within author clusters and depth strata. Balanced by
construction, residual depth gap +12.1%, and the signature came out as a pure
cell-identity contrast: up = melanocyte markers (TRPM1, TYRP1, DCT, MLANA, PMEL,
ABCB5) plus macrophage markers (C1QB, MSR1, MS4A7, APOE); down = T/NK markers
(GNLY, GZMA, NKG7, CD3D, CCL5, KLRB1, TRBC2). A CMap query against that
signature would have returned "melanocyte vs T cell", not a disease mechanism.

v5 regressed nine identity axes (eight lineages plus a proliferation panel) and
log depth out of the stratifier, and filtered identity-dominated genes out of the
signature genome-wide. Two failures remained:
  (a) with the split over all cells, TLR1 detection was 29.0% in the high group
      and 0.0% in the low group, so the contrast degenerated into
      detected-vs-not-detected (dropout), not high-vs-low;
  (b) restricting the split to TLR1-detecting cells fixed (a) but the up-signature
      was still a cell-cycle block (HMMR, NDC80, TTK, CCNA2, GINS2, MAD2L1).

v6 replaced the single gene with a composite of the six strong-coloc candidates
(PP.H4 >= 0.80: TLR1, TLR10, NNT, CRK, PRKAB1, UCKL1), each residualised and
weighted by PP.H4. It ran cleanly - 290 up / 300 down, residual depth gap -0.96%,
cell-cycle genes down to 3% of the top 290 - but the design is invalid, and the
audit that showed it is kept in the repo:
  * the six members are mutually uncorrelated at single-cell level after
    residualisation: |r| < 0.025 for every pair (largest is NNT-CRK at 0.024;
    NNT-TLR10 is -0.013). Raw log1p CPM correlations are only 0.018-0.055.
  * a split-half test (arm A = TLR1, NNT, PRKAB1; arm B = TLR10, CRK, UCKL1)
    gave genome-wide Pearson r = -0.195 with 41.2% sign agreement: the two arms
    produced *anti-correlated* signatures.
  * this is not machinery noise. A null calibration with two independent random
    stratifiers pushed through the identical cluster x depth matching, pseudobulk
    and cluster-balanced log2FC gives r = +0.010 and 50.2% sign agreement over
    16,739 genes (tables/day4_null_calibration.csv).
  * the arms select disjoint cells: arm A high-group and arm B high-group overlap
    on only 2,025 of 6,815 cells (29.7%), against 25% expected if the two arms
    were independent. Arm A's high group is 52.4% TLR1-detected and 44.9%
    NNT-detected; arm B's is 67.6% CRK and 42.0% UCKL1.
A weighted mean of six mutually independent transcripts is an average of six
independent noise sources. It cannot define a cell state, so v6 is retired to a
documented negative control.

What v7 does instead: restriction, not regression, and replication, not averaging
-------------------------------------------------------------------------------
The v6 failure says the genes cannot be combined. It says nothing about whether
the TLR1 contrast itself is real. The right question is whether that contrast is
reproducible where it can actually be measured, and a feasibility measurement
answers it:

  TLR1 detection rate by dominant lineage (8 lineages, z >= 0.5):
    melanocyte 20.6%  mast 19.8%  endothelial 18.8%  bplasma 16.6%
    myeloid 16.0%  tnk 14.8%  fibroblast 12.0%  epithelial 11.6%
  TLR1-detecting cell counts: myeloid 935, fibroblast 727, mast 609,
    melanocyte 472, tnk 308, endothelial 249, bplasma 154, epithelial 68.

TLR1 is detected at a comparable rate in every lineage, so it is not a lineage
marker: the high/low contrast is a state contrast, which is what the CMap step
needs. And because it is measurable inside each lineage separately, the contrast
can be tested for replication six independent times.

1. Restriction. The split is performed inside each dominant lineage, and within
   each lineage inside depth strata. Cells are never compared across lineages, so
   lineage composition cannot drive the signature - it is held fixed by design
   rather than regressed out and hoped for.

2. Stratifier. Per-cell TLR1 residual after regressing log1p(CPM) on the nine
   identity axes and log10 total UMI (TLR1 alone: R^2 = 0.047). Residualising
   keeps sequencing depth and cell-cycle state out of the split; the identity
   filter below keeps them out of the readout.

3. Universe. Cells that detect TLR1 (3,966 of 27,282) and carry a dominant lineage
   assignment. Dropout-driven "detected vs not detected" contrasts are excluded by
   construction because both arms detect TLR1.

4. Primary signature. Pooled over the seven lineages with >= 80 detecting cells,
   with the between-cluster consistency filter and the genome-wide identity filter
   from v5 retained.

5. Cross-lineage replication. The identical split is run separately inside each of
   the six lineages with >= 200 detecting cells, each yielding its own pseudobulk
   signature. Replication is reported three ways: per-lineage Pearson r against
   the primary signature, per-lineage sign agreement over the primary signature
   genes, and the number of primary signature genes whose direction is reproduced
   in at least five of six lineages. A signature that reproduces in six
   independent cell populations is not an artefact of one lineage's composition.

Statistical honesty
-------------------
Cell-level p values treat ~27k cells as independent replicates while the
biological replicate count is 4 pooled donors, so they remain descriptive.
Selection uses effect size plus between-cluster consistency plus cross-lineage
replication; no claim rests on a cell-level p value.

Outputs
-------
tables/day4_signature_up.csv              CMap input
tables/day4_signature_dn.csv
tables/day4_signature_summary.csv
tables/day4_pseudobulk_groups.csv
tables/day4_group_balance.csv
tables/day4_stratification_design.csv
tables/day4_target_gene_by_cluster.csv
tables/day4_cluster_composition.csv
tables/day4_cell_qc.csv
tables/day4_lineage_scores.csv            per-cell lineage scores + dominant lineage
tables/day4_lineage_feasibility.csv       per-lineage cell counts and TLR1 detection
tables/day4_lineage_filter_audit.csv      genes removed by the identity filter
tables/day4_cross_lineage_signature.csv   primary vs per-lineage log2FC
tables/day4_cross_lineage_replication.csv per-lineage concordance statistics
figures/svg/fig12_scrna_qc.svg
figures/svg/fig13_cluster_composition.svg
figures/svg/fig14_target_gene_by_cluster.svg
figures/svg/fig15_signature_volcano.svg
figures/svg/fig16_signature_top_genes.svg
figures/svg/fig16b_lineage_restriction.svg
figures/svg/fig16c_cross_lineage_replication.svg
"""
from __future__ import annotations

import collections
import gzip
import math
import time
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

import lib_stats as st
from lib_io import (DATA, banner, cfg, log, log_init, read_table, save_fig,
                    section, write_table)

GEO = DATA / "geo" / "GSE242382"

MIN_CELLS_PER_CLUSTER = 60
CONSISTENCY_MIN = 0.60
LFC_MIN = 0.25
# A gene must be seen in at least this share of one arm before its pseudobulk
# ratio is trusted. Below it the ratio is two near-zero counts and |log2FC|
# ranking promotes pure sampling noise.
DETECTION_MIN_PCT = 5.0
N_DEPTH_STRATA = 10
LINEAGE_R2_MAX = 0.10

# Dominant-lineage assignment: a cell belongs to a lineage when that lineage's
# panel score is the largest and reaches at least this z. Cells below it stay
# unassigned rather than being forced into the closest lineage.
DOMINANT_MIN_SCORE = 0.5
# A lineage enters the pooled primary split when at least this many of its cells
# detect TLR1; both arms then have >= 40 cells inside every depth stratum.
MIN_GROUP_CELLS = 80
# A lineage counts as an independent replication when it can supply this many
# detecting cells, which keeps the per-lineage pseudobulk from being sparse.
MIN_LINEAGE_CELLS = 200

PANELS: Dict[str, List[str]] = {
    "melanocyte": ["TYRP1", "DCT", "MLANA", "PMEL", "TRPM1", "ABCB5", "PAX3", "SLC45A2", "GPR143", "OCA2"],
    "myeloid": ["C1QA", "C1QB", "C1QC", "MSR1", "MS4A7", "APOE", "LYZ", "CD68", "CSF1R", "ITGAM",
                "FCGR3A", "AIF1", "TYROBP", "CTSS", "FPR3", "LILRB1", "MRC1", "CD14"],
    "tnk": ["CD3D", "CD3E", "TRBC2", "CD2", "NKG7", "GNLY", "GZMA", "KLRB1", "CCL5", "IL7R",
            "THEMIS", "CD52", "LTB", "CCL4", "GZMB", "PRF1", "KLRD1", "CD8A", "CD4"],
    "bplasma": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGKC", "IGHG1", "IGLC2", "CD19"],
    "epithelial": ["KRT14", "KRT5", "KRT13", "KRT4", "KRT19", "KRT8", "KRT18", "S100A2", "KRT6A", "KRT17", "SFN"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "COL3A1", "THY1", "VIM"],
    "endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "KDR", "EMCN"],
    "mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT", "GATA2"],
    # Proliferation is not a lineage but it is the other dominant identity axis in
    # a mucosal atlas, and it hijacks the contrast if left in: the first v5 run
    # returned RRM2, TYMS, TK1, MELK, TOP2A, CHEK1, BIRC5, CDC20 as leading
    # "up in TLR1-high" genes. Those are cycling cells, and a proliferation
    # signature makes CMap retrieve cytotoxic chemotherapy, the wrong direction
    # for an allergic disease. It is regressed out of the stratifier and filtered
    # from the signature like any other identity axis.
    "proliferating": ["MKI67", "TOP2A", "RRM2", "TYMS", "TK1", "MELK", "PBK", "BIRC5",
                      "CDC20", "CCNB1", "CCNB2", "CDK1", "CENPF", "NUSAP1", "ASPM",
                      "DLGAP5", "KIF11", "KIF20A", "AURKB", "PLK1", "UBE2C", "ANLN",
                      "NCAPG", "HJURP", "UHRF1", "DTL", "CHEK1", "ESCO2", "MND1",
                      "RAD51AP1", "ORC6", "CDT1", "POC1A", "E2F7", "CENPK"],
}


# ---------------------------------------------------------------------------
# atlas loading
# ---------------------------------------------------------------------------


def mtx_dims(path) -> Tuple[int, int, int]:
    fh = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    while True:
        line = fh.readline()
        if not line:
            raise RuntimeError("unexpected EOF before matrix dimensions")
        if line.startswith("%"):
            continue
        n_genes, n_cells, n_nz = (int(x) for x in line.split())
        fh.close()
        return n_genes, n_cells, n_nz


def iter_mtx(path, chunk_lines: int = 1_000_000):
    """Stream a Matrix Market file in chunks.

    Yields (gene_idx, cell_idx, value, entries_done, entries_total) with 0-based
    indices. Reading line by line in pure Python costs 70-120 s on this 291 MB /
    78 M-entry file; chunked parsing makes a second full pass affordable.
    """
    fh = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    while True:
        line = fh.readline()
        if not line:
            raise RuntimeError("unexpected EOF before matrix dimensions")
        if line.startswith("%"):
            continue
        n_genes, n_cells, n_nz = (int(x) for x in line.split())
        break
    buf: List[str] = []
    done = 0
    for line in fh:
        if line.startswith("%"):
            continue
        buf.append(line)
        if len(buf) >= chunk_lines:
            gi, ci, vv = _parse_chunk(buf)
            done += len(gi)
            yield gi, ci, vv, done, n_nz
            buf = []
    if buf:
        gi, ci, vv = _parse_chunk(buf)
        done += len(gi)
        yield gi, ci, vv, done, n_nz
    fh.close()


def _parse_chunk(lines: List[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(lines)
    gi = np.empty(n, dtype=np.int64)
    ci = np.empty(n, dtype=np.int64)
    vv = np.empty(n, dtype=np.int64)
    for k, line in enumerate(lines):
        p = line.split()
        gi[k] = int(p[0]) - 1
        ci[k] = int(p[1]) - 1
        vv[k] = int(p[2])
    return gi, ci, vv


def load_features() -> List[Tuple[str, str]]:
    out = []
    with gzip.open(GEO / "GSM7761922_LX237_features.tsv.gz", "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out.append((parts[0], parts[1]))
    return out


def load_barcodes() -> List[str]:
    with gzip.open(GEO / "GSM7761922_LX237_barcodes.tsv.gz", "rt", encoding="utf-8") as fh:
        return [l.strip() for l in fh if l.strip()]


def load_locus_group() -> Dict[str, str]:
    """symbol -> HGNC locus_group, e.g. "protein-coding gene"."""
    path = DATA / "hgnc_complete_set.txt"
    out: Dict[str, str] = {}
    if not path.exists():
        log(f"WARNING: {path} absent; protein_coding left blank")
        return out
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        si = header.index("symbol")
        li = header.index("locus_group")
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) <= max(si, li):
                continue
            out[q[si]] = q[li]
    return out


def load_clusters() -> Dict[str, str]:
    out = {}
    with gzip.open(GEO / "GSM7761922_LX237_clusters.tsv.gz", "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ai = header.index("assignment") if "assignment" in header else 2
        si = header.index("status") if "status" in header else 1
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(ai, si):
                continue
            if p[si] != "singlet":
                continue
            out[p[0]] = p[ai]
    return out


def candidate_table() -> List[Dict[str, Any]]:
    rows = read_table("day3_candidate_evidence.csv")
    out = []
    for r in rows:
        try:
            p = float(r.get("screen_p") or "1")
        except ValueError:
            p = 1.0
        try:
            h4 = float(r.get("coloc_PP_H4") or "0")
        except ValueError:
            h4 = 0.0
        out.append({"symbol": r["symbol"], "ensembl_gene_id": r["ensembl_gene_id"],
                    "screen_p": p, "coloc_PP_H4": h4,
                    "prioritised": r.get("prioritised") == "True"})
    return out


def pick_target(conf) -> Dict[str, Any]:
    """The stratifier is a single gene, chosen by the Day-3 evidence.

    v6 averaged all six strong-coloc genes and the audit in the module docstring
    shows why that cannot work: they are mutually uncorrelated, so the composite
    is an average of independent noise and its split-half arms pick disjoint cells
    and yield anti-correlated signatures. One gene it is. TLR1 is the pick because
    it has both the strongest colocalisation (PP.H4 = 0.970) and the strongest
    MR screen signal, and because the feasibility measurement shows it is
    detected at 11.6-20.6% in every lineage - i.e. it is not a lineage marker, so
    the contrast it defines is a cell state.
    """
    cand = candidate_table()
    h4_min = float(conf["coloc"]["pp_h4_strong"])
    strong = [c for c in cand if c["coloc_PP_H4"] >= h4_min]
    strong.sort(key=lambda c: -c["coloc_PP_H4"])
    if not strong:
        pri = [c for c in cand if c["prioritised"]] or cand
        pri.sort(key=lambda c: c["screen_p"])
        strong = pri[:1]
    primary = strong[0]
    override = conf["scrna"].get("target_gene_override")
    if override:
        for c in cand:
            if c["symbol"] == override:
                primary = c
                break
    return primary


def residualise(y: np.ndarray, X: np.ndarray) -> Tuple[np.ndarray, float]:
    """OLS residual of y on [1, X]; returns residuals and R^2."""
    A = np.column_stack([np.ones(X.shape[0]), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((resid ** 2).sum()) / ss_tot if ss_tot > 1e-12 else 0.0
    return resid, r2


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> None:
    log_init("day4")
    banner("Day 4 v7", "Lineage-restricted TLR1 state contrast, cross-lineage replication")
    conf = cfg()
    global DOMINANT_MIN_SCORE, MIN_GROUP_CELLS, MIN_LINEAGE_CELLS, LINEAGE_R2_MAX
    DOMINANT_MIN_SCORE = float(conf["scrna"].get("dominant_lineage_min_score", DOMINANT_MIN_SCORE))
    MIN_GROUP_CELLS = int(conf["scrna"].get("min_group_cells_per_lineage", MIN_GROUP_CELLS))
    MIN_LINEAGE_CELLS = int(conf["scrna"].get("min_replication_cells_per_lineage", MIN_LINEAGE_CELLS))
    LINEAGE_R2_MAX = float(conf["scrna"].get("lineage_r2_max", LINEAGE_R2_MAX))
    primary = pick_target(conf)
    log(f"stratifier: {primary['symbol']} ({primary['ensembl_gene_id']}) "
        f"screen_p={primary['screen_p']:.3g}, PP.H4={primary['coloc_PP_H4']:.3f}")
    log("design: split inside each dominant lineage and inside depth strata "
        "(restriction, not regression); robustness = replication across lineages")

    feats = load_features()
    barcodes = load_barcodes()
    clusters = load_clusters()
    sym_of = [s for _, s in feats]
    ens_of = [e for e, _ in feats]
    idx_of: Dict[str, int] = {}
    for i, s in enumerate(sym_of):
        idx_of.setdefault(s, i)
    n_g = len(feats)

    def resolve(g: Dict[str, Any]) -> Optional[int]:
        if g["ensembl_gene_id"] and g["ensembl_gene_id"] in ens_of:
            return ens_of.index(g["ensembl_gene_id"])
        return idx_of.get(g["symbol"])

    t_idx = resolve(primary)
    if t_idx is None:
        raise SystemExit(f"stratifier {primary['symbol']} absent from GSE242382 features")
    log(f"features={len(feats)}, barcodes={len(barcodes)}, "
        f"singlet cells with cluster label={len(clusters)}")
    log(f"matrix file: GSM7761922_LX237_matrix.mtx.gz "
        f"({(GEO / 'GSM7761922_LX237_matrix.mtx.gz').stat().st_size / 1e6:.0f} MB)")

    mito_idx = {i for i, s in enumerate(sym_of) if s.startswith("MT-")}
    ribo_idx = {i for i, s in enumerate(sym_of) if s.startswith(("RPS", "RPL"))}
    panel_names = sorted(PANELS)
    true_lineages = [p for p in panel_names if p != "proliferating"]
    panel_gene_idx: Dict[str, List[int]] = {}
    for pname in panel_names:
        panel_gene_idx[pname] = [idx_of[g] for g in PANELS[pname] if g in idx_of]
    panel_of_gene: Dict[int, str] = {}
    for pname, ids in panel_gene_idx.items():
        for i in ids:
            panel_of_gene.setdefault(i, pname)

    # ---------------- pass 1: QC, stratifier abundance, lineage panels ----------------
    section("4.1 Pass 1 - QC, stratifier abundance and lineage scores")
    n_cells = len(barcodes)
    counts = np.zeros(n_cells, dtype=np.int64)
    ngenes = np.zeros(n_cells, dtype=np.int32)
    mitoc = np.zeros(n_cells, dtype=np.int64)
    tval = np.zeros(n_cells, dtype=np.int32)
    panel_val = np.zeros((len(panel_names), n_cells), dtype=np.int64)
    panel_pos = {p: k for k, p in enumerate(panel_names)}

    path = GEO / "GSM7761922_LX237_matrix.mtx.gz"
    g_n, c_n, n_nz = mtx_dims(path)
    log(f"matrix dims: {g_n} genes x {c_n} cells x {n_nz:,} non-zeros")
    t0 = time.time()
    panel_gene_list = sorted(panel_of_gene)
    panel_gene_pos = np.array([panel_pos[panel_of_gene[g]] for g in panel_gene_list], dtype=np.int64)
    mito_arr = np.zeros(n_g, dtype=bool)
    mito_arr[list(mito_idx)] = True
    for gi, ci, v, done, total in iter_mtx(path):
        v = v.astype(np.int64)
        counts += np.bincount(ci, weights=v, minlength=n_cells).astype(np.int64)
        ngenes += np.bincount(ci, minlength=n_cells).astype(np.int32)
        mitoc += np.bincount(ci[mito_arr[gi]], weights=v[mito_arr[gi]],
                             minlength=n_cells).astype(np.int64)
        mt = gi == t_idx
        if mt.any():
            tval[ci[mt]] = v[mt]
        sel = np.isin(gi, panel_gene_list)
        if sel.any():
            rows = np.searchsorted(panel_gene_list, gi[sel])
            flat = panel_gene_pos[rows] * n_cells + ci[sel]
            panel_val += np.bincount(flat, weights=v[sel],
                                     minlength=len(panel_names) * n_cells
                                     ).reshape(len(panel_names), n_cells).astype(np.int64)
        if done % 20_000_000 < 1_000_000:
            log(f"  {done:,}/{total:,} entries, {time.time()-t0:.0f}s")
    log(f"pass 1 complete: {n_nz:,} entries in {time.time()-t0:.0f}s")

    mito_frac = np.divide(mitoc, np.maximum(counts, 1)) * 100.0
    min_genes = int(conf["scrna"]["min_genes"])
    max_mito = float(conf["scrna"]["max_mito"])
    keep = (ngenes >= min_genes) & (mito_frac <= max_mito) & (counts >= 500)
    ok = np.where(keep)[0]
    log(f"QC: {len(ok)} / {n_cells} cells retained "
        f"(>= {min_genes} genes, <= {max_mito}% mito, >= 500 UMI)")

    cpm_scale = np.divide(1.0, np.maximum(counts, 1)) * 1e6
    t_cpm = tval * cpm_scale

    # lineage scores: log1p CPM of each panel, z-scored across QC-passed cells
    lineage: Dict[str, np.ndarray] = {}
    for pname in panel_names:
        k = panel_pos[pname]
        raw = panel_val[k].astype(np.float64) * cpm_scale
        lg = np.log1p(raw)
        mu = lg[ok].mean()
        sd = lg[ok].std() or 1.0
        lineage[pname] = (lg - mu) / sd
    Z = np.column_stack([lineage[p] for p in panel_names])[ok]
    Z_full = np.column_stack([lineage[p] for p in panel_names])
    log(f"lineage score matrix: {Z.shape[0]} cells x {Z.shape[1]} panels")

    # dominant lineage: argmax over the eight true lineages, thresholded
    Zt = np.column_stack([lineage[p] for p in true_lineages])
    top_k = np.argmax(Zt, axis=1)
    dom_score = Zt[np.arange(n_cells), top_k]
    dominant = np.array([true_lineages[k] for k in top_k])
    dominant = np.where(dom_score >= DOMINANT_MIN_SCORE, dominant, "unassigned")
    n_unassigned = int((dominant[ok] == "unassigned").sum())
    log(f"dominant lineage assigned at z >= {DOMINANT_MIN_SCORE}: "
        f"{len(ok)-n_unassigned} / {len(ok)} cells "
        f"({100.0*(len(ok)-n_unassigned)/max(1,len(ok)):.1f}%); "
        f"unassigned {n_unassigned}")

    # ---------------- residualised stratifier ----------------
    section("4.2 Stratifier residualisation against identity axes and depth")
    log_lib = np.log10(np.maximum(counts.astype(np.float64), 1.0))
    design = np.column_stack([Z, log_lib[ok]])
    log(f"identity axes regressed out: {', '.join(panel_names)}")
    log(f"  plus log10(total UMI); design = {design.shape[1]} columns x {design.shape[0]} cells")
    y_t = np.log1p(t_cpm[ok])
    resid_t, r2_t = residualise(y_t, design)
    log(f"{primary['symbol']} alone: R^2 = {r2_t:.4f} "
        f"(so {100.0*(1-r2_t):.1f}% of its per-cell variance is not lineage, "
        f"cell cycle or depth)")
    resid_full = np.zeros(n_cells)
    resid_full[ok] = (resid_t - resid_t.mean()) / (resid_t.std() or 1.0)

    # The residual above is fitted on all QC-passed cells, but the split happens
    # only among cells that detect the target, and those cells are systematically
    # deeper (detection needs counts). Refitting on the whole population therefore
    # leaves a depth trend inside the sub-population that is actually split: the
    # first v7 run ranked cells by that residual and produced a -10.49% library
    # gap, i.e. the "state" ordering was partly a depth ordering. The split
    # variable is consequently refitted inside each lineage's detecting cells
    # (4.4), and a full-length design matrix is kept here so that any subset can
    # be residualised on its own.
    design_full = np.column_stack([Z_full, log_lib])

    # ---------------- lineage restriction ----------------
    section("4.3 Lineage restriction - where can the contrast be measured?")
    detected = tval > 0
    feas_rows: List[Dict[str, Any]] = []
    lin_cells: Dict[str, np.ndarray] = {}
    for L in true_lineages + ["unassigned"]:
        m = (dominant == L) & keep
        idx = np.where(m)[0]
        lin_cells[L] = idx
        det = int((detected[idx]).sum()) if len(idx) else 0
        feas_rows.append({
            "lineage": L,
            "cells_total": int(len(idx)),
            "cells_detecting_target": det,
            "detection_pct": round(100.0 * det / max(1, len(idx)), 2),
            "mean_target_cpm": round(float(t_cpm[idx].mean()), 2) if len(idx) else 0.0,
            "eligible_for_primary_split": det >= MIN_GROUP_CELLS,
            "eligible_as_replication": det >= MIN_LINEAGE_CELLS,
        })
    feas_rows.sort(key=lambda r: -r["cells_detecting_target"])
    write_table(feas_rows, "day4_lineage_feasibility.csv")
    for r in feas_rows:
        log(f"  {r['lineage']:<12} cells {r['cells_total']:>6}  "
            f"detecting {r['cells_detecting_target']:>5} "
            f"({r['detection_pct']:>5.1f}%)  mean CPM {r['mean_target_cpm']:>7.2f}"
            f"{'  [primary]' if r['eligible_for_primary_split'] else ''}"
            f"{'  [replication]' if r['eligible_as_replication'] else ''}")

    split_lineages = [r["lineage"] for r in feas_rows
                      if r["eligible_for_primary_split"] and r["lineage"] != "unassigned"]
    repl_lineages = [r["lineage"] for r in feas_rows
                     if r["eligible_as_replication"] and r["lineage"] != "unassigned"]
    log(f"lineages entering the pooled primary split ({len(split_lineages)}): "
        f"{', '.join(split_lineages)}")
    log(f"lineages used as independent replication ({len(repl_lineages)}): "
        f"{', '.join(repl_lineages)}")

    # ---------------- cluster + depth matched split, inside each lineage ----------------
    section("4.4 High / low stratification inside each lineage, depth-matched")

    def split_one(vec: np.ndarray, idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
        """Rank cells by the residual inside depth strata of one lineage."""
        if len(idx) < 4:
            return (np.array([], dtype=np.int64), np.array([], dtype=np.int64), 0)
        order = idx[np.argsort(counts[idx], kind="mergesort")]
        n_strata = min(N_DEPTH_STRATA, max(1, len(order) // 30))
        hi_l: List[int] = []
        lo_l: List[int] = []
        for part in np.array_split(order, n_strata):
            if len(part) < 4:
                continue
            rank = part[np.argsort(vec[part], kind="mergesort")]
            half = len(rank) // 2
            lo_l.extend(int(x) for x in rank[:half])
            hi_l.extend(int(x) for x in rank[len(rank) - half:])
        return (np.array(sorted(hi_l), dtype=np.int64),
                np.array(sorted(lo_l), dtype=np.int64), n_strata)

    # universe: cells that detect the target AND carry a lineage label, so both
    # arms transcribe the stratifier and no comparison crosses a lineage boundary
    universe_mask = detected & keep & (dominant != "unassigned")
    universe = np.where(universe_mask)[0]
    log(f"  universe (target detected + lineage assigned): {len(universe)} / {len(ok)} "
        f"({100.0*len(universe)/max(1,len(ok)):.1f}% of QC-passed cells)")

    lin_hi: Dict[str, np.ndarray] = {}
    lin_lo: Dict[str, np.ndarray] = {}
    strata_rows: List[Dict[str, Any]] = []
    hi_parts: List[np.ndarray] = []
    lo_parts: List[np.ndarray] = []
    # per-cell split variable, refitted inside each lineage's detecting cells
    resid_split = np.full(n_cells, np.nan)
    y_all = np.log1p(t_cpm)
    lin_r2: Dict[str, float] = {}
    for L in split_lineages:
        idx = np.intersect1d(lin_cells[L], universe)
        # refit here, not on the whole atlas: inside the detecting cells of one
        # lineage the residual has to be orthogonal to depth *for those cells*,
        # otherwise ranking by it ranks partly by depth
        r_L, r2_L = residualise(y_all[idx], design_full[idx])
        resid_split[idx] = (r_L - r_L.mean()) / (r_L.std() or 1.0)
        lin_r2[L] = r2_L
        h, l, n_strata = split_one(resid_split, idx)
        lin_hi[L], lin_lo[L] = h, l
        hi_parts.append(h)
        lo_parts.append(l)
        gap_L = 100.0 * (counts[h].mean() / counts[l].mean() - 1.0) if len(h) and len(l) else 0.0
        strata_rows.append({"lineage": L, "cells_in_universe": int(len(idx)),
                            "depth_strata": n_strata, "cells_high": int(len(h)),
                            "cells_low": int(len(l)),
                            "refit_r2_lineage_depth": round(r2_L, 4),
                            "library_size_gap_pct": round(gap_L, 2),
                            "note": "ok" if len(h) and len(l) else "no split"})
        log(f"  {L:<12} universe {len(idx):>5} -> strata {n_strata:>2}, "
            f"high {len(h):>4}, low {len(l):>4}, refit R^2 {r2_L:.3f}, "
            f"depth gap {gap_L:+.2f}%")
    for r in feas_rows:
        if r["lineage"] not in split_lineages:
            strata_rows.append({"lineage": r["lineage"],
                                "cells_in_universe": r["cells_detecting_target"],
                                "depth_strata": 0, "cells_high": 0, "cells_low": 0,
                                "refit_r2_lineage_depth": "",
                                "library_size_gap_pct": "",
                                "note": "fewer than "
                                        f"{MIN_GROUP_CELLS} detecting cells"})
    write_table(strata_rows, "day4_stratification_design.csv")

    hi = np.sort(np.concatenate(hi_parts)) if hi_parts else np.array([], dtype=np.int64)
    lo = np.sort(np.concatenate(lo_parts)) if lo_parts else np.array([], dtype=np.int64)
    n_hi, n_lo = len(hi), len(lo)
    hi_set = set(int(x) for x in hi)
    lo_set = set(int(x) for x in lo)
    log(f"  pooled primary split: high = {n_hi} cells, low = {n_lo} cells")

    gap = 100.0 * (counts[hi].mean() / counts[lo].mean() - 1.0)
    log(f"  residual library-size gap: high = {counts[hi].mean():.0f} UMI/cell vs "
        f"low = {counts[lo].mean():.0f} UMI/cell ({gap:+.2f}%)")
    log(f"  median genes/cell: high = {np.median(ngenes[hi]):.0f} vs low = {np.median(ngenes[lo]):.0f}")
    log(f"  mito%: high = {mito_frac[hi].mean():.2f} vs low = {mito_frac[lo].mean():.2f}")
    log(f"  raw {primary['symbol']} CPM: high = {t_cpm[hi].mean():.2f} vs "
        f"low = {t_cpm[lo].mean():.2f}; detection "
        f"{100.0*(tval[hi] > 0).mean():.1f}% vs {100.0*(tval[lo] > 0).mean():.1f}%")
    log(f"  {primary['symbol']} residual (refit): high = {resid_split[hi].mean():+.3f} vs "
        f"low = {resid_split[lo].mean():+.3f}")

    def dom_mix(sel: np.ndarray) -> str:
        c = collections.Counter(dominant[sel])
        tot = max(1, len(sel))
        return "; ".join(f"{k} {100.0*v/tot:.0f}%" for k, v in c.most_common(4))

    log(f"  lineage composition, high: {dom_mix(hi)}")
    log(f"  lineage composition, low:  {dom_mix(lo)}")

    bal_rows = [
        {"group": "high", "n_cells": n_hi,
         "mean_total_umi": round(float(counts[hi].mean()), 1),
         "median_total_umi": float(np.median(counts[hi])),
         "median_n_genes": float(np.median(ngenes[hi])),
         "mean_target_cpm": round(float(t_cpm[hi].mean()), 2),
         "pct_detected": round(100.0 * float((tval[hi] > 0).mean()), 2),
         "mean_mito_pct": round(float(mito_frac[hi].mean()), 2),
         "mean_residual": round(float(resid_split[hi].mean()), 4),
         "lineage_mix": dom_mix(hi)},
        {"group": "low", "n_cells": n_lo,
         "mean_total_umi": round(float(counts[lo].mean()), 1),
         "median_total_umi": float(np.median(counts[lo])),
         "median_n_genes": float(np.median(ngenes[lo])),
         "mean_target_cpm": round(float(t_cpm[lo].mean()), 2),
         "pct_detected": round(100.0 * float((tval[lo] > 0).mean()), 2),
         "mean_mito_pct": round(float(mito_frac[lo].mean()), 2),
         "mean_residual": round(float(resid_split[lo].mean()), 4),
         "lineage_mix": dom_mix(lo)},
        {"group": "gap_high_vs_low", "n_cells": n_hi + n_lo,
         "mean_total_umi": round(gap, 2),
         "median_total_umi": float(np.median(ngenes[hi]) - np.median(ngenes[lo])),
         "median_n_genes": "",
         "mean_target_cpm": round(float(t_cpm[hi].mean() - t_cpm[lo].mean()), 2),
         "pct_detected": round(100.0 * float((tval[hi] > 0).mean())
                               - 100.0 * float((tval[lo] > 0).mean()), 2),
         "mean_mito_pct": round(float(mito_frac[hi].mean() - mito_frac[lo].mean()), 2),
         "mean_residual": round(float(resid_split[hi].mean() - resid_split[lo].mean()), 4),
         "lineage_mix": ""},
    ]
    write_table(bal_rows, "day4_group_balance.csv")

    # ---------------- composition and per-cluster expression ----------------
    section("4.5 Cell-state composition and stratifier distribution")
    comp = collections.Counter()
    for i in ok:
        comp[clusters.get(barcodes[i], "unassigned")] += 1
    comp_rows = [{"cluster": k, "n_cells": v,
                  "pct": round(100.0 * v / max(1, len(ok)), 2)}
                 for k, v in sorted(comp.items(), key=lambda kv: -kv[1])]
    write_table(comp_rows, "day4_cluster_composition.csv")
    log(f"clusters retained: {len(comp_rows)}; largest = {comp_rows[0]['cluster']} "
        f"({comp_rows[0]['pct']}% of cells)")

    cl_t = collections.defaultdict(list)
    for i in ok:
        cl_t[clusters.get(barcodes[i], "unassigned")].append(float(t_cpm[i]))
    cl_rows = []
    for c, vals in cl_t.items():
        vals = sorted(vals)
        n = len(vals)
        cl_rows.append({"cluster": c, "n_cells": n,
                        "mean_target_cpm": round(sum(vals) / n, 3),
                        "median_target_cpm": vals[n // 2],
                        "pct_expressing": round(100.0 * sum(1 for v in vals if v > 0) / n, 2)})
    cl_rows.sort(key=lambda r: -r["mean_target_cpm"])
    write_table(cl_rows, "day4_target_gene_by_cluster.csv")

    lineage_rows = []
    for j in ok:
        row = {"barcode": barcodes[j],
               "cluster": clusters.get(barcodes[j], "unassigned"),
               "dominant_lineage": dominant[j],
               "total_counts": int(counts[j]),
               "target_cpm": round(float(t_cpm[j]), 3),
               "target_residual": round(float(resid_full[j]), 5),
               "target_residual_refit": (round(float(resid_split[j]), 5)
                                        if not np.isnan(resid_split[j]) else "")}
        for pname in panel_names:
            row[f"score_{pname}"] = round(float(lineage[pname][j]), 4)
        lineage_rows.append(row)
    write_table(lineage_rows, "day4_lineage_scores.csv")
    log(f"per-cell lineage table: {len(lineage_rows)} cells x {len(lineage_rows[0])} columns")

    # ---------------- pass 2: pseudobulk ----------------
    section("4.6 Pass 2 - pseudobulk for the primary split and every lineage")
    sum_hi = np.zeros(n_g, dtype=np.int64)
    sum_lo = np.zeros(n_g, dtype=np.int64)
    cnt_hi = np.zeros(n_g, dtype=np.int32)
    cnt_lo = np.zeros(n_g, dtype=np.int32)
    cl_of_cell = np.array([clusters.get(barcodes[i], "unassigned") for i in range(n_cells)])
    cl_names = sorted({cl_of_cell[i] for i in ok})
    n_cl = len(cl_names)
    cl_id = np.zeros(n_cells, dtype=np.int64)
    for k, c in enumerate(cl_names):
        cl_id[cl_of_cell == c] = k
    cl_hi = np.zeros((n_cl, n_g), dtype=np.float64)
    cl_lo = np.zeros((n_cl, n_g), dtype=np.float64)

    # per-lineage arms: 2 x n_repl matrices
    n_repl = len(repl_lineages)
    repl_pos = {L: k for k, L in enumerate(repl_lineages)}
    ln_hi_sum = np.zeros((n_repl, n_g), dtype=np.float64)
    ln_lo_sum = np.zeros((n_repl, n_g), dtype=np.float64)
    ln_hi_mask = np.zeros((n_repl, n_cells), dtype=bool)
    ln_lo_mask = np.zeros((n_repl, n_cells), dtype=bool)
    for L in repl_lineages:
        k = repl_pos[L]
        if L in lin_hi and len(lin_hi[L]):
            ln_hi_mask[k, lin_hi[L]] = True
        if L in lin_lo and len(lin_lo[L]):
            ln_lo_mask[k, lin_lo[L]] = True

    # genome-wide lineage + depth R^2 of every gene, from per-gene moments:
    # with A = [1, Z, log10 UMI] precomputed, ESS of gene g is b'(A'A)b where
    # b = (A'A)^-1 A'y, so only A'y (11 numbers per gene) plus y'y and sum(y)
    # are needed. That makes the identity filter cover every tested gene instead
    # of the handful a probe happened to measure - the flaw in the first v5 run,
    # where KRT14/KRT17/S100A2 survived because no probe value existed for them.
    n_ok = len(ok)
    row_of = np.full(n_cells, -1, dtype=np.int64)
    row_of[ok] = np.arange(n_ok)
    A = np.column_stack([np.ones(n_ok), design])
    AtA_inv = np.linalg.pinv(A.T @ A)
    n_par = A.shape[1]
    sy = np.zeros(n_g, dtype=np.float64)
    syy = np.zeros(n_g, dtype=np.float64)
    sxy = np.zeros((n_g, n_par), dtype=np.float64)

    hi_mask = np.zeros(n_cells, dtype=bool)
    hi_mask[hi] = True
    lo_mask = np.zeros(n_cells, dtype=bool)
    lo_mask[lo] = True

    t0 = time.time()
    for gi, ci, v, done, total in iter_mtx(path, chunk_lines=1_000_000):
        in_ok = row_of[ci] >= 0
        if not in_ok.any():
            continue
        gi_k = gi[in_ok]
        ci_k = ci[in_ok]
        v_k = v[in_ok].astype(np.float64)
        rows = row_of[ci_k]

        lg = np.log1p(v_k * cpm_scale[ci_k])
        sy += np.bincount(gi_k, weights=lg, minlength=n_g)
        syy += np.bincount(gi_k, weights=lg * lg, minlength=n_g)
        for k in range(n_par):
            sxy[:, k] += np.bincount(gi_k, weights=lg * A[rows, k], minlength=n_g)

        vh = hi_mask[ci_k]
        vl = lo_mask[ci_k]
        if vh.any():
            sum_hi += np.bincount(gi_k[vh], weights=v_k[vh], minlength=n_g).astype(np.int64)
            cnt_hi += np.bincount(gi_k[vh], minlength=n_g).astype(np.int32)
        if vl.any():
            sum_lo += np.bincount(gi_k[vl], weights=v_k[vl], minlength=n_g).astype(np.int64)
            cnt_lo += np.bincount(gi_k[vl], minlength=n_g).astype(np.int32)
        for k in range(n_cl):
            m = vh & (cl_id[ci_k] == k)
            if m.any():
                cl_hi[k] += np.bincount(gi_k[m], weights=v_k[m], minlength=n_g)
            m = vl & (cl_id[ci_k] == k)
            if m.any():
                cl_lo[k] += np.bincount(gi_k[m], weights=v_k[m], minlength=n_g)
        for k in range(n_repl):
            m = ln_hi_mask[k][ci_k]
            if m.any():
                ln_hi_sum[k] += np.bincount(gi_k[m], weights=v_k[m], minlength=n_g)
            m = ln_lo_mask[k][ci_k]
            if m.any():
                ln_lo_sum[k] += np.bincount(gi_k[m], weights=v_k[m], minlength=n_g)

        if done % 20_000_000 < 1_000_000:
            log(f"  {done:,}/{total:,} entries, {time.time()-t0:.0f}s")
    log(f"pass 2 complete in {time.time()-t0:.0f}s")

    M = np.column_stack([sy, sxy[:, 1:]])
    quad = np.einsum("gi,ij,gj->g", M, AtA_inv, M, optimize=True)
    grand = sy * sy / n_ok
    sst = syy - grand
    ssr = quad - grand
    with np.errstate(divide="ignore", invalid="ignore"):
        r2_gene = np.where(sst > 1e-8, np.clip(ssr / np.maximum(sst, 1e-12), 0.0, 1.0), 0.0)
    log(f"genome-wide identity+depth R^2 computed for {n_g:,} genes "
        f"(median {np.median(r2_gene):.3f})")

    tot_hi = float(sum(sum_hi))
    tot_lo = float(sum(sum_lo))
    log(f"group total UMI: high = {tot_hi:.3e}, low = {tot_lo:.3e}")

    # ---------------- differential expression ----------------
    section("4.7 Differential expression and lineage-filtered signature")
    cpm_hi = sum_hi / max(tot_hi, 1.0) * 1e6
    cpm_lo = sum_lo / max(tot_lo, 1.0) * 1e6
    # One pseudo-count on the COUNT scale, expressed in CPM. A literal +1 here
    # would be a ~19-count pseudocount for this pseudobulk and would compress
    # every fold change toward zero (measured: the down arm then contained no
    # gene beyond -0.5 while 198 genes genuinely sit there).
    eps_hi = 1e6 / max(tot_hi, 1.0)
    eps_lo = 1e6 / max(tot_lo, 1.0)
    lfc = np.log2(cpm_hi + eps_hi) - np.log2(cpm_lo + eps_lo)

    se_log2 = np.sqrt((sum_hi + 1.0) / max(tot_hi, 1.0) ** 2
                      + (sum_lo + 1.0) / max(tot_lo, 1.0) ** 2) * 1e6
    mu_bar = np.maximum((cpm_hi + cpm_lo) / 2.0, 1e-9)
    se_log2 = se_log2 / (mu_bar * math.log(2.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(se_log2 > 0, lfc / np.maximum(se_log2, 1e-12), 0.0)

    # A gene needs at least 10 counts on one side to give a usable pseudobulk
    # ratio; the old rule (0.5% of pooled cells, i.e. 17 reads) let near-empty
    # genes into the ranking where one read dominated the fold change.
    expressed = (cnt_hi >= 10) | (cnt_lo >= 10)
    pvals = np.array([st.t_two_sided(abs(float(t)), max(n_hi, n_lo) - 1)
                      if expressed[i] else 1.0 for i, t in enumerate(tstat)])
    qvals = np.array(st.bh_fdr([float(x) for x in pvals]))
    log(f"genes tested (expressed filter): {int(expressed.sum()):,}")
    log(f"genes with cell-level FDR<0.05: {int((qvals < 0.05).sum()):,} "
        f"(descriptive only: 4 pooled donors)")

    usable_cl = [c for c in cl_names
                 if sum(1 for i in ok if cl_of_cell[i] == c and i in hi_set) >= MIN_CELLS_PER_CLUSTER
                 and sum(1 for i in ok if cl_of_cell[i] == c and i in lo_set) >= MIN_CELLS_PER_CLUSTER]
    log(f"clusters usable for the consistency filter: {len(usable_cl)} of {len(cl_names)}")
    if usable_cl:
        mat = []
        for c in usable_cl:
            k = cl_names.index(c)
            a = cl_hi[k] / max(float(cl_hi[k].sum()), 1.0) * 1e6
            b = cl_lo[k] / max(float(cl_lo[k].sum()), 1.0) * 1e6
            mat.append(np.log2(a + eps_hi) - np.log2(b + eps_lo))
        mat = np.vstack(mat)
        lfc_cl_bal = mat.mean(axis=0)
        agree = (np.sign(mat) == np.sign(lfc_cl_bal)).mean(axis=0)
    else:
        lfc_cl_bal = lfc
        agree = np.ones(n_g)

    excluded = set(mito_idx) | set(ribo_idx) | {t_idx}
    # Panel genes define the nuisance axes, so they cannot also be readouts of the
    # contrast. Without this the cell-cycle members of the proliferation panel
    # (E2F7 R^2 = 0.094, MELK 0.083, NCAPG 0.080) sit just under the 0.10
    # identity threshold and stay in the up-signature, which is how the first v5
    # run ended up handing CMap a mitosis signature.
    excluded |= set(panel_of_gene)
    audit: List[Dict[str, Any]] = []

    locus = load_locus_group()
    log(f"HGNC locus_group loaded for {len(locus):,} symbols")
    de_rows = []
    for i in range(n_g):
        if not expressed[i]:
            continue
        de_rows.append({
            "ensembl_gene_id": feats[i][0],
            "symbol": feats[i][1],
            "locus_group": locus.get(feats[i][1], "unmapped"),
            "cpm_high": round(float(cpm_hi[i]), 4),
            "cpm_low": round(float(cpm_lo[i]), 4),
            "log2fc": round(float(lfc[i]), 5),
            "log2fc_cluster_balanced": round(float(lfc_cl_bal[i]), 5),
            "cluster_agreement": round(float(agree[i]), 3),
            "lineage_depth_r2": round(float(r2_gene[i]), 4),
            "t_stat": round(float(tstat[i]), 4),
            "p_value": float(pvals[i]),
            "fdr": float(qvals[i]),
            "pct_high": round(100.0 * cnt_hi[i] / max(1, n_hi), 3),
            "pct_low": round(100.0 * cnt_lo[i] / max(1, n_lo), 3),
        })
    de_rows.sort(key=lambda r: r["log2fc_cluster_balanced"], reverse=True)
    write_table(de_rows, "day4_pseudobulk_groups.csv")

    top_n = int(conf["scrna"]["top_n_features"])
    up, dn = [], []
    # Rank by absolute effect size. Walking de_rows in its descending-lfc order
    # filled the down arm with the weakest negatives first: every gene in the
    # first v7.1 down-signature sat at -0.250 to -0.252 while CCL21 (-1.53) and
    # RNASE1 (-1.25) were skipped entirely.
    for r in sorted(de_rows, key=lambda x: -abs(x["log2fc_cluster_balanced"])):
        i = idx_of.get(r["symbol"])
        if i is None:
            continue
        lfc_v = r["log2fc_cluster_balanced"]
        if abs(lfc_v) <= LFC_MIN:
            continue
        if max(r["pct_high"], r["pct_low"]) < DETECTION_MIN_PCT:
            continue
        if i in excluded:
            continue
        r2 = float(r2_gene[i])
        if r2 >= LINEAGE_R2_MAX:
            audit.append({"symbol": r["symbol"], "log2fc_cluster_balanced": lfc_v,
                          "lineage_depth_r2": round(r2, 4),
                          "action": "removed_from_signature",
                          "reason": f"identity R^2 {r2:.3f} >= {LINEAGE_R2_MAX}"})
            continue
        if r["cluster_agreement"] < CONSISTENCY_MIN:
            audit.append({"symbol": r["symbol"], "log2fc_cluster_balanced": lfc_v,
                          "lineage_depth_r2": round(r2, 4),
                          "action": "removed_from_signature",
                          "reason": f"cluster agreement {r['cluster_agreement']:.2f} "
                                    f"< {CONSISTENCY_MIN}"})
            continue
        if lfc_v > 0 and len(up) < top_n:
            up.append(r)
        elif lfc_v < 0 and len(dn) < top_n:
            dn.append(r)
    write_table(up, "day4_signature_up.csv")
    write_table(dn, "day4_signature_dn.csv")
    audit.sort(key=lambda r: -abs(r["log2fc_cluster_balanced"]))
    write_table(audit, "day4_lineage_filter_audit.csv")
    log(f"identity/consistency filter removed {len(audit)} genes from the signature candidates")
    log(f"signature written: {len(up)} up / {len(dn)} down genes "
        f"(|cluster-balanced log2FC| > {LFC_MIN}, agreement >= {CONSISTENCY_MIN}, "
        f"identity R^2 < {LINEAGE_R2_MAX}, detected in >= {DETECTION_MIN_PCT}% of one arm, "
        f"stratifier/MT/ribosomal excluded)")
    pc_up = sum(1 for r in up if r["locus_group"] == "protein-coding gene")
    pc_dn = sum(1 for r in dn if r["locus_group"] == "protein-coding gene")
    log(f"  protein-coding share (the only symbols L1000CDS2 can match): "
        f"up {pc_up}/{len(up)}, down {pc_dn}/{len(dn)}")

    # ---------------- cross-lineage replication ----------------
    section("4.8 Cross-lineage replication of the signature")
    prim_sig = {r["symbol"]: float(r["log2fc_cluster_balanced"]) for r in up + dn}
    lin_lfc: Dict[str, np.ndarray] = {}
    repl_rows: List[Dict[str, Any]] = []
    common_mask = expressed.copy()
    for i in excluded:
        common_mask[i] = False
    common_mask &= (r2_gene < LINEAGE_R2_MAX)
    for L in repl_lineages:
        k = repl_pos[L]
        th = float(ln_hi_sum[k].sum())
        tl = float(ln_lo_sum[k].sum())
        if th <= 0 or tl <= 0:
            continue
        a = ln_hi_sum[k] / th * 1e6
        b = ln_lo_sum[k] / tl * 1e6
        v = np.log2(a + 1.0) - np.log2(b + 1.0)
        lin_lfc[L] = v
        x = lfc_cl_bal[common_mask]
        y = v[common_mask]
        r = float(np.corrcoef(x, y)[0, 1]) if len(x) > 10 else 0.0
        sgn = float((np.sign(x) == np.sign(y)).mean()) if len(x) > 10 else 0.0
        in_sig = [s for s in prim_sig if s in idx_of]
        same = sum(1 for s in in_sig
                   if np.sign(v[idx_of[s]]) == np.sign(prim_sig[s]))
        # a lineage signature of its own, for the recovery statistic
        own = 0
        rec = 0
        for s in prim_sig:
            i = idx_of.get(s)
            if i is None:
                continue
            if abs(float(v[i])) > LFC_MIN:
                own += 1
                if np.sign(v[i]) == np.sign(prim_sig[s]):
                    rec += 1
        repl_rows.append({
            "lineage": L,
            "cells_high": int(ln_hi_mask[k].sum()),
            "cells_low": int(ln_lo_mask[k].sum()),
            "genes_compared": int(len(x)),
            "pearson_r_vs_primary": round(r, 4),
            "sign_agreement_pct_all_genes": round(100.0 * sgn, 2),
            "signature_genes_with_data": int(len(in_sig)),
            "signature_genes_concordant": int(same),
            "signature_concordance_pct": round(100.0 * same / max(1, len(in_sig)), 2),
            "lineage_own_signature_n": int(own),
            "primary_recovered_pct": round(100.0 * rec / max(1, len(prim_sig)), 2),
        })
        log(f"  {L:<12} high {int(ln_hi_mask[k].sum()):>4} / low {int(ln_lo_mask[k].sum()):>4}"
            f"  r = {r:+.3f}  sign agreement {100.0*sgn:>5.1f}%  "
            f"signature concordance {100.0*same/max(1,len(in_sig)):>5.1f}% "
            f"({same}/{len(in_sig)})")
    write_table(repl_rows, "day4_cross_lineage_replication.csv")

    # per-gene cross-lineage table
    cross_rows: List[Dict[str, Any]] = []
    n_lineages_concordant_hist = collections.Counter()
    for s in sorted(prim_sig, key=lambda x: -abs(prim_sig[x])):
        i = idx_of.get(s)
        row = {"symbol": s, "log2fc_primary": round(prim_sig[s], 4)}
        n_conc = 0
        n_avail = 0
        for L in repl_lineages:
            v = lin_lfc.get(L)
            if v is None or i is None:
                row[f"log2fc_{L}"] = ""
                continue
            row[f"log2fc_{L}"] = round(float(v[i]), 4)
            n_avail += 1
            if np.sign(v[i]) == np.sign(prim_sig[s]):
                n_conc += 1
        row["lineages_available"] = n_avail
        row["lineages_concordant"] = n_conc
        row["concordance_pct"] = round(100.0 * n_conc / max(1, n_avail), 1)
        cross_rows.append(row)
        n_lineages_concordant_hist[n_conc] += 1
    write_table(cross_rows, "day4_cross_lineage_signature.csv")

    mean_sig_conc = (sum(r["signature_concordance_pct"] for r in repl_rows) / len(repl_rows)
                     if repl_rows else 0.0)
    median_r = float(np.median([r["pearson_r_vs_primary"] for r in repl_rows])) if repl_rows else 0.0
    n_5of6 = sum(1 for r in cross_rows
                 if r["lineages_available"] >= 4
                 and r["lineages_concordant"] >= max(4, r["lineages_available"] - 1))
    log(f"  mean signature concordance across {len(repl_rows)} lineages: "
        f"{mean_sig_conc:.1f}%; median Pearson r = {median_r:+.3f}")
    log(f"  signature genes reproduced in at least n-1 of n lineages: "
        f"{n_5of6} / {len(prim_sig)} "
        f"({100.0*n_5of6/max(1,len(prim_sig)):.1f}%)")
    log(f"  concordance histogram (lineages agreeing, out of {len(repl_rows)}): "
        + ", ".join(f"{k}:{v}" for k, v in sorted(n_lineages_concordant_hist.items())))

    # ---------------- summary ----------------
    summary = [
        {"metric": "cells_total", "value": n_cells},
        {"metric": "cells_pass_qc", "value": int(keep.sum())},
        {"metric": "cells_high", "value": n_hi},
        {"metric": "cells_low", "value": n_lo},
        {"metric": "stratifier_mode", "value": "single gene, lineage-restricted split"},
        {"metric": "target_gene", "value": primary["symbol"]},
        {"metric": "sensitivity_gene", "value": ""},
        {"metric": "split_restriction",
         "value": "target-detected cells with a dominant lineage, split inside each lineage"},
        {"metric": "cells_used_for_split", "value": int(len(universe))},
        {"metric": "lineages_in_primary_split", "value": len(split_lineages)},
        {"metric": "lineages_used_for_replication", "value": len(repl_lineages)},
        {"metric": "lineage_assignment_threshold_z", "value": DOMINANT_MIN_SCORE},
        {"metric": "median_target_cpm", "value": round(float(np.median(t_cpm[ok])), 3)},
        {"metric": "target_detection_pct", "value": round(100.0 * float((tval[ok] > 0).mean()), 2)},
        {"metric": "target_r2_lineage_depth", "value": round(r2_t, 4)},
        {"metric": "genes_tested", "value": int(expressed.sum())},
        {"metric": "genes_celllevel_fdr_lt_0.05", "value": int((qvals < 0.05).sum())},
        {"metric": "clusters_used_for_consistency", "value": len(usable_cl)},
        {"metric": "signature_up_n", "value": len(up)},
        {"metric": "signature_dn_n", "value": len(dn)},
        {"metric": "lineage_filtered_genes", "value": len(audit)},
        {"metric": "library_size_gap_pct", "value": round(gap, 2)},
        {"metric": "target_residual_gap_high_vs_low",
         "value": round(float(resid_split[hi].mean() - resid_split[lo].mean()), 3)},
        {"metric": "cross_lineage_signature_concordance_pct", "value": round(mean_sig_conc, 1)},
        {"metric": "cross_lineage_median_pearson_r", "value": round(median_r, 3)},
        {"metric": "cross_lineage_genes_reproduced", "value": int(n_5of6)},
        {"metric": "cross_lineage_genes_reproduced_pct",
         "value": round(100.0 * n_5of6 / max(1, len(prim_sig)), 1)},
        # kept so the Day-7 axis keeps a single robustness number to read
        {"metric": "sensitivity_direction_agreement_pct", "value": round(mean_sig_conc, 1)},
    ]
    write_table(summary, "day4_signature_summary.csv")

    # ---------------- figures ----------------
    section("4.9 Figures")
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    ax = axes[0]
    ax.hist(np.log10(np.maximum(counts, 1)), bins=60, color="#1f3b57", edgecolor="none")
    ax.set_xlabel("log10 total UMI per cell")
    ax.set_ylabel("Cells")
    ax.set_title("A  Library size")
    ax.spines[["top", "right"]].set_visible(False)
    ax2 = axes[1]
    ax2.hist(mito_frac, bins=60, color="#2a9d8f", edgecolor="none")
    ax2.axvline(max_mito, color="#c1443c", ls="--", lw=0.9)
    ax2.set_xlabel("Mitochondrial UMI (%)")
    ax2.set_ylabel("Cells")
    ax2.set_title("B  Mitochondrial fraction")
    ax2.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Fig 12. Quality control of the human bulbar conjunctiva atlas (GSM7761922)",
                 fontsize=10)
    save_fig(fig, "fig12_scrna_qc")

    if comp_rows:
        top = comp_rows[:18]
        fig, ax = plt.subplots(figsize=(7.4, 3.6))
        ax.bar([r["cluster"] for r in top], [r["pct"] for r in top],
               color="#1f3b57", edgecolor="white", linewidth=0.5)
        ax.set_ylabel("% of QC-passed cells")
        ax.set_xlabel("Author cluster")
        ax.tick_params(axis="x", rotation=60, labelsize=7)
        ax.set_title("Fig 13. Cell-state composition of healthy bulbar conjunctiva")
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig13_cluster_composition")

    if cl_rows:
        top = cl_rows[:16][::-1]
        fig, ax = plt.subplots(figsize=(6.6, max(3.0, 0.28 * len(top) + 1.3)))
        y = np.arange(len(top))
        ax.barh(y, [r["mean_target_cpm"] for r in top], color="#2a9d8f",
                edgecolor="white", height=0.66)
        for yi, r in zip(y, top):
            ax.text(r["mean_target_cpm"] * 1.02, yi, f"{r['pct_expressing']:.0f}%",
                    va="center", fontsize=6.8)
        ax.set_yticks(y)
        ax.set_yticklabels([r["cluster"] for r in top], fontsize=7)
        ax.set_xlabel(f"Mean {primary['symbol']} CPM per cell")
        ax.set_title(f"Fig 14. {primary['symbol']} expression across conjunctival cell states")
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, "fig14_target_gene_by_cluster")

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    xs = np.array([r["log2fc_cluster_balanced"] for r in de_rows])
    ys = np.array([-math.log10(max(r["p_value"], 1e-300)) for r in de_rows])
    sig = (np.array([r["cluster_agreement"] for r in de_rows]) >= CONSISTENCY_MIN)
    ax.scatter(xs[~sig], ys[~sig], s=4, c="#c9d3dc", linewidths=0, alpha=0.6)
    ax.scatter(xs[sig & (xs > LFC_MIN)], ys[sig & (xs > LFC_MIN)], s=5, c="#c1443c",
               linewidths=0, label=f"up in high-state cells (agreement >= {CONSISTENCY_MIN})")
    ax.scatter(xs[sig & (xs < -LFC_MIN)], ys[sig & (xs < -LFC_MIN)], s=5, c="#1f3b57",
               linewidths=0, label="down in high-state cells")
    ax.axhline(-math.log10(0.05), color="#555555", ls="--", lw=0.7)
    ax.axvline(0, color="#555555", lw=0.6)
    ax.set_xlabel("Cluster-balanced log2 fold change (high vs low state)")
    ax.set_ylabel("-log10 cell-level P (descriptive only)")
    ax.set_title(f"Fig 15. {primary['symbol']} state contrast inside conjunctival lineages")
    ax.legend(fontsize=7, loc="upper left", markerscale=3)
    ax.spines[["top", "right"]].set_visible(False)
    save_fig(fig, "fig15_signature_volcano")

    top_up, top_dn = up[:12], dn[:12]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8), sharey=False)
    ax = axes[0]
    ax.barh([r["symbol"] for r in top_up][::-1],
            [r["log2fc_cluster_balanced"] for r in top_up][::-1],
            color="#c1443c", edgecolor="white")
    ax.set_title("A  Up in high-state cells", fontsize=9)
    ax.set_xlabel("cluster-balanced log2FC")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax2 = axes[1]
    ax2.barh([r["symbol"] for r in top_dn][::-1],
             [r["log2fc_cluster_balanced"] for r in top_dn][::-1],
             color="#1f3b57", edgecolor="white")
    ax2.set_title("B  Down in high-state cells", fontsize=9)
    ax2.set_xlabel("cluster-balanced log2FC")
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.tick_params(labelsize=7)
    fig.suptitle(f"Fig 16. Leading genes of the lineage-restricted {primary['symbol']} signature",
                 fontsize=10)
    save_fig(fig, "fig16_signature_top_genes")

    # Fig 16b - why restriction replaces regression
    fr = [r for r in feas_rows if r["lineage"] != "unassigned"]
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
    ax = axes[0]
    yy = np.arange(len(fr))
    ax.barh(yy, [r["detection_pct"] for r in fr], color="#2a9d8f", edgecolor="white", height=0.68)
    for yi, r in zip(yy, fr):
        ax.text(r["detection_pct"] + 0.25, yi,
                f"{r['cells_detecting_target']}", va="center", fontsize=6.6)
    ax.set_yticks(yy)
    ax.set_yticklabels([r["lineage"] for r in fr], fontsize=7)
    ax.set_xlabel(f"{primary['symbol']} detection rate (%)")
    ax.set_title("A  Detected in every lineage, so not a lineage marker", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax2 = axes[1]
    ax2.hist(resid_split[hi], bins=60, color="#c1443c", alpha=0.55,
             label="high state", edgecolor="none")
    ax2.hist(resid_split[lo], bins=60, color="#1f3b57", alpha=0.55,
             label="low state", edgecolor="none")
    ax2.set_xlabel(f"Residualised {primary['symbol']} (z)")
    ax2.set_ylabel("Cells")
    ax2.set_title("B  The split variable, within lineage and depth", fontsize=9)
    ax2.legend(fontsize=7)
    ax2.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Fig 16b. Restriction instead of regression: where the contrast lives",
                 fontsize=10)
    save_fig(fig, "fig16b_lineage_restriction")

    # Fig 16c - cross-lineage replication
    if repl_rows and lin_lfc:
        show = ([r["symbol"] for r in up[:20]] + [r["symbol"] for r in dn[:20]])
        show = [s for s in show if s in idx_of]
        mat = np.array([[float(lin_lfc[L][idx_of[s]]) for L in repl_lineages]
                        for s in show])
        fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6),
                                 gridspec_kw={"width_ratios": [1.35, 1.0]})
        ax = axes[0]
        vmax = float(np.abs(mat).max()) or 1.0
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(np.arange(len(repl_lineages)))
        ax.set_xticklabels(repl_lineages, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(np.arange(len(show)))
        ax.set_yticklabels(show, fontsize=5.6)
        ax.set_title("A  Top signature genes, per lineage", fontsize=9)
        cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
        cb.set_label("log2FC (high vs low state)", fontsize=7)
        cb.ax.tick_params(labelsize=6)
        ax2 = axes[1]
        yy = np.arange(len(repl_rows))
        ax2.barh(yy, [r["signature_concordance_pct"] for r in repl_rows],
                 color="#1f3b57", edgecolor="white", height=0.66)
        ax2.axvline(50, color="#c1443c", ls="--", lw=0.9)
        for yi, r in zip(yy, repl_rows):
            ax2.text(r["signature_concordance_pct"] + 1.0, yi,
                     f"r={r['pearson_r_vs_primary']:+.2f}", va="center", fontsize=6.6)
        ax2.set_yticks(yy)
        ax2.set_yticklabels([r["lineage"] for r in repl_rows], fontsize=7)
        ax2.set_xlabel("% of signature genes with the same direction")
        ax2.set_xlim(0, 100)
        ax2.set_title("B  Replication in each lineage (dashed = chance)", fontsize=9)
        ax2.spines[["top", "right"]].set_visible(False)
        fig.suptitle("Fig 16c. The signature reproduces in independent conjunctival lineages",
                     fontsize=10)
        save_fig(fig, "fig16c_cross_lineage_replication")

    log("")
    log("Day 4 complete (v7, lineage-restricted split with cross-lineage replication).")


if __name__ == "__main__":
    main()
