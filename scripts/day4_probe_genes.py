"""Probe - choose the Day 4 stratification gene on evidence instead of by default.

Why
---
Day 4 v4 stratified cells by TLR1 and produced a signature whose leading up genes
were melanocyte and macrophage markers and whose leading down genes were T/NK
markers. That is the signature of cell identity, not of a transcriptional state:
TLR1 is itself a myeloid lineage gene, so "TLR1 high vs low" re-orders cell types
even after cluster and depth matching.

Before redesigning the day, the candidate genes that Day 3 prioritised are scored
on the two properties that matter for stratification in this atlas:
  detection rate   - too sparse and the split degenerates into detected/undetected
  lineage R^2      - how much of the gene's per-cell variance is explained by
                     lineage marker scores; high R^2 means the split is a lineage
                     split and cannot be used as a state contrast

Outputs
-------
tables/day4_probe_cells.csv     per-cell QC, lineage scores, candidate gene CPM
tables/day4_probe_gene_scores.csv  per-gene detection, lineage R^2, verdict
"""
from __future__ import annotations

import collections
import gzip
import math
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

import lib_stats as st
from lib_io import DATA, TABLES, log, log_init, read_table, section, write_table

GEO = DATA / "geo" / "GSE242382"

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
}


def open_mtx(path: Path):
    fh = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    while True:
        line = fh.readline()
        if not line:
            raise RuntimeError("unexpected EOF before matrix dimensions")
        if line.startswith("%"):
            continue
        n_genes, n_cells, n_nz = (int(x) for x in line.split())
        return fh, n_genes, n_cells, n_nz


def load_features() -> List[Tuple[str, str]]:
    out = []
    with gzip.open(GEO / "GSM7761922_LX237_features.tsv.gz", "rt", encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                out.append((p[0], p[1]))
    return out


def load_barcodes() -> List[str]:
    with gzip.open(GEO / "GSM7761922_LX237_barcodes.tsv.gz", "rt", encoding="utf-8") as fh:
        return [l.strip() for l in fh if l.strip()]


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


def main() -> None:
    log_init("day4probe")
    section("Probe 1 - candidate genes and lineage panels in GSE242382")

    feats = load_features()
    barcodes = load_barcodes()
    clusters = load_clusters()
    sym_of = [s for _, s in feats]
    ens_of = [e for e, _ in feats]
    idx_of: Dict[str, int] = {}
    for i, s in enumerate(sym_of):
        idx_of.setdefault(s, i)

    cand_rows = [r for r in read_table("day3_candidate_evidence.csv") if r.get("prioritised") == "True"]
    cand = [(r["symbol"], r["ensembl_gene_id"], float(r["coloc_PP_H4"] or 0), float(r["screen_p"] or 1))
            for r in cand_rows]
    log(f"prioritised candidates from Day 3: {len(cand)} -> {', '.join(c[0] for c in cand)}")

    goi: Dict[int, str] = {}          # gene row index -> symbol
    for sym, ens, _, _ in cand:
        i = None
        if ens and ens in ens_of:
            i = ens_of.index(ens)
        elif sym in idx_of:
            i = idx_of[sym]
        if i is not None:
            goi[i] = sym
    panel_gene_idx: Dict[str, List[int]] = {}
    for pname, genes in PANELS.items():
        ids = [idx_of[g] for g in genes if g in idx_of]
        panel_gene_idx[pname] = ids
    panel_of_gene: Dict[int, str] = {}
    for pname, ids in panel_gene_idx.items():
        for i in ids:
            panel_of_gene.setdefault(i, pname)

    log(f"genes of interest matched to features: {len(goi)}")
    log(f"panel genes matched: {sum(len(v) for v in panel_gene_idx.values())}")

    n_cells = len(barcodes)
    counts = np.zeros(n_cells, dtype=np.int64)
    ngenes = np.zeros(n_cells, dtype=np.int32)
    mitoc = np.zeros(n_cells, dtype=np.int64)
    mito_idx = {i for i, s in enumerate(sym_of) if s.startswith("MT-")}
    goi_keys = sorted(goi)
    goi_pos = {g: k for k, g in enumerate(goi_keys)}
    goi_val = np.zeros((len(goi_keys), n_cells), dtype=np.int32)
    panel_names = sorted(panel_gene_idx)
    panel_val = np.zeros((len(panel_names), n_cells), dtype=np.int64)
    panel_pos = {p: k for k, p in enumerate(panel_names)}

    path = GEO / "GSM7761922_LX237_matrix.mtx.gz"
    fh, g_n, c_n, n_nz = open_mtx(path)
    log(f"matrix: {g_n} genes x {c_n} cells x {n_nz:,} non-zeros")
    t0 = time.time()
    step = max(1, n_nz // 10)
    nxt = step
    read = 0
    for line in fh:
        p = line.split()
        if len(p) != 3:
            continue
        gi = int(p[0]) - 1
        ci = int(p[1]) - 1
        v = int(p[2])
        if ci >= n_cells:
            continue
        counts[ci] += v
        ngenes[ci] += 1
        if gi in mito_idx:
            mitoc[ci] += v
        k = goi_pos.get(gi)
        if k is not None:
            goi_val[k, ci] = v
        pn = panel_of_gene.get(gi)
        if pn is not None:
            panel_val[panel_pos[pn], ci] += v
        read += 1
        if read >= nxt:
            log(f"  {read:,}/{n_nz:,}, {time.time()-t0:.0f}s")
            nxt += step
    fh.close()
    log(f"pass complete: {read:,} entries in {time.time()-t0:.0f}s")

    mito_frac = np.divide(mitoc, np.maximum(counts, 1)) * 100.0
    keep = (ngenes >= 200) & (mito_frac <= 20.0) & (counts >= 500)
    log(f"QC: {int(keep.sum())} / {n_cells} cells")

    cpm = np.divide(1.0, np.maximum(counts, 1)) * 1e6
    goi_cpm = goi_val * cpm[None, :]

    # lineage scores: mean log1p CPM of the panel, then z-scored across kept cells
    lib = counts.astype(np.float64)
    scores = {}
    for pname in panel_names:
        ids = panel_gene_idx[pname]
        s = np.zeros(n_cells)
        for i in ids:
            pass
        scores[pname] = s
    # panel_val already holds summed counts per panel; convert to CPM of panel size
    for pname in panel_names:
        k = panel_pos[pname]
        s = panel_val[k].astype(np.float64) / np.maximum(lib, 1.0) * 1e6
        s = np.log1p(s)
        mu = s[keep].mean()
        sd = s[keep].std() or 1.0
        scores[pname] = (s - mu) / sd
    Z = np.column_stack([scores[p] for p in panel_names])[keep]
    log(f"lineage score matrix: {Z.shape[0]} cells x {Z.shape[1]} panels")

    # per-gene: detection, lineage R^2
    gene_rows = []
    for k, gi in enumerate(goi_keys):
        sym = goi[gi]
        y = np.log1p(goi_cpm[k][keep])
        det = float((goi_val[k][keep] > 0).mean())
        # OLS of y on [1, Z]
        A = np.column_stack([np.ones(Z.shape[0]), Z])
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = y - A @ beta
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1.0 - float((resid ** 2).sum()) / ss_tot if ss_tot > 0 else 0.0
        # which panel carries the signal
        cors = {}
        for pname in panel_names:
            zz = scores[pname][keep]
            if zz.std() == 0:
                cors[pname] = 0.0
                continue
            cors[pname] = float(np.corrcoef(y, zz)[0, 1]) if y.std() > 0 else 0.0
        top = sorted(cors.items(), key=lambda kv: -abs(kv[1]))[:2]
        rec = next(c for c in cand if c[0] == sym)
        gene_rows.append({
            "symbol": sym,
            "ensembl_gene_id": rec[1],
            "coloc_PP_H4": rec[2],
            "screen_p": rec[3],
            "pct_detected": round(100.0 * det, 2),
            "mean_cpm": round(float(goi_cpm[k][keep].mean()), 3),
            "median_cpm": round(float(np.median(goi_cpm[k][keep])), 3),
            "lineage_r2": round(r2, 4),
            "top_panel_1": top[0][0], "corr_1": round(top[0][1], 3),
            "top_panel_2": top[1][0], "corr_2": round(top[1][1], 3),
        })
    gene_rows.sort(key=lambda r: (r["lineage_r2"], -r["coloc_PP_H4"]))
    write_table(gene_rows, "day4_probe_gene_scores.csv")

    section("Probe 2 - ranking")
    log(f"{'symbol':<10} {'H4':>6} {'det%':>6} {'lineageR2':>9} {'top panel':<12} {'corr':>6}")
    for r in gene_rows:
        log(f"{r['symbol']:<10} {r['coloc_PP_H4']:>6.3f} {r['pct_detected']:>6.2f} "
            f"{r['lineage_r2']:>9.4f} {r['top_panel_1']:<12} {r['corr_1']:>6.3f}")

    cell_rows = []
    for j in range(n_cells):
        if not keep[j]:
            continue
        row = {
            "barcode": barcodes[j],
            "cluster": clusters.get(barcodes[j], "unassigned"),
            "total_counts": int(counts[j]),
            "n_genes": int(ngenes[j]),
            "mito_pct": round(float(mito_frac[j]), 3),
        }
        for pname in panel_names:
            row[f"score_{pname}"] = round(float(scores[pname][j]), 4)
        for k, gi in enumerate(goi_keys):
            row[f"cpm_{goi[gi]}"] = round(float(goi_cpm[k][j]), 4)
        cell_rows.append(row)
    write_table(cell_rows, "day4_probe_cells.csv")
    log(f"per-cell probe table: {len(cell_rows)} cells, {len(cell_rows[0])} columns")


if __name__ == "__main__":
    main()
