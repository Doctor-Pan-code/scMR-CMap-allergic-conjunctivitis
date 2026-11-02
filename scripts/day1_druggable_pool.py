"""Day 1 - Druggable genome universe construction.

Goal
----
Build the analysis starting point: a protein-coding, drug-annotated gene
universe that is (a) instrumentable by cis-eQTL in OpenGWAS/eQTLGen and
(b) linked to known pharmacological agents.

Inputs
------
data/hgnc_complete_set.txt     HGNC full set (symbol, ensembl id, locus group)
data/dgidb_genes.tsv           DGIdb gene claims with source databases
data/dgidb_interactions.tsv    DGIdb drug-gene interactions
data/opengwas_catalog.json     full OpenGWAS dataset catalogue

Outputs
-------
tables/day1_hgnc_protein_coding.csv
tables/day1_druggable_pool.csv
tables/day1_pool_summary.csv
tables/day1_druggable_sources.csv
figures/svg/fig01_druggable_source_composition.svg
figures/svg/fig02_gene_pool_funnel.svg
figures/svg/fig03_pool_by_chromosome.svg
"""
from __future__ import annotations

import collections
import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lib_io import (DATA, ROOT, TABLES, banner, cfg, log, log_init, save_fig,
                    section, write_table)

CLINICAL_SOURCES = {
    "DrugBank", "ChEMBL", "GuideToPharmacology", "TTD", "Pharos", "dGene",
    "FoundationOneGenes", "TdgClinicalTrial", "NCI", "CarisMolecularIntelligence",
    "ClearityFoundationBiomarkers", "ClearityFoundationClinicalTrial",
}
LITERATURE_SOURCES = {
    "HopkinsGroom", "HingoraniCasas", "RussLampel", "DGI", "GO", "Wolf",
    "OncoKB", "CancerCommons", "MyCancerGenome", "CGI", "CIViC",
}
IDENTITY_SOURCES = {"NCBI", "HGNC", "Ensembl", "MESH", "UMLS", "Entrez"}


def load_hgnc() -> list[dict]:
    path = DATA / "hgnc_complete_set.txt"
    out = []
    with open(path, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            if row.get("status") != "Approved":
                continue
            if row.get("locus_group") != "protein-coding gene":
                continue
            ensg = (row.get("ensembl_gene_id") or "").strip()
            sym = (row.get("symbol") or "").strip()
            if not sym or not ensg:
                continue
            out.append({
                "symbol": sym,
                "ensembl_gene_id": ensg,
                "hgnc_id": row.get("hgnc_id", ""),
                "name": row.get("name", ""),
                "location": row.get("location", ""),
                "gene_group": row.get("gene_group", ""),
                "alias_symbol": row.get("alias_symbol", ""),
                "prev_symbol": row.get("prev_symbol", ""),
                "entrez_id": row.get("entrez_id", ""),
                "uniprot_ids": row.get("uniprot_ids", ""),
            })
    return out


def load_dgidb_sources() -> dict[str, set[str]]:
    path = DATA / "dgidb_genes.tsv"
    by_symbol: dict[str, set[str]] = collections.defaultdict(set)
    by_claim: dict[str, set[str]] = collections.defaultdict(set)
    with open(path, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            src = (row.get("source_db_name") or "").strip()
            if not src:
                continue
            gname = (row.get("gene_name") or "").strip()
            claim = (row.get("gene_claim_name") or "").strip()
            if gname and gname != "NULL":
                by_symbol[gname].add(src)
            if claim:
                by_claim[claim].add(src)
    for k, v in by_claim.items():
        by_symbol[k] |= v
    return by_symbol


def load_dgidb_interactions() -> tuple[dict[str, int], dict[str, set[str]]]:
    path = DATA / "dgidb_interactions.tsv"
    counts: dict[str, int] = collections.Counter()
    drugs: dict[str, set[str]] = collections.defaultdict(set)
    with open(path, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            g = (row.get("gene_name") or "").strip()
            if not g or g == "NULL":
                continue
            counts[g] += 1
            d = (row.get("drug_name") or "").strip()
            if d and d != "NULL":
                drugs[g].add(d)
    return dict(counts), {k: v for k, v in drugs.items()}


def main() -> None:
    log_init("day1")
    banner("Day 1", "Druggable genome universe construction")
    conf = cfg()

    section("1.1 HGNC protein-coding universe")
    hgnc = load_hgnc()
    log(f"approved protein-coding genes with Ensembl ID: {len(hgnc)}")
    write_table(hgnc, "day1_hgnc_protein_coding.csv")

    section("1.2 OpenGWAS cis-eQTL dataset availability")
    catalog = json.loads((DATA / "opengwas_catalog.json").read_text(encoding="utf-8"))
    eqtl_ids = {k[len("eqtl-a-"):]: k for k in catalog if k.startswith("eqtl-a-")}
    log(f"eqtl-a datasets in OpenGWAS: {len(eqtl_ids)}")

    section("1.3 DGIdb annotation")
    src_map = load_dgidb_sources()
    n_drugs, drug_names = load_dgidb_interactions()
    log(f"DGIdb genes with >=1 source claim: {len(src_map)}")
    log(f"DGIdb genes with >=1 recorded interaction: {len(n_drugs)}")

    section("1.4 Build tiered druggable pool")
    pool = []
    for g in hgnc:
        sym = g["symbol"]
        ensg = g["ensembl_gene_id"]
        srcs = src_map.get(sym, set()) | src_map.get(ensg, set())
        if not srcs:
            # try aliases
            for alias in (g["alias_symbol"] or "").split("|"):
                alias = alias.strip()
                if alias and alias in src_map:
                    srcs |= src_map[alias]
        clin = sorted(srcs & CLINICAL_SOURCES)
        lit = sorted(srcs & LITERATURE_SOURCES)
        ident = sorted(srcs & IDENTITY_SOURCES)
        n_int = n_drugs.get(sym, 0)
        if not (clin or lit or n_int):
            continue
        if clin and n_int >= 1:
            tier = "T1_clinical_with_drug"
        elif clin:
            tier = "T2_clinical_target"
        elif lit or n_int:
            tier = "T3_literature_or_experimental"
        else:
            continue
        eqtl_id = eqtl_ids.get(ensg)
        pool.append({
            **g,
            "dgidb_sources": ";".join(sorted(srcs)),
            "clinical_sources": ";".join(clin),
            "n_sources": len(srcs),
            "n_drug_interactions": n_int,
            "n_drugs": len(drug_names.get(sym, set())),
            "example_drugs": ";".join(sorted(drug_names.get(sym, set()))[:8]),
            "druggable_tier": tier,
            "eqtl_dataset": eqtl_id or "",
            "has_cis_eqtl": bool(eqtl_id),
            "chromosome": _chrom(g["location"]),
        })
    log(f"druggable pool size: {len(pool)}")
    tier_counts = collections.Counter(p["druggable_tier"] for p in pool)
    for t, c in sorted(tier_counts.items()):
        log(f"  {t}: {c}")
    n_inst = sum(1 for p in pool if p["has_cis_eqtl"])
    log(f"  with cis-eQTL instrument source: {n_inst}")

    write_table(pool, "day1_druggable_pool.csv")

    section("1.5 Pool summary")
    summary = []
    for t in sorted(tier_counts):
        sub = [p for p in pool if p["druggable_tier"] == t]
        summary.append({
            "tier": t,
            "n_genes": len(sub),
            "n_with_cis_eqtl": sum(1 for p in sub if p["has_cis_eqtl"]),
            "median_n_drugs": _median([p["n_drugs"] for p in sub]),
            "median_n_interactions": _median([p["n_drug_interactions"] for p in sub]),
        })
    summary.append({
        "tier": "ALL",
        "n_genes": len(pool),
        "n_with_cis_eqtl": n_inst,
        "median_n_drugs": _median([p["n_drugs"] for p in pool]),
        "median_n_interactions": _median([p["n_drug_interactions"] for p in pool]),
    })
    write_table(summary, "day1_pool_summary.csv")
    for s in summary:
        log(f"  {s['tier']:<28} genes={s['n_genes']:<6} instrumentable={s['n_with_cis_eqtl']}")

    section("1.6 Druggable source composition")
    src_counter = collections.Counter()
    for p in pool:
        for s in p["dgidb_sources"].split(";"):
            if s:
                src_counter[s] += 1
    src_rows = [{"source_db": k, "n_genes": v,
                 "class": "clinical" if k in CLINICAL_SOURCES else
                          ("literature" if k in LITERATURE_SOURCES else
                           ("identity" if k in IDENTITY_SOURCES else "other"))}
                for k, v in src_counter.most_common()]
    write_table(src_rows, "day1_druggable_sources.csv")

    section("1.7 Figures")
    # Fig 01: source composition
    top = src_rows[:14]
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    names = [r["source_db"] for r in top][::-1]
    vals = [r["n_genes"] for r in top][::-1]
    colors = {"clinical": "#1f3b57", "literature": "#2a9d8f",
              "identity": "#8d99ae", "other": "#c9d3dc"}
    bars = ax.barh(names, vals, color=[colors[r["class"]] for r in top][::-1],
                   edgecolor="white", linewidth=0.6)
    for b, v in zip(bars, vals):
        ax.text(v + max(vals) * 0.01, b.get_y() + b.get_height() / 2, f"{v:,}",
                va="center", fontsize=7.5)
    ax.set_xlabel("Number of genes in the tiered druggable pool")
    ax.set_title("Fig 1. DGIdb source databases supporting the druggable gene pool")
    ax.spines[["top", "right"]].set_visible(False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=colors[k], label=k) for k in
                       ["clinical", "literature", "identity", "other"]],
              loc="lower right", fontsize=7)
    save_fig(fig, "fig01_druggable_source_composition")

    # Fig 02: funnel
    stages = ["HGNC protein-coding\n(approved, with ENSG)",
              "Any DGIdb evidence",
              "Tiered druggable pool",
              "Instrumentable by\ncis-eQTL (OpenGWAS)"]
    values = [len(hgnc),
              len({p["symbol"] for p in pool}) + max(0, len(src_map) - len(hgnc)),
              len(pool),
              n_inst]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    y = np.arange(len(stages))[::-1]
    ax.barh(y, values, color=["#8d99ae", "#6b8fa8", "#2a9d8f", "#1f3b57"],
            edgecolor="white", height=0.62)
    for yi, v in zip(y, values):
        ax.text(v * 1.02, yi, f"{v:,}", va="center", fontsize=8.5)
    ax.set_yticks(y)
    ax.set_yticklabels(stages, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("Gene count (log scale)")
    ax.set_title("Fig 2. Attrition from the protein-coding genome to the instrumentable drug-target pool")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(1, max(values) * 4)
    save_fig(fig, "fig02_gene_pool_funnel")

    # Fig 03: pool by chromosome
    chrom_counts = collections.Counter(p["chromosome"] for p in pool)
    order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    xs = [c for c in order if c in chrom_counts]
    ys = [chrom_counts[c] for c in xs]
    fig, ax = plt.subplots(figsize=(8.0, 3.4))
    ax.bar(xs, ys, color="#1f3b57", edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Druggable genes")
    ax.set_title("Fig 3. Chromosomal distribution of the instrumentable druggable pool")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)
    save_fig(fig, "fig03_pool_by_chromosome")

    log("")
    log("Day 1 complete.")


def _chrom(location: str) -> str:
    if not location:
        return "NA"
    s = location.split("q")[0].split("p")[0]
    return s or "NA"


def _median(vals):
    v = sorted(x for x in vals if x is not None)
    if not v:
        return 0
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


if __name__ == "__main__":
    main()
