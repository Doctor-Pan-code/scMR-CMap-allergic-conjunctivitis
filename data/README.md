# Upstream data

The five upstream inputs are public but are **not** redistributed in this repository: they total
about 330 MB, each has its own licence and update cycle, and two of them are re-issued regularly
under the same URL, so a frozen copy would go stale silently.

Download them into this directory (`data/`) before running Day 1.

---

## 1. HGNC complete set → `data/hgnc_complete_set.txt`

```bash
curl -L -o data/hgnc_complete_set.txt \
  https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
```

Tab-separated, one row per HGNC record. The pipeline keeps `locus_group == "protein-coding gene"`
and records that carry an `ensembl_gene_id`; that filter yields the 19,255-gene universe. The
archive is stamped with a download date; if the file is re-fetched later, the count can shift by a
few dozen genes as symbols are merged or withdrawn, so record the date you used.

Alternative source (same content, different host):
<https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=md_eg_id&col=md_ensembl_id&status=Approved&status_opt=2&locus_group=protein-coding+gene&fmt=tsv&hdr=1&submit=submit>

## 2. DGIdb genes → `data/dgidb_genes.tsv`
## 3. DGIdb interactions → `data/dgidb_interactions.tsv`

```bash
curl -L -o data/dgidb_genes.tsv        https://dgidb.org/data/latest/genes.tsv
curl -L -o data/dgidb_interactions.tsv https://dgidb.org/data/latest/interactions.tsv
```

`genes.tsv` carries the gene-level claim counts per source database; `interactions.tsv` carries the
drug-gene pairs. Both are the `latest` release, which DGIdb rebuilds periodically — record the
download date. Tier assignment depends only on which sources appear for a gene, so a rebuild
changes tier counts only for genes whose source support changed.

## 4. OpenGWAS dataset catalogue → `data/opengwas_catalog.json`

This one is produced by the pipeline, not downloaded by hand: `scripts/day1_druggable_pool.py`
queries the OpenGWAS API for every `eqtl-a-ENSG*` dataset and writes the response to
`data/opengwas_catalog.json`. It needs a token first (see below). If you already have the file, the
script reuses it instead of re-querying.

To fetch it standalone:

```bash
curl -s -H "Authorization: Bearer $OPENGWAS_TOKEN" \
  https://api.opengwas.io/api/gwasinfo > data/opengwas_catalog.json
```

## 5. Single-cell atlas GSE242382 → `data/geo/GSE242382/`

Four files from the human bulbar conjunctiva atlas (GSM7761922, 4 pooled healthy biopsies):

```
GSM7761922_LX237_matrix.mtx.gz      290,974,030 B
GSM7761922_LX237_features.tsv.gz        333,437 B
GSM7761922_LX237_barcodes.tsv.gz        130,938 B
GSM7761922_LX237_clusters.tsv.gz      1,601,527 B
```

From GEO:

```bash
mkdir -p data/geo/GSE242382
cd data/geo/GSE242382
for f in matrix.mtx features.tsv barcodes.tsv clusters.tsv; do
  curl -L -O "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7761nnn/GSM7761922/suppl/GSM7761922_LX237_${f}.gz"
done
```

Or from the series page: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242382>

`day4_scrna_signature.py` streams the `.mtx.gz` in chunks — it never expands the file to disk and
never holds the full matrix in memory, so the 290 MB archive is read in place.

---

## The L1000 landmark list is included

`data/l1000_landmark_genes.txt` (14,631 B, 978 genes) is small, has no licence encumbrance, and is
needed by Day 5, so it ships with the repository rather than being fetched.

---

## OpenGWAS token

Register at <https://api.opengwas.io> (free for academic use), then:

```bash
cp config/opengwas_token.txt.example config/opengwas_token.txt
# paste the JWT on a single line, no quotes
```

`config/opengwas_token.txt` is excluded by `.gitignore`. Never commit it.

---

## Verifying the inputs

```bash
python - <<'PY'
from pathlib import Path
need = {
    "data/hgnc_complete_set.txt": 16_982_321,
    "data/dgidb_genes.tsv": 4_356_295,
    "data/dgidb_interactions.tsv": 12_178_745,
    "data/opengwas_catalog.json": 20_917_301,
    "data/l1000_landmark_genes.txt": 14_631,
    "data/geo/GSE242382/GSM7761922_LX237_matrix.mtx.gz": 290_974_030,
    "data/geo/GSE242382/GSM7761922_LX237_features.tsv.gz": 333_437,
    "data/geo/GSE242382/GSM7761922_LX237_barcodes.tsv.gz": 130_938,
    "data/geo/GSE242382/GSM7761922_LX237_clusters.tsv.gz": 1_601_527,
}
for name, size in need.items():
    p = Path(name)
    got = p.stat().st_size if p.exists() else None
    flag = "ok" if got == size else ("MISSING" if got is None else "size differs")
    print(f"{flag:14s} {name}  {got}")
PY
```

Sizes are those of the archive used for the manuscript. A newer `latest` release will differ and
that is expected; only the single-cell files must match byte for byte, because the cell counts in
the paper (27,600 cells, 18,007 genes) depend on that exact archive.
