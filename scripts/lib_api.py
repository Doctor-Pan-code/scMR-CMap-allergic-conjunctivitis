"""API client wrappers with on-disk caching for the scMR + CMap project."""
from __future__ import annotations

import json
import math
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence

import requests
import urllib3

from lib_io import cache_get, cache_put, cfg, log, token

urllib3.disable_warnings()

_SESSION = requests.Session()
_SESSION.verify = False
_SESSION.headers.update({"User-Agent": "scMR-CMap-project/1.0"})


def _sleep() -> None:
    time.sleep(float(cfg()["api"].get("request_sleep", 1.0)))


def _post(url: str, *, json_body: Optional[dict] = None, data=None, files=None,
          headers: Optional[dict] = None, timeout: Optional[int] = None,
          retries: Optional[int] = None) -> requests.Response:
    r_max = retries if retries is not None else int(cfg()["api"].get("max_retry", 4))
    to = timeout or int(cfg()["api"].get("timeout", 300))
    last = None
    for attempt in range(r_max):
        try:
            r = _SESSION.post(url, json=json_body, data=data, files=files,
                              headers=headers, timeout=to)
            if r.status_code in (429, 500, 502, 503, 504):
                last = r
                time.sleep(2.0 * (attempt + 1))
                continue
            return r
        except Exception as exc:  # network hiccup
            last = exc
            time.sleep(2.0 * (attempt + 1))
    if isinstance(last, requests.Response):
        return last
    raise RuntimeError(f"POST failed after {r_max} attempts: {url}: {last}")


def _get(url: str, *, params: Optional[dict] = None, headers: Optional[dict] = None,
         timeout: Optional[int] = None, retries: Optional[int] = None) -> requests.Response:
    r_max = retries if retries is not None else int(cfg()["api"].get("max_retry", 4))
    to = timeout or int(cfg()["api"].get("timeout", 300))
    last = None
    for attempt in range(r_max):
        try:
            r = _SESSION.get(url, params=params, headers=headers, timeout=to)
            if r.status_code in (429, 500, 502, 503, 504):
                last = r
                time.sleep(2.0 * (attempt + 1))
                continue
            return r
        except Exception as exc:
            last = exc
            time.sleep(2.0 * (attempt + 1))
    if isinstance(last, requests.Response):
        return last
    raise RuntimeError(f"GET failed after {r_max} attempts: {url}: {last}")


# ---------------------------------------------------------------------------
# OpenGWAS
# ---------------------------------------------------------------------------


def _og_headers() -> Dict[str, str]:
    return {"Authorization": "Bearer " + token(), "Content-Type": "application/json"}


def og_gwasinfo(ids: Sequence[str]) -> List[Dict[str, Any]]:
    base = cfg()["api"]["opengwas_base"]
    key = "gwasinfo_" + "_".join(sorted(ids)[:6]) + f"_n{len(ids)}"
    hit = cache_get(key)
    if hit is not None:
        return hit
    r = _post(f"{base}/gwasinfo", json_body={"id": list(ids)}, headers=_og_headers())
    if r.status_code != 200:
        raise RuntimeError(f"gwasinfo {r.status_code}: {r.text[:200]}")
    out = r.json()
    cache_put(key, out)
    return out


def og_tophits(ids: Sequence[str], pval: float = 5e-8) -> List[Dict[str, Any]]:
    base = cfg()["api"]["opengwas_base"]
    r = _post(f"{base}/tophits", json_body={"id": list(ids), "pval": float(pval)},
              headers=_og_headers())
    if r.status_code != 200:
        raise RuntimeError(f"tophits {r.status_code}: {r.text[:200]}")
    return r.json()


def og_associations(ids: Sequence[str], variants: Sequence[str]) -> List[Dict[str, Any]]:
    """Cartesian product query. Constraint: len(ids) * len(variants) <= 64."""
    base = cfg()["api"]["opengwas_base"]
    if len(ids) * len(variants) > 64:
        raise ValueError("OpenGWAS associations requires N(id)*N(variant) <= 64")
    r = _post(f"{base}/associations", json_body={"id": list(ids), "variant": list(variants)},
              headers=_og_headers())
    if r.status_code != 200:
        raise RuntimeError(f"associations {r.status_code}: {r.text[:200]}")
    return r.json()


def og_associations_batched(variant: str, ids: Sequence[str], batch: int = 60) -> Dict[str, Dict[str, Any]]:
    """One variant against many datasets, respecting the 64 product cap.

    The OpenGWAS /api/associations endpoint enforces N(dataset) * N(variant) <= 64
    and answers a violation with a bare 401 {"message": "Unknown error."}, which
    reads like an authentication failure but is not. A batch of exactly 64 is
    therefore one id away from silently losing a whole chunk (Day 6 lost 64 of
    1400 endpoints for TLR10 that way), so a failed chunk is retried in halves
    down to a single dataset before it is given up on.
    """
    out: Dict[str, Dict[str, Any]] = {}
    ids = list(dict.fromkeys(ids))

    def fetch(chunk: Sequence[str]) -> None:
        try:
            rows = og_associations(list(chunk), [variant])
        except Exception as exc:
            if len(chunk) > 1:
                mid = len(chunk) // 2
                log(f"    associations batch of {len(chunk)} failed ({exc}); "
                    f"retrying as {mid}+{len(chunk)-mid}")
                fetch(chunk[:mid])
                fetch(chunk[mid:])
                return
            log(f"    associations failed for single dataset {chunk[0]}: {exc}")
            return
        for row in rows:
            out[row["id"]] = row
        _sleep()

    for i in range(0, len(ids), batch):
        fetch(ids[i:i + batch])
    return out


def og_associations_pair(variants: Sequence[str], dataset: str, batch: int = 60) -> Dict[str, Dict[str, Any]]:
    """One dataset against many variants, respecting the 64 product cap."""
    out: Dict[str, Dict[str, Any]] = {}
    variants = list(dict.fromkeys(variants))
    for i in range(0, len(variants), batch):
        chunk = variants[i:i + batch]
        try:
            rows = og_associations([dataset], chunk)
        except Exception as exc:
            log(f"    associations batch failed ({i}-{i+len(chunk)}): {exc}")
            continue
        for row in rows:
            out[row["rsid"]] = row
        _sleep()
    return out


# ---------------------------------------------------------------------------
# Ensembl REST
# ---------------------------------------------------------------------------


def ensembl_overlap_variants(chrom: str, start: int, end: int, assembly: str = "grch38") -> List[Dict[str, Any]]:
    host = cfg()["api"]["ensembl_rest"] if assembly == "grch38" else cfg()["api"]["ensembl_rest_grch37"]
    url = f"{host}/overlap/region/human/{chrom}:{start}-{end}"
    r = _get(url, params={"feature": "variation"}, headers={"Content-Type": "application/json"},
             timeout=180)
    if r.status_code != 200:
        return []
    try:
        return r.json()
    except Exception:
        return []


def ensembl_ld_pairwise(rs_a: str, rs_b: str, population: str = "1000GENOMES:phase_3:EUR") -> Optional[Dict[str, Any]]:
    url = f"{cfg()['api']['ensembl_rest']}/ld/human/pairwise/{rs_a}/{rs_b}"
    r = _get(url, params={"population_name": population},
             headers={"Content-Type": "application/json"}, timeout=90)
    if r.status_code != 200:
        return None
    try:
        data = r.json()
    except Exception:
        return None
    if not isinstance(data, list):
        return None
    for rec in data:
        if rec.get("population_name") == population:
            return rec
    return data[0] if data else None


def ensembl_ld_matrix(rsids: Sequence[str], population: str = "1000GENOMES:phase_3:EUR") -> Optional[Dict[str, Any]]:
    url = f"{cfg()['api']['ensembl_rest']}/ld/human/pairwise"
    r = _post(url, json_body={"variants": list(rsids), "population_name": population},
              headers={"Content-Type": "application/json"}, timeout=120)
    if r.status_code != 200:
        return None
    try:
        return r.json()
    except Exception:
        return None


def ensembl_lookup_symbol(symbol: str, assembly: str = "grch37") -> Optional[Dict[str, Any]]:
    host = cfg()["api"]["ensembl_rest"] if assembly == "grch38" else cfg()["api"]["ensembl_rest_grch37"]
    key = f"ensembl_sym_{assembly}_{symbol}"
    hit = cache_get(key)
    if hit is not None:
        return hit
    r = _get(f"{host}/lookup/symbol/homo_sapiens/{symbol}",
             headers={"Content-Type": "application/json"}, timeout=60)
    if r.status_code != 200:
        return None
    j = r.json()
    cache_put(key, j)
    return j


def ensembl_lookup_id(ensg: str, assembly: str = "grch37") -> Optional[Dict[str, Any]]:
    host = cfg()["api"]["ensembl_rest"] if assembly == "grch38" else cfg()["api"]["ensembl_rest_grch37"]
    key = f"ensembl_id_{assembly}_{ensg}"
    hit = cache_get(key)
    if hit is not None:
        return hit
    r = _get(f"{host}/lookup/id/{ensg}", headers={"Content-Type": "application/json"}, timeout=60)
    if r.status_code != 200:
        return None
    j = r.json()
    cache_put(key, j)
    return j


# ---------------------------------------------------------------------------
# L1000CDS2 (CMap connectivity)
# ---------------------------------------------------------------------------


def cmap_query(up: Sequence[str], dn: Sequence[str], aggravate: bool = False,
               db_version: Optional[str] = None) -> Dict[str, Any]:
    cfgc = cfg()["cmap"]
    body = {
        "data": {"upGenes": list(up), "dnGenes": list(dn)},
        "config": {
            "aggravate": bool(aggravate),
            "searchMethod": "geneSet",
            "share": False,
            "combination": False,
            "db-version": db_version or cfgc["db_version"],
        },
    }
    r = _post(cfg()["api"]["l1000cds2"], json_body=body, timeout=300)
    if r.status_code != 200:
        raise RuntimeError(f"l1000cds2 {r.status_code}: {r.text[:200]}")
    return r.json()


# ---------------------------------------------------------------------------
# Enrichr
# ---------------------------------------------------------------------------


def enrichr_add_list(genes: Sequence[str], description: str = "scMR_CMap") -> Optional[int]:
    r = _post(cfg()["api"]["enrichr_addlist"],
              files={"list": (None, "\n".join(genes)), "description": (None, description)},
              timeout=120)
    if r.status_code != 200:
        return None
    try:
        return int(r.json()["userListId"])
    except Exception:
        return None


def enrichr_enrich(user_list_id: int, libraries, timeout: int = 240) -> Dict[str, Any]:
    if isinstance(libraries, (list, tuple)):
        libraries = ",".join(libraries)
    r = _get(cfg()["api"]["enrichr_enrich"],
             params={"userListId": user_list_id, "backgroundType": libraries},
             timeout=timeout, retries=1)
    if r.status_code != 200:
        return {}
    try:
        return r.json()
    except Exception:
        return {}
