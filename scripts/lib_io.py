"""Shared project utilities: paths, IO, logging, plotting (SVG only)."""
from __future__ import annotations

import csv
import json
import math
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import rcParams  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "config" / "config.json"
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIG_SVG = ROOT / "figures" / "svg"
FIG_PNG = ROOT / "figures" / "png"
REPORT = ROOT / "report"
LOGS = ROOT / "logs"
INTER = DATA / "intermediate"

for _p in (DATA, TABLES, FIG_SVG, FIG_PNG, REPORT, LOGS, INTER):
    _p.mkdir(parents=True, exist_ok=True)

_CFG: Optional[Dict[str, Any]] = None


def cfg() -> Dict[str, Any]:
    global _CFG
    if _CFG is None:
        with open(CONFIG_PATH, encoding="utf-8") as fh:
            _CFG = json.load(fh)
    return _CFG


def token() -> str:
    p = ROOT / cfg()["api"]["token_file"]
    return p.read_text(encoding="ascii").strip()


# ---------------------------------------------------------------------------
# logging
# ---------------------------------------------------------------------------

_LOG_STATE: Dict[str, Any] = {"path": None}


def log_init(day: str) -> Path:
    p = LOGS / f"{day}_{datetime.now():%Y%m%d_%H%M%S}.log"
    _LOG_STATE["path"] = p
    return p


def log(msg: str) -> None:
    line = f"[{datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    p = _LOG_STATE.get("path")
    if p:
        with open(p, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def section(title: str) -> None:
    log("")
    log("=" * 78)
    log(title)
    log("=" * 78)


# ---------------------------------------------------------------------------
# tabular IO
# ---------------------------------------------------------------------------


def write_table(rows: Iterable[Dict[str, Any]], name: str, columns: Optional[Sequence[str]] = None) -> Path:
    rows = list(rows)
    path = TABLES / name
    if not rows:
        path.write_text("", encoding="utf-8")
        return path
    if columns is None:
        columns = list(rows[0].keys())
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(columns), extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt_cell(r.get(k)) for k in columns})
    return path


def _fmt_cell(v: Any) -> Any:
    if v is None:
        return ""
    if isinstance(v, float):
        if math.isnan(v):
            return ""
        if math.isinf(v):
            return "Inf" if v > 0 else "-Inf"
        return f"{v:.6g}"
    if isinstance(v, (list, tuple, set)):
        return ";".join(str(x) for x in v)
    return v


def read_table(name: str) -> List[Dict[str, str]]:
    p = TABLES / name
    if not p.exists():
        return []
    with open(p, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_json(obj: Any, name: str) -> Path:
    p = INTER / name
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=1), encoding="utf-8")
    return p


def read_json(name: str, default=None):
    p = INTER / name
    if not p.exists():
        return default
    return json.loads(p.read_text(encoding="utf-8"))


def cache_path(key: str) -> Path:
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in key)
    return INTER / "cache" / (safe + ".json")


def cache_get(key: str):
    p = cache_path(key)
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return None
    return None


def cache_put(key: str, obj) -> None:
    p = cache_path(key)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False), encoding="utf-8")


# ---------------------------------------------------------------------------
# plotting (SVG primary)
# ---------------------------------------------------------------------------

PALETTE = {
    "navy": "#1f3b57",
    "teal": "#2a9d8f",
    "amber": "#e9a13b",
    "red": "#c1443c",
    "grey": "#8d99ae",
    "light": "#dde5ec",
    "purple": "#6b4e8f",
    "green": "#4b7f52",
}

CAT_PALETTE = ["#1f3b57", "#2a9d8f", "#e9a13b", "#c1443c", "#6b4e8f",
               "#4b7f52", "#8d99ae", "#b5651d", "#3d6d8f", "#a05195"]


def style_init() -> None:
    rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "axes.edgecolor": "#3a3a3a",
        "axes.linewidth": 0.8,
        "axes.grid": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
        "legend.frameon": False,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "out",
        "ytick.direction": "out",
    })


def save_fig(fig, stem: str) -> List[Path]:
    """Save figure as SVG (primary deliverable) plus a PNG preview."""
    style_init()
    out = []
    p_svg = FIG_SVG / f"{stem}.svg"
    fig.savefig(p_svg, format="svg", bbox_inches="tight")
    out.append(p_svg)
    p_png = FIG_PNG / f"{stem}.png"
    fig.savefig(p_png, format="png", dpi=220, bbox_inches="tight")
    out.append(p_png)
    plt.close(fig)
    log(f"  figure -> {p_svg.relative_to(ROOT)}")
    return out


def fig_placeholder_note(ax, text: str) -> None:
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes,
            fontsize=9, color="#666666", wrap=True)
    ax.set_axis_off()


def banner(day: str, title: str) -> None:
    log("")
    log("#" * 78)
    log(f"# {day}: {title}")
    log("#" * 78)
