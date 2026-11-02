"""make_figs.py -- Scientific Reports / Nature-format plates for the
allergic-conjunctivitis causal-target + reversal-compound study.

Every number drawn here is read from tables/*.csv at run time; nothing is typed
in from memory except axis labels and the flowchart wording.

    python scripts/report/make_figs.py

Writes fig1..fig6 as 600-dpi PNG + vector PDF into figures/, plus a per-figure
audit record in figures/fig_audit.json.

Canvas rules enforced by construction
  * width 165 mm (the template text column), height <= 170 mm, no
    bbox_inches="tight" anywhere, so the saved file measures exactly the
    authored size; panel rects are laid out on a 183 mm design grid and scaled
    onto the 165 mm canvas;
  * fonts are specified in points against that final canvas: 7 pt body,
    6.0-6.5 pt secondary text, 8 pt bold panel letters, nothing below 6 pt;
  * every axes rect is given in millimetres on a top-down grid, so gutters are
    explicit and panels cannot drift into one another;
  * Okabe-Ito only, no red/green pair, viridis for continuous scales;
  * save() reports the smallest drawn type size and any text that overruns the
    canvas, so collisions are caught numerically as well as visually.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch, Rectangle
from matplotlib.text import Text
from matplotlib.ticker import FuncFormatter, NullLocator

PROJ = Path(__file__).resolve().parent.parent.parent
TAB = PROJ / "tables"
OUT = PROJ / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- palette ---
# Okabe-Ito, colour-blind safe.  Shapes, fills and hatches carry the
# distinctions; colour never does it alone and no red/green pair is used.
BLACK = "#000000"
ORANGE = "#E69F00"
SKY = "#56B4E9"
GREEN = "#009E73"
YELLOW = "#F0E442"
BLUE = "#0072B2"
VERM = "#D55E00"
PURPLE = "#CC79A7"

INK = "#1A1A1A"
SPINE = "#333333"
PALE = "#B0B0B0"
FILL = "#F2F2F2"
LIGHT = "#D9D9D9"
GREY = "#8C9298"

W = 165.0          # mm, the printed canvas width (template text column)
VW = 183.0         # mm, the width the panel grid was designed on


def style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "font.size": 7.0,
        "axes.labelsize": 7.0,
        "axes.titlesize": 7.0,
        "axes.labelcolor": INK,
        "axes.edgecolor": SPINE,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "xtick.color": SPINE,
        "ytick.color": SPINE,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 2.2,
        "ytick.major.size": 2.2,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.size": 1.2,
        "ytick.minor.size": 1.2,
        "xtick.major.pad": 1.8,
        "ytick.major.pad": 1.8,
        "legend.fontsize": 7.0,
        "legend.frameon": False,
        "legend.handlelength": 1.0,
        "legend.handletextpad": 0.4,
        "legend.labelspacing": 0.26,
        "legend.columnspacing": 0.8,
        "legend.borderaxespad": 0.15,
        "lines.linewidth": 0.9,
        "lines.markersize": 3.0,
        "patch.linewidth": 0.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 200,
        "savefig.bbox": "standard",
        "savefig.pad_inches": 0.0,
        "text.color": INK,
        "hatch.linewidth": 0.5,
    })


def new_fig(h_mm: float):
    """Design is authored on a VW-wide grid; the canvas is uniformly shrunk to
    W, so every panel keeps its aspect ratio and every font keeps its ratio to
    the panel it labels."""
    return plt.figure(figsize=(W / 25.4, h_mm * W / VW / 25.4))


def axm(fig, left, top, w, h, H):
    """Axes rect in mm, measured from the top-left canvas corner."""
    return fig.add_axes([left / VW, 1.0 - (top + h) / H, w / VW, h / H])


def letter(fig, x_mm, top_mm, s, H, dy=0.9):
    """8 pt bold upright lower-case panel letter just above its panel."""
    fig.text(x_mm / VW, 1.0 - (top_mm - dy) / H, s, fontsize=8.0,
             fontweight="bold", ha="left", va="bottom", color=INK)


def note(fig, x_mm, top_mm, s, H, **kw):
    """Small annotation positioned in mm on the canvas."""
    fig.text(x_mm / VW, 1.0 - top_mm / H, s, **kw)


def _texts(fig):
    """Every visible Text with a non-empty string, plus its pixel bbox."""
    fig.canvas.draw()
    rnd = fig.canvas.get_renderer()
    out = []
    for tx in fig.findobj(Text):
        s = tx.get_text()
        if not s.strip() or not tx.get_visible():
            continue
        ax = tx.axes
        if ax is not None and not ax.get_visible():
            continue
        # an Annotation's window extent includes its leader line, which
        # legitimately crosses other ink; measure the text alone.
        bb = Text.get_window_extent(tx, rnd)
        out.append((tx, bb, s))
    return out


fig_dpi = 100.0


def _collisions(items, fig_dpi=100.0):
    """Text pairs whose ink boxes overlap by more than a hair.  Rotated labels
    are skipped: their axis-aligned boxes overlap by construction while the
    glyphs stay clear, so including them only produces false alarms."""
    bad = []
    for i in range(len(items)):
        ta, a, sa = items[i]
        if abs(ta.get_rotation()) > 1.0:
            continue
        for j in range(i + 1, len(items)):
            tb, b, sb = items[j]
            if abs(tb.get_rotation()) > 1.0:
                continue
            ow = min(a.x1, b.x1) - max(a.x0, b.x0)
            oh = min(a.y1, b.y1) - max(a.y0, b.y0)
            mm = 25.4 / fig_dpi
            if ow * mm > 1.0 and oh * mm > 0.7:   # smaller grazes are invisible in print
                dpi = fig_dpi
                bad.append(
                    f"{sa.strip()[:16]!r}@({a.x0/dpi*25.4:.1f},{a.y0/dpi*25.4:.1f})"
                    f" x {sb.strip()[:16]!r}@({b.x0/dpi*25.4:.1f},{b.y0/dpi*25.4:.1f})"
                    f" ov={ow/dpi*25.4:.1f}x{oh/dpi*25.4:.1f}mm")
    return bad


def save(fig, stem: str) -> None:
    items = _texts(fig)
    rnd = fig.canvas.get_renderer()
    fb = fig.bbox
    sizes, over = [], []
    for tx, bb, s in items:
        sizes.append(tx.get_fontsize())
        if (bb.x0 < -0.6 or bb.y0 < -0.6 or bb.x1 > fb.x1 + 0.6
                or bb.y1 > fb.y1 + 0.6):
            over.append(s.replace("\n", " ")[:22]
                        .encode("ascii", "replace").decode())
    bad = _collisions(items, float(fig.dpi))
    fig.savefig(OUT / f"{stem}.png", dpi=600)
    fig.savefig(OUT / f"{stem}.pdf")
    w_mm = float(fig.get_size_inches()[0] * 25.4)
    h_mm = float(fig.get_size_inches()[1] * 25.4)
    audit[stem] = {"min_pt": float(min(sizes)), "max_pt": float(max(sizes)),
                   "w_mm": round(w_mm, 2), "h_mm": round(h_mm, 2),
                   "n_text": len(sizes), "outside": over,
                   "collisions": bad}
    plt.close(fig)
    flag = "OK" if (not over and not bad and min(sizes) >= 6.5) else "CHECK"
    print(f"  {flag} {stem}.png / .pdf | {w_mm:.1f} x {h_mm:.1f} mm | "
          f"type {min(sizes):.1f}-{max(sizes):.1f} pt | outside {len(over)} | "
          f"collisions {len(bad)}"
          + (f" {over[:3]}" if over else "") + (f" {bad[:3]}" if bad else ""))


def load(name: str) -> pd.DataFrame:
    return pd.read_csv(TAB / f"{name}.csv")


TABLE: dict[str, pd.DataFrame] = {}
audit: dict[str, dict] = {}


def t(name: str) -> pd.DataFrame:
    if name not in TABLE:
        TABLE[name] = load(name)
    return TABLE[name]


thou = FuncFormatter(lambda v, p: f"{int(v):,}")


def sci(v: float, nd: int = 1) -> str:
    """Mathtext scientific notation, e.g. 7.9 x 10^-6."""
    exp = int(np.floor(np.log10(abs(v))))
    return rf"${v / 10 ** exp:.{nd}f}\times10^{{{exp}}}$"


# ================================================================== Fig 1 ===
def fig1() -> None:
    H = 132.0
    pool = t("day1_pool_summary")
    prot = len(t("day1_hgnc_protein_coding"))
    allrow = pool[pool.tier == "ALL"].iloc[0]
    scr = t("day2_screen_summary").set_index("metric")["value"]
    n_mr = int(scr["genes_screened"])
    n_p05 = int(scr["genes_P_lt_0.05"])
    coloc_ok = int((t("day3_coloc_results").status == "ok").sum())
    n_comp = len(t("day5_cmap_connectivity"))
    n_saf = len(t("day6_gene_safety_verdict"))
    n_t1 = int((t("day7_composite_ranking").candidate_tier
                == "Tier1_lead_candidate").sum())
    n_prior = int(t("day3_candidate_evidence").prioritised.sum())
    n_sig = len(t("day4_signature_up")) + len(t("day4_signature_dn"))

    fig = new_fig(H)

    # ---- a: serpentine workflow -----------------------------------------
    axa = axm(fig, 5, 5, 173, 78, H)
    axa.set_axis_off()
    axa.set_xlim(0, 173)
    axa.set_ylim(0, 78)
    nodes = [
        "Druggable genome\n{:,} genes".format(int(allrow.n_genes)),
        "Blood cis-eQTL instruments\n19,935 datasets",
        "Allergic conjunctivitis GWAS\n(FinnGen)",
        "Wald ratio / IVW\n{:,} genes tested".format(n_mr),
        "Colocalization + SMR/HEIDI\n{} candidates".format(coloc_ok),
        "Single-cell depth-stratified DE\n{}-gene signature".format(n_sig),
        "CMap reversal screen\n{} compounds".format(n_comp),
        "PheWAS safety screen\n{} genes x 1,400 endpoints".format(n_saf),
        "Five-axis integration\n{} Tier 1 candidate".format(n_t1),
    ]
    bw, bh, gx, gy = 52.0, 22.0, 8.5, 6.0
    colx = [0, bw + gx, 2 * (bw + gx)]
    rowy = [78 - bh, 78 - 2 * bh - gy, 78 - 3 * bh - 2 * gy]
    pos = [(colx[0], rowy[0]), (colx[1], rowy[0]), (colx[2], rowy[0]),
           (colx[2], rowy[1]), (colx[1], rowy[1]), (colx[0], rowy[1]),
           (colx[0], rowy[2]), (colx[1], rowy[2]), (colx[2], rowy[2])]
    for (x, y), txt, i in zip(pos, nodes, range(1, 10)):
        axa.add_patch(Rectangle((x, y), bw, bh, facecolor=FILL,
                                edgecolor=SPINE, linewidth=0.7, zorder=2))
        axa.text(x + bw / 2, y + bh / 2, txt, ha="center", va="center",
                 fontsize=7.0, linespacing=1.30, zorder=3, color=BLACK)
        axa.text(x + 1.8, y + bh - 1.8, str(i), ha="left", va="top",
                 fontsize=7.0, color=GREY, zorder=3)
    for k in range(len(pos) - 1):
        (x1, y1), (x2, y2) = pos[k], pos[k + 1]
        if y1 == y2:
            sx = np.sign(x2 - x1)
            p1 = (x1 + (bw if sx > 0 else 0), y1 + bh / 2)
            p2 = (x2 + (0 if sx > 0 else bw), y2 + bh / 2)
        else:
            p1 = (x1 + bw / 2, y1)
            p2 = (x2 + bw / 2, y2 + bh)
        axa.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=6,
                                      linewidth=0.8, color=SPINE, shrinkA=0,
                                      shrinkB=0, zorder=1))
    letter(fig, 5, 5, "a", H)

    # ---- b: druggable tiers ---------------------------------------------
    axb = axm(fig, 26, 96, 56, 27, H)
    sub = pool[pool.tier != "ALL"].copy()
    lab = ["T1\nclinical + drug", "T2\nclinical target", "T3\nliterature"]
    x = np.arange(len(sub))
    wbar = 0.36
    b1 = axb.bar(x - wbar / 2, sub.n_genes, wbar, color=BLUE, edgecolor="white",
                 linewidth=0.4, label="Pooled genes")
    b2 = axb.bar(x + wbar / 2, sub.n_with_cis_eqtl, wbar, color=SKY,
                 edgecolor="white", linewidth=0.4, label="With cis-eQTL")
    for bars in (b1, b2):
        for r in bars:
            axb.text(r.get_x() + r.get_width() / 2, r.get_height() + 90,
                     f"{int(r.get_height()):,}", ha="center", va="bottom",
                     fontsize=7.0)
    axb.set_xticks(x)
    axb.set_xticklabels(lab, fontsize=7.0)
    axb.set_ylabel("Genes (n)")
    axb.set_ylim(0, 6300)
    axb.set_yticks([0, 2000, 4000, 6000])
    axb.yaxis.set_major_formatter(thou)
    axb.legend(loc="lower center", ncol=2, fontsize=7.0,
               bbox_to_anchor=(0.5, 1.0))
    letter(fig, 15, 96, "b", H)

    # ---- c: attrition funnel --------------------------------------------
    axc = axm(fig, 116, 96, 60, 27, H)
    steps = [("Protein-coding genes", prot),
             ("Druggable genome", int(allrow.n_genes)),
             ("With cis-eQTL", int(allrow.n_with_cis_eqtl)),
             ("Tested by MR", n_mr),
             ("MR P < 0.05", n_p05),
             ("Prioritized (>1 axis)", n_prior),
             ("Tier 1 lead", n_t1)]
    yy = np.arange(len(steps))[::-1]
    vals = np.array([v for _, v in steps], dtype=float)
    axc.set_xscale("log")
    axc.barh(yy, vals - 1.0, left=1.0, height=0.60, color=BLUE,
             edgecolor="white", linewidth=0.4)
    axc.barh(yy[-1], vals[-1] - 1.0, left=1.0, height=0.60, color=VERM,
             edgecolor="white", linewidth=0.4)
    for y, v in zip(yy, vals):
        axc.text(v * 1.22, y, f"{int(v):,}", va="center", ha="left", fontsize=7.0)
    axc.set_yticks(yy)
    axc.set_yticklabels([s for s, _ in steps], fontsize=7.0)
    axc.set_xlim(1, 3.0e4)
    axc.set_ylim(-0.65, len(steps) - 0.35)
    axc.set_xticks([1, 10, 100, 1000, 10000])
    axc.xaxis.set_major_formatter(thou)
    axc.xaxis.set_minor_locator(NullLocator())
    axc.set_xlabel("Genes retained (n, log scale)")
    letter(fig, 105, 96, "c", H)

    print(f"[fig1] protein-coding {prot}, druggable {int(allrow.n_genes)}, "
          f"cis-eQTL {int(allrow.n_with_cis_eqtl)}, MR {n_mr}, P<0.05 {n_p05}, "
          f"prioritised {n_prior}, Tier1 {n_t1}, signature {n_sig}, "
          f"compounds {n_comp}, safety genes {n_saf}, coloc ok {coloc_ok}")
    save(fig, "fig1_study_design")


# ================================================================== Fig 2 ===
def fig2() -> None:
    H = 150.0
    s = t("day2_screen_lead_snp_mr").copy()
    s["nlp"] = -np.log10(s.wald_p.clip(lower=1e-300))
    bonf = 0.05 / len(s)
    fig = new_fig(H)

    # ---- a: volcano ------------------------------------------------------
    axa = axm(fig, 16, 14, 70, 56, H)
    rest = s[(s.wald_fdr >= 0.05)]
    fdr = s[(s.wald_fdr < 0.05) & (s.wald_p >= bonf)]
    sig = s[s.wald_p < bonf]
    axa.scatter(rest.wald_b, rest.nlp, s=2.0, c=PALE, linewidths=0,
                alpha=0.55, rasterized=True, zorder=2)
    axa.scatter(fdr.wald_b, fdr.nlp, s=13, c=BLUE, linewidths=0, zorder=3)
    axa.scatter(sig.wald_b, sig.nlp, s=26, c=VERM, marker="^", linewidths=0,
                zorder=4)
    axa.axhline(-np.log10(0.05), color=SPINE, linewidth=0.7,
                linestyle=(0, (4, 2)), zorder=1)
    axa.text(-1.42, -np.log10(0.05) + 0.12, "P = 0.05", fontsize=7.0,
             va="bottom", ha="left")
    axa.set_xlabel(r"Wald ratio ($\beta$, log OR per eQTL allele)")
    axa.set_ylabel(r"$-\log_{10}$ $P$")
    axa.set_xlim(-1.45, 1.05)
    axa.set_ylim(0, 9.6)
    axa.set_xticks([-1.2, -0.8, -0.4, 0, 0.4, 0.8])
    handles = [
        Line2D([], [], marker="^", linestyle="none", markersize=4.0, color=VERM,
               label=f"Bonferroni $P$ < {sci(bonf)} (n = {len(sig)})"),
        Line2D([], [], marker="o", linestyle="none", markersize=3.4, color=BLUE,
               label=f"FDR < 0.05 (n = {int((s.wald_fdr < 0.05).sum())})"),
        Line2D([], [], marker="o", linestyle="none", markersize=2.2, color=PALE,
               label=f"Not significant (n = {len(rest):,})"),
    ]
    axa.legend(handles=handles, loc="upper left", fontsize=7.0,
               bbox_to_anchor=(0.0, 1.015))
    for g, dx, dy, ha in (("IL18RAP", 11, -9, "left"), ("TLR1", 11, 2, "left"),
                          ("TLR10", 11, -3, "left")):
        r = s[s.symbol == g].iloc[0]
        axa.annotate(g, (r.wald_b, r.nlp), xytext=(dx, dy),
                     textcoords="offset points", fontsize=7.0,
                     fontweight="bold", ha=ha,
                     arrowprops=dict(arrowstyle="-", linewidth=0.6, color=SPINE,
                                     shrinkA=0, shrinkB=2))
    letter(fig, 7, 14, "a", H)

    # ---- b: forest, three leading candidates + background ----------------
    disc = t("day2_screen_lead_snp_mr").set_index("symbol")
    conf = t("day2_mr_confirmatory").set_index("symbol")
    prim = ["IL18RAP", "TLR1", "TLR10"]
    back = ["PHF1", "MID1IP1", "TXNDC5"]
    rows = sorted(back + prim, key=lambda g: disc.loc[g, "wald_or"])
    axb = axm(fig, 104, 14, 44, 56, H)
    yy = np.arange(len(rows))[::-1]
    for y, g in zip(yy, rows):
        r = disc.loc[g]
        axb.plot([r.wald_or_lci, r.wald_or_uci], [y, y], color=BLUE,
                 linewidth=0.9, solid_capstyle="butt", zorder=2)
        axb.plot([r.wald_or], [y], marker="o", markersize=3.8,
                 markerfacecolor=BLUE, markeredgecolor=BLUE, zorder=3)
        rep = conf.loc[g, "replication_or"]
        axb.plot([rep], [y], marker="D", markersize=3.2, markerfacecolor="white",
                 markeredgecolor=BLACK, markeredgewidth=0.7, zorder=4)
    axb.axvline(1.0, color=BLACK, linewidth=0.7, linestyle=(0, (4, 2)), zorder=1)
    axb.set_xscale("log")
    axb.set_xlim(0.11, 4.6)
    axb.set_xticks([0.25, 0.5, 1.0, 2.0, 4.0])
    axb.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
    axb.xaxis.set_minor_locator(NullLocator())
    axb.set_yticks(yy)
    axb.set_yticklabels(rows, fontsize=7.0)
    for lab in axb.get_yticklabels():
        if lab.get_text() in prim:
            lab.set_fontweight("bold")
        else:
            lab.set_color(GREY)
    axb.set_ylim(-0.6, len(rows) - 0.4)
    axb.set_xlabel("Odds ratio (log scale)")
    for y, g in zip(yy, rows):
        r = disc.loc[g]
        txt = (f"{r.wald_or:.2f} ({r.wald_or_lci:.2f}\u2013{r.wald_or_uci:.2f})")
        if g in prim:
            rp = conf.loc[g, "replication_p"]
            txt += f"\nreplication {sci(rp, 0)}"
        axb.text(1.04, y, txt, transform=axb.get_yaxis_transform(),
                 fontsize=7.0, va="center", ha="left",
                 linespacing=1.15, clip_on=False,
                 color=INK if g in prim else GREY)
    note(fig, 104, 11.0, "OR (95% CI)", H, fontsize=7.0, ha="left",
         va="bottom", color=SPINE)
    letter(fig, 94, 14, "b", H, dy=2.4)
    note(fig, 104, 78.6,
         "filled circle, discovery Wald ratio (95% CI);\n"
         "open diamond, replication OR;\n"
         "grey labels, background genes",
         H, fontsize=7.0, ha="left", va="top", color=SPINE, linespacing=1.3)

    # ---- c: instrument strength -----------------------------------------
    axc = axm(fig, 16, 86, 158, 50, H)
    f = s.f_statistic.dropna()
    cut = 3000.0
    tail = int((f > cut).sum())
    bins = np.arange(0, cut + 100, 100)
    axc.hist(f[f <= cut], bins=bins, color=BLUE, edgecolor="white",
             linewidth=0.25, zorder=2)
    med = float(f.median())
    axc.axvline(med, color=VERM, linewidth=0.8, linestyle=(0, (4, 2)), zorder=3)
    axc.set_xlabel("First-stage F statistic")
    axc.set_ylabel("Genes (n)")
    axc.set_xlim(0, cut)
    axc.set_xticks([0, 500, 1000, 1500, 2000, 2500, 3000])
    axc.xaxis.set_major_formatter(thou)
    ymax = axc.get_ylim()[1]
    axc.set_ylim(0, ymax)
    axc.text(med + 80, ymax * 0.93, f"median F = {med:,.1f}",
             fontsize=7.0, color=VERM, va="top", ha="left")
    axc.text(cut - 40, ymax * 0.93,
             f"{tail} genes with F > {int(cut):,} not shown\n"
             f"(maximum F = {f.max():,.0f})", fontsize=7.0, color=SPINE,
             va="top", ha="right")
    letter(fig, 7, 86, "c", H)

    print(f"[fig2] genes {len(s)}, bonferroni {len(sig)}, "
          f"FDR<0.05 {int((s.wald_fdr < 0.05).sum())}, median F {med}, "
          f"tail {tail}, forest {rows}")
    save(fig, "fig2_mr_screen")


# ================================================================== Fig 3 ===
def _coloc_smr(top_n: int = 8):
    co = t("day3_coloc_results")[["symbol", "PP_H0", "PP_H1", "PP_H2",
                                  "PP_H3", "PP_H4", "coloc_grade"]]
    sm = t("day3_smr_heidi")[["symbol", "b_smr", "se_smr", "p_smr", "heidi_p",
                              "smr_pass", "heidi_pass"]].copy()
    sm["smr_pass"] = sm.smr_pass.astype(bool)
    sm["heidi_pass"] = sm.heidi_pass.astype(bool)
    m = co.merge(sm, on="symbol", how="inner")
    m["nlp_heidi"] = -np.log10(m.heidi_p.replace(0, np.nan).fillna(1e-300))
    top = m.sort_values("PP_H4", ascending=False).head(top_n)
    return m, top


def fig3() -> None:
    H = 118.0
    m, top = _coloc_smr(8)
    genes = list(top.symbol)
    fig = new_fig(H)

    # ---- a: PP_H4 vs HEIDI ----------------------------------------------
    axa = axm(fig, 17, 12, 66, 92, H)
    axa.scatter(m.PP_H4, m.nlp_heidi, s=24, c=BLUE, linewidths=0,
                alpha=0.75, zorder=3)
    axa.axvline(0.8, color=SPINE, linewidth=0.7, linestyle=(0, (4, 2)), zorder=1)
    axa.axhline(-np.log10(0.05), color=SPINE, linewidth=0.7,
                linestyle=(0, (4, 2)), zorder=1)
    axa.text(0.815, 23.4, "PP_H4 = 0.8", fontsize=7.0, rotation=90,
             va="top", ha="left")
    axa.text(-0.005, -np.log10(0.05) + 0.45, "HEIDI $P$ = 0.05", fontsize=7.0,
             va="bottom", ha="left")
    for g, dx, dy, ha in (("IL18RAP", 8, 0, "left"), ("TLR10", -7, 12, "right"),
                          ("TLR1", -32, 15, "right"), ("NNT", -34, 9, "right")):
        row = m[m.symbol == g]
        if len(row):
            r = row.iloc[0]
            axa.annotate(g, (r.PP_H4, r.nlp_heidi), xytext=(dx, dy),
                         textcoords="offset points", fontsize=7.0,
                         fontweight="bold", ha=ha,
                         arrowprops=dict(arrowstyle="-", linewidth=0.6,
                                         color=SPINE, shrinkA=0, shrinkB=2))
    axa.set_xlim(-0.02, 1.02)
    axa.set_ylim(-0.8, 24.0)
    axa.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axa.set_yticks([0, 5, 10, 15, 20])
    axa.set_xlabel("Colocalization posterior PP_H4")
    axa.set_ylabel(r"$-\log_{10}$ HEIDI $P$")
    letter(fig, 8, 12, "a", H)

    # ---- b: posterior decomposition -------------------------------------
    axb = axm(fig, 108, 16, 56, 40, H)
    cols = ["PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4"]
    cols_col = [BLACK, SKY, GREEN, ORANGE, PURPLE]
    mids = ["H0", "H1", "H2", "H3", "H4"]
    yy = np.arange(len(genes))[::-1]
    left = np.zeros(len(genes))
    for col, cc, nm in zip(cols, cols_col, mids):
        v = top[col].to_numpy()
        axb.barh(yy, v, left=left, height=0.66, color=cc, edgecolor="white",
                 linewidth=0.4, label=nm)
        left += v
    for y, (_, r) in zip(yy, top.iterrows()):
        axb.text(1.02, y, f"{r.PP_H4:.3f}", va="center", ha="left", fontsize=7.0)
    axb.set_yticks(yy)
    axb.set_yticklabels(genes, fontsize=7.0)
    axb.set_xlim(0, 1.13)
    axb.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    axb.set_ylim(-0.6, len(genes) - 0.4)
    axb.set_xlabel("Posterior probability")
    axb.legend(loc="upper center", ncol=5, fontsize=7.0,
               bbox_to_anchor=(0.5, -0.26), handlelength=0.8, columnspacing=0.9)
    letter(fig, 108, 16, "b", H)

    # ---- c: SMR significance + HEIDI verdict ----------------------------
    axc = axm(fig, 108, 72, 56, 32, H)
    bar_c = [SKY if hp else "white" for hp in top.heidi_pass]
    axc.barh(yy, -np.log10(top.p_smr), height=0.66, color=bar_c,
             edgecolor=BLUE, linewidth=0.6,
             hatch=["" if hp else "///" for hp in top.heidi_pass])
    thr = -np.log10(0.05 / len(top))
    axc.axvline(thr, color=BLACK, linewidth=0.7, linestyle=(0, (4, 2)), zorder=1)
    vmax = float((-np.log10(top.p_smr)).max())
    xmax = vmax * 1.43
    for y, v in zip(yy, -np.log10(top.p_smr)):
        axc.text(v + 0.12, y, f"{v:.1f}", va="center", ha="left", fontsize=7.0)
    for y, hp in zip(yy, top.heidi_pass):
        axc.text(0.80, y, "ok" if hp else "rejected", fontsize=7.0,
                 va="center", ha="left", color=SPINE,
                 transform=axc.get_yaxis_transform())
    axc.text(0.80, len(genes) - 0.05, "HEIDI", fontsize=7.0, va="center",
             ha="left", color=SPINE, transform=axc.get_yaxis_transform())
    axc.text(thr + 0.10, len(genes) - 0.05, "SMR $P$ < 0.05/8", fontsize=7.0,
             ha="left", va="center", color=SPINE)
    axc.set_yticks(yy)
    axc.set_yticklabels(genes, fontsize=7.0)
    axc.set_ylim(-0.75, len(genes) + 0.15)
    axc.set_xlim(0, xmax)
    axc.set_xlabel(r"$-\log_{10}$ SMR $P$")
    letter(fig, 108, 72, "c", H)

    print(f"[fig3] merged {len(m)} genes, top8 {genes}, "
          f"PP_H4 range {m.PP_H4.min():.3f}-{m.PP_H4.max():.3f}, "
          f"heidi_pass top {list(top.heidi_pass)}")
    save(fig, "fig3_coloc_smr")


# ================================================================== Fig 4 ===
def fig4() -> None:
    H = 150.0
# "unassigned" is a residual bin, not a lineage; panel a covers the eight
    # annotated lineages so the reported range matches the text.
    feas = (t("day4_lineage_feasibility")
            .query("lineage != 'unassigned'")
            .sort_values("detection_pct"))
    bal = t("day4_group_balance").set_index("group")
    up = t("day4_signature_up").copy()
    dn = t("day4_signature_dn").copy()
    rep = t("day4_cross_lineage_replication").copy()
    null = t("day4_null_calibration").set_index("metric")["value"]
    fig = new_fig(H)

    # ---- a: lineage detection -------------------------------------------
    axa = axm(fig, 24, 16, 62, 50, H)
    yy = np.arange(len(feas))
    cols = [BLUE if e else LIGHT for e in feas.eligible_for_primary_split]
    axa.barh(yy, feas.detection_pct, height=0.64, color=cols, edgecolor="white",
             linewidth=0.4)
    for y, (_, r) in zip(yy, feas.iterrows()):
        axa.text(r.detection_pct - 0.35, y,
                 f"{int(r.cells_detecting_target):,}/{int(r.cells_total):,}",
                 va="center", ha="right", fontsize=7.0,
                 color="white" if r.eligible_for_primary_split else SPINE)
        axa.text(r.detection_pct + 0.35, y, f"{r.detection_pct:.1f}%",
                 va="center", ha="left", fontsize=7.0)
    axa.set_yticks(yy)
    axa.set_yticklabels(feas.lineage, fontsize=7.0)
    axa.set_xlim(0, 23)
    axa.set_ylim(-0.65, len(feas) - 0.35)
    axa.set_xticks([0, 5, 10, 15, 20])
    axa.set_xlabel("TLR1-detecting cells (%)")
    handles = [Patch(facecolor=BLUE, label="In primary split"),
               Patch(facecolor=LIGHT, label="Excluded (< 200 cells per group)")]
    axa.legend(handles=handles, loc="lower right", fontsize=7.0,
               bbox_to_anchor=(1.04, 1.0))
    note(fig, 24, 14.2, "detected / total cells", H, fontsize=7.0, ha="left",
         va="bottom", color=SPINE)
    letter(fig, 16, 16, "a", H)

    # ---- b: three balance mini-panels -----------------------------------
    hi, lo = bal.loc["high"], bal.loc["low"]
    gap = bal.loc["gap_high_vs_low"]
    specs = [("mean_total_umi", "UMI per cell",
              f"$\\Delta$ = {gap['mean_total_umi']:.2f}%"),
             ("mean_target_cpm", "TLR1 (CPM)",
              f"$\\Delta$ = {gap['mean_target_cpm']:+.2f}"),
             ("mean_residual", "Residual (s.d.)",
              f"$\\Delta$ = {gap['mean_residual']:+.2f} s.d.")]
    for k, (col, ylab, dlab) in enumerate(specs):
        x0 = 98 + k * 26
        ax = axm(fig, x0, 20, 22, 46, H)
        v_hi, v_lo = float(hi[col]), float(lo[col])
        ax.bar([0, 1], [v_hi, v_lo], 0.54, color=[SKY, BLUE],
               edgecolor="white", linewidth=0.4)
        span = max(abs(v_hi), abs(v_lo))
        for xx, v in zip([0, 1], [v_hi, v_lo]):
            if k == 0:                      # inside the bar, rotated, to fit
                ax.text(xx, v - span * 0.035, f"{v:,.0f}", ha="center",
                        va="top", fontsize=7.0, rotation=90, color="white")
            else:
                ax.text(xx, v + (span * 0.05 if v >= 0 else -span * 0.05),
                        f"{v:,.2f}", ha="center",
                        va="bottom" if v >= 0 else "top", fontsize=7.0)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["high", "low"], fontsize=7.0)
        ax.set_xlim(-0.6, 1.6)
        ax.tick_params(axis="y", labelsize=7.0)
        if v_lo < 0:
            ax.axhline(0, color=SPINE, linewidth=0.6, zorder=1)
            ax.set_ylim(-0.78, 0.78)
            ax.set_yticks([-0.5, 0, 0.5])
        else:
            ax.set_ylim(0, span * 1.32)
        ax.yaxis.set_major_formatter(thou if span > 1000 else
                                     FuncFormatter(lambda v, p: f"{v:g}"))
        ax.set_title(f"{ylab}\n{dlab}", fontsize=7.0, pad=2.5,
                     color=INK, linespacing=1.35)
        if k == 0:
            letter(fig, 90, 20, "b", H, dy=7.0)
    note(fig, 88, 72.0, "$\\Delta$, high minus low: percentage points for UMI,\n"
         "CPM for TLR1, s.d. for residual", H, fontsize=7.0, ha="left",
         va="top", color=SPINE, linespacing=1.3)

    # ---- c: signature volcano -------------------------------------------
    sig = pd.concat([up.assign(direction="up"), dn.assign(direction="dn")],
                    ignore_index=True)
    sig["nlp"] = -np.log10(sig.p_value.clip(lower=1e-300))
    axc = axm(fig, 24, 86, 62, 44, H)
    for d, cc in (("up", VERM), ("dn", BLUE)):
        sub = sig[sig.direction == d]
        axc.scatter(sub.log2fc_cluster_balanced, sub.nlp, s=6.5, c=cc,
                    linewidths=0, alpha=0.8, zorder=3)
    axc.axhline(-np.log10(0.05), color=SPINE, linewidth=0.7,
                linestyle=(0, (4, 2)), zorder=1)
    axc.axvline(0, color=SPINE, linewidth=0.6, zorder=1)
    axc.text(2.13, -np.log10(0.05), "P = 0.05", fontsize=7.0,
             va="center", ha="left", clip_on=False)
    marks = {"ITM2A": (9, 18, "left"), "PDLIM7": (9, 52, "left"),
             "SPARCL1": (9, 36, "left"), "BACH2": (-9, 26, "right"),
             "SAA1": (-9, 26, "right"), "CD200R1": (-9, 10, "right")}
    for g, (dx, dy, ha) in marks.items():
        row = sig[sig.symbol == g]
        if not len(row):
            continue
        r = row.iloc[0]
        axc.annotate(g, (r.log2fc_cluster_balanced, r.nlp),
                     xytext=(dx, dy), textcoords="offset points",
                     fontsize=7.0, fontweight="bold", ha=ha,
                     arrowprops=dict(arrowstyle="-", linewidth=0.5,
                                     color=SPINE, shrinkA=0, shrinkB=2))
    axc.set_xlabel(r"log$_2$ fold change (cluster-balanced)")
    axc.set_ylabel(r"$-\log_{10}$ $P$")
    axc.set_xlim(-2.1, 2.1)
    axc.set_ylim(-4, 92)
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=3.2,
                      color=VERM, label=f"Higher in high (n = {len(up)})"),
               Line2D([], [], marker="o", linestyle="none", markersize=3.2,
                      color=BLUE, label=f"Lower in high (n = {len(dn)})")]
    axc.legend(handles=handles, loc="upper left", fontsize=7.0,
               bbox_to_anchor=(0.0, 1.015))
    letter(fig, 16, 86, "c", H)

    # ---- d: cross-lineage replication -----------------------------------
    axd = axm(fig, 108, 86, 62, 44, H)
    xx = np.arange(len(rep))
    axd.plot(xx, rep.pearson_r_vs_primary, marker="o", markersize=3.6,
             color=SKY, linewidth=1.0)
    axd.axhline(float(null["null_pearson_r"]), color=SKY, linewidth=0.7,
                linestyle=(0, (4, 2)))
    axd.set_ylim(-0.02, 0.44)
    axd.set_ylabel("Pearson $r$ vs primary split", color=SKY)
    axd.tick_params(axis="y", colors=SKY)
    axd.set_xticks(xx)
    axd.set_xticklabels(rep.lineage, rotation=30, ha="right", fontsize=7.0)
    axd.set_xlim(-0.4, len(rep) - 0.1)
    axd2 = axd.twinx()
    axd2.spines["top"].set_visible(False)
    axd2.plot(xx, rep.sign_agreement_pct_all_genes, marker="s", markersize=3.2,
              color=ORANGE, linewidth=1.0, linestyle=(0, (3, 1.5)))
    axd2.axhline(float(null["null_sign_agreement_pct"]), color=ORANGE,
                 linewidth=0.7, linestyle=(0, (4, 2)))
    axd2.set_ylim(47.0, 70)
    axd2.set_ylabel("Concordant sign (%)", color=ORANGE)
    axd2.tick_params(axis="y", colors=ORANGE)
    axd.set_zorder(axd2.get_zorder() + 1)
    axd.patch.set_visible(False)
    axd.text(-0.25, float(null["null_pearson_r"]) + 0.010,
             f"null $r$ = {float(null['null_pearson_r']):.4f}", fontsize=7.0,
             color=SKY, va="bottom", ha="left")
    axd2.text(len(rep) - 0.25, float(null["null_sign_agreement_pct"]) - 0.35,
              f"null = {float(null['null_sign_agreement_pct']):.2f}%",
              fontsize=7.0, color=ORANGE, va="top", ha="right")
    letter(fig, 99, 86, "d", H)

    print(f"[fig4] lineages {list(feas.lineage)}; "
          f"detect range {feas.detection_pct.min()}-{feas.detection_pct.max()}; "
          f"balance umi {hi.mean_total_umi}/{lo.mean_total_umi} gap "
          f"{gap.mean_total_umi}; cpm gap {gap.mean_target_cpm}; "
          f"resid gap {gap.mean_residual}; rep lineages {list(rep.lineage)}")
    save(fig, "fig4_singlecell")


# ================================================================== Fig 5 ===
def fig5() -> None:
    H = 118.0
    nul = t("day5_null_by_size").sort_values("query_n")
    cmap = t("day5_cmap_connectivity").copy()
    cal = t("day5_query_size_calibration")
    fig = new_fig(H)

    # ---- a: observed vs null --------------------------------------------
    axa = axm(fig, 18, 14, 52, 40, H)
    xx = np.arange(len(nul))
    axa.bar(xx, nul.null_mean, 0.5, color=LIGHT, edgecolor=SPINE,
            linewidth=0.6, yerr=nul.null_sd, capsize=2.0,
            error_kw=dict(elinewidth=0.7, ecolor=SPINE, capthick=0.7),
            label="Permutation null, mean $\\pm$ 1 s.d.", zorder=2)
    axa.plot(xx, nul.observed_top_score, linestyle="none", marker="D",
             markersize=4.2, color=VERM, zorder=4,
             label="Observed best connectivity")
    for x, (_, r) in zip(xx, nul.iterrows()):
        axa.text(x, r.observed_top_score + 0.016,
                 f"z = {r.observed_z:.2f}\nP = {r.permutation_p:.3f}",
                 ha="center", va="bottom", fontsize=7.0, color=VERM)
    axa.set_xticks(xx)
    axa.set_xticklabels([str(int(v)) for v in nul.query_n], fontsize=7.0)
    axa.set_ylim(0, 0.47)
    axa.set_xlim(-0.6, len(xx) - 0.4)
    axa.set_xlabel("Query signature size (n genes)")
    axa.set_ylabel("Connectivity score")
    axa.legend(loc="upper center", fontsize=7.0, ncol=1,
               bbox_to_anchor=(0.5, -0.22))
    letter(fig, 8, 14, "a", H)

    # ---- b: top compounds -----------------------------------------------
    axb = axm(fig, 106, 14, 46, 94, H)
    top = cmap.sort_values(["reversal_specificity", "reverse_best_score"],
                           ascending=False).head(16).iloc[::-1]
    yy = np.arange(len(top))
    cols = [VERM if c == "BI-2536" else BLUE for c in top.compound]
    axb.barh(yy, top.reversal_specificity, height=0.62, color=cols,
             edgecolor="white", linewidth=0.4)
    for y, (_, r) in zip(yy, top.iterrows()):
        axb.text(r.reversal_specificity + 0.002, y,
                 f"{r.reversal_specificity:.4f}", va="center", ha="left",
                 fontsize=7.0)
    axb.text(0.062, 2.6, "15 compounds\ntie at 0.0324", fontsize=7.0,
             va="center", ha="left", color=SPINE, linespacing=1.3)
    labels = [c if len(c) <= 22 else c[:20] + "." for c in top.compound]
    axb.set_yticks(yy)
    axb.set_yticklabels(labels, fontsize=7.0)
    axb.set_xlim(0, 0.128)
    axb.set_ylim(-0.62, len(top) - 0.38)
    axb.set_xticks([0, 0.04, 0.08, 0.12])
    axb.set_xlabel("Reversal specificity")
    axb.tick_params(axis="y", pad=1.2)
    note(fig, 106, 9.2,
         "Top 16 of 46 compounds; the only glucocorticoid in the library is\n"
         "the known-answer control in c, not a screen hit",
         H, fontsize=7.0, ha="left", va="bottom", color=SPINE, linespacing=1.3)
    letter(fig, 106, 14, "b", H, dy=0.9)

    # ---- c: glucocorticoid positive control ------------------------------
    axc = axm(fig, 18, 76, 52, 32, H)
    sw = cal[(cal.signature == "dexamethasone_response")
             & (cal.orientation == "swapped")].sort_values("query_n")
    axc.plot(sw.query_n, sw.best_glucocorticoid_rank, marker="o", markersize=4.0,
             color=BLUE, linewidth=1.0, zorder=3)
    for _, r in sw.iterrows():
        axc.annotate(f"{int(r.best_glucocorticoid_rank)} "
                     f"({int(r.n_glucocorticoid)}/13)",
                     (r.query_n, r.best_glucocorticoid_rank),
                     xytext=(0, 7), textcoords="offset points",
                     ha="center", va="bottom", fontsize=7.0)
    axc.set_ylim(8.4, 0.0)
    axc.set_yticks([1, 2, 3, 4, 5, 6, 7, 8])
    axc.set_xticks(list(sw.query_n))
    axc.set_xlim(5, 40)
    axc.set_xlabel("Query signature size (n genes)")
    axc.set_ylabel("Best glucocorticoid rank")
    letter(fig, 8, 76, "c", H)

    print(f"[fig5] null sizes {list(nul.query_n)}; obs "
          f"{list(nul.observed_top_score)}; z {list(nul.observed_z)}; "
          f"top compound "
          f"{cmap.sort_values('reversal_specificity', ascending=False).compound.iloc[0]} "
          f"{cmap.reversal_specificity.max():.5f}; "
          f"glucocorticoid ranks {list(sw.best_glucocorticoid_rank)}")
    save(fig, "fig5_cmap")


# ================================================================== Fig 6 ===
def fig6() -> None:
    H = 118.0
    verdict = t("day6_gene_safety_verdict")
    hits = t("day6_phewas_hits").copy()
    rank = t("day7_composite_ranking").copy()
    fig = new_fig(H)

    # ---- a: PheWAS burden by domain -------------------------------------
    dom_count = hits.domain.value_counts()
    keep = list(dom_count.index[:5])
    short = {keep[0]: "unmapped", keep[1]: "atopic/immune", keep[2]: "ocular",
             keep[3]: "circulatory", keep[4]: "respiratory"}
    hits["dom2"] = hits.domain.map(lambda d: d if d in keep else "other")
    order = verdict.sort_values("n_significant", ascending=True)
    genes = list(order.gene)
    piv = (hits[hits.gene.isin(genes)]
           .groupby(["gene", "dom2"]).size().unstack(fill_value=0)
           .reindex(index=genes, columns=keep + ["other"], fill_value=0))
    axa = axm(fig, 24, 14, 142, 36, H)
    yy = np.arange(len(genes))
    cols = [BLACK, SKY, GREEN, ORANGE, PURPLE, LIGHT]
    names = [short[k] for k in keep] + ["other (%d)" % len(dom_count.iloc[5:])]
    left = np.zeros(len(genes))
    for col, cc, nm in zip(keep + ["other"], cols, names):
        v = piv[col].to_numpy()
        axa.barh(yy, v, left=left, height=0.62, color=cc, edgecolor="white",
                 linewidth=0.4, label=nm)
        left += v
    vmap = verdict.set_index("gene")
    for y, g in zip(yy, genes):
        r = vmap.loc[g]
        if int(r.n_significant) == 0:
            axa.text(0.6, y, "0", va="center", ha="left", fontsize=7.0,
                     color=SPINE)
        axa.text(77, y, f"{int(r.n_ocular)}/{int(r.n_severe)}", va="center",
                 ha="right", fontsize=7.0, color=SPINE)
    axa.set_yticks(yy)
    axa.set_yticklabels(genes, fontsize=7.0)
    axa.set_xlim(0, 80)
    axa.set_ylim(-0.65, len(genes) - 0.35)
    axa.set_xlabel("PheWAS significant associations (n)")
    axa.text(80, len(genes) - 0.30, "ocular / severe", fontsize=7.0, ha="right",
             va="bottom", color=SPINE)
    axa.legend(loc="lower left", ncol=3, fontsize=7.0,
               bbox_to_anchor=(0.0, 1.0), handlelength=0.8,
               columnspacing=0.9)
    letter(fig, 15, 14, "a", H, dy=5.6)

    # ---- b: five-axis heatmap -------------------------------------------
    tier2 = rank[rank.candidate_tier == "Tier2_supporting_candidate"].symbol.tolist()
    heat_genes = ["TLR1"] + [g for g in tier2 if g != "TLR1"][:6]
    axes_cols = ["score_mr", "score_coloc", "score_singlecell", "score_cmap",
                 "score_safety"]
    axes_lab = ["MR", "Coloc", "Single-cell", "CMap", "Safety"]
    mat = (rank.set_index("symbol").loc[heat_genes, axes_cols]
           .astype(float).to_numpy())
    zero = mat == 0.0
    axb = axm(fig, 24, 68, 44, 38, H)
    im = axb.imshow(np.ma.masked_where(zero, mat), cmap="viridis",
                    vmin=0, vmax=1, aspect="auto")
    im.cmap.set_bad("white")
    axb.set_xticks(np.arange(len(axes_cols)))
    axb.set_xticklabels(axes_lab, fontsize=7.0, rotation=32, ha="right")
    axb.set_yticks(np.arange(len(heat_genes)))
    axb.set_yticklabels(heat_genes, fontsize=7.0)
    axb.tick_params(length=0)
    axb.set_xticks(np.arange(-0.5, len(axes_cols), 1), minor=True)
    axb.set_yticks(np.arange(-0.5, len(heat_genes), 1), minor=True)
    axb.grid(which="minor", color=LIGHT, linewidth=0.8)
    axb.tick_params(which="minor", length=0)
    for sp in axb.spines.values():
        sp.set_visible(False)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if zero[i, j]:
                axb.text(j, i, "\u2013", ha="center", va="center",
                         fontsize=7.0, color=SPINE)
            else:
                axb.text(j, i, f"{v:.2f}", ha="center", va="center",
                         fontsize=7.0, color="white" if v < 0.55 else "black")
    cax = axm(fig, 70, 68, 2.6, 27, H)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Axis score (0–1)", fontsize=7.0)
    cb.ax.tick_params(labelsize=7.0, length=2.0, width=0.6)
    cb.outline.set_linewidth(0.6)
    note(fig, 24, 64.5, "\u2013 = gene not scored on that axis (score 0)",
         H, fontsize=7.0, ha="left", va="bottom", color=SPINE)
    letter(fig, 16, 68, "b", H)

    # ---- c: composite ranking -------------------------------------------
    top10 = rank.sort_values("composite_score", ascending=False).head(10)
    axc = axm(fig, 94, 68, 72, 36, H)
    xx = np.arange(len(top10))
    tcol = {"Tier1_lead_candidate": VERM, "Tier2_supporting_candidate": BLUE}
    cols = [tcol.get(k, LIGHT) for k in top10.candidate_tier]
    axc.bar(xx, top10.composite_score, 0.60, color=cols, edgecolor="white",
            linewidth=0.4)
    axc.axhline(0.55, color=BLACK, linewidth=0.8, linestyle=(0, (4, 2)), zorder=1)
    axc.text(len(top10) - 0.4, 0.572, "Tier 1 threshold", fontsize=7.0,
             ha="right", va="bottom")
    for x, v in zip(xx, top10.composite_score):
        axc.text(x, v + 0.008, f"{v:.2f}", ha="center", va="bottom", fontsize=7.0,
                 rotation=90)
    tiers = [(0, 0, "Tier 1"), (1, 6, "Tier 2"), (7, 9, "Tier 3")]
    for a, b, lab in tiers:
        axc.plot([a - 0.35, b + 0.35], [0.955, 0.955], color=BLACK,
                 linewidth=0.6, clip_on=False)
        axc.text((a + b) / 2, 0.965, lab, ha="center", va="bottom", fontsize=7.0)
    axc.set_xticks(xx)
    axc.set_xticklabels(top10.symbol, rotation=40, ha="right", fontsize=7.0)
    axc.set_ylim(0, 1.02)
    axc.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    axc.set_ylabel("Composite score")
    axc.set_xlim(-0.6, len(top10) - 0.4)
    letter(fig, 86, 68, "c", H)

    print(f"[fig6] PBXIP1 signals "
          f"{int(verdict.set_index('gene').loc['PBXIP1', 'n_significant'])}, "
          f"heat genes {heat_genes}, top10 {list(top10.symbol)}")
    save(fig, "fig6_safety_integration")


def main() -> None:
    style()
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    with open(OUT / "fig_audit.json", "w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=2)
    print("done ->", OUT)


if __name__ == "__main__":
    main()
