#!/usr/bin/env python3
"""
NTC GO GSEA heatmaps — v5
- Rows ordered by hierarchical clustering (NO dendrogram drawn, NO theme grouping)
- Chemical color bar flush on top + vertical separators + compact chemical legend
- No title; short term labels (same cap for A5 and article); tight right margin
"""
import pandas as pd
import numpy as np
import re, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]   # repository root (this file lives in Scripts/figures/)
INPUT_CSV = ROOT / "Output" / "GSEA_NTC" / "ALL_FILES_NTC_GO_GSEA_SIG.csv"
OUT_DIR   = str(ROOT / "Output" / "Figures")
FDR_THRESH = 0.05
TOP_N      = 25
MAXLEN     = 30          # same short cap for both figures
DPI        = 300
os.makedirs(OUT_DIR, exist_ok=True)

CHEM_COLORS = {
    "2-4D": "#e41a1c", "CAFF": "#ff7f00", "CYCL": "#984ea3",
    "GLYP": "#4daf4a", "IMID": "#377eb8", "RACI": "#a65628",
    "TEBU": "#f781bf", "THIR": "#999999",
}

def shorten_term(s, maxlen=MAXLEN):
    s = str(s).strip()
    for a, b in [
        ("positive regulation of", "pos. reg."), ("negative regulation of", "neg. reg."),
        ("regulation of", "reg."), ("involved in", "in"),
        ("via plasma membrane", "via PM"), ("cell adhesion molecules", "CAMs"),
        ("signaling pathway", "signaling"), ("morphogenesis", "morph."),
        ("establishment or maintenance of", "est./maint."),
        ("polymerization or depolymerization", "polym./depolym."),
        ("extracellular matrix", "ECM"), ("epithelial to mesenchymal transition", "EMT"),
        ("epithelial-mesenchymal transition", "EMT"), ("biosynthetic process", "biosynth."),
        ("development", "dev."), ("differentiation", "diff."), ("organization", "org."),
        ("specification", "spec."), ("canonical", "canon."), ("receptor", "rec."),
    ]:
        s = s.replace(a, b)
    if len(s) > maxlen:
        s = s[:maxlen-1].rstrip() + "…"
    return s

# ── Load & parse ──────────────────────────────────────────────────────────────
df = pd.read_csv(INPUT_CSV)
df.columns = df.columns.str.strip()

def parse_file(f):
    f = f.strip(); chem = f.split("_")[0]
    m_s = re.search(r"_S(\d+)_", f); m_t = re.search(r"_T(\d+)_", f)
    return chem, int(m_s.group(1)) if m_s else 0, int(m_t.group(1)) if m_t else 0

df[["chem", "conc", "tp"]] = pd.DataFrame(df["file"].apply(parse_file).tolist(), index=df.index)
df["score"] = np.sign(df["NES"]) * -np.log10(df["p.adjust"].clip(lower=1e-300))
df["sig"]   = df["p.adjust"] < FDR_THRESH

df_sorted = df.sort_values(["chem", "conc", "tp"])
ci = df_sorted[["file", "chem", "conc", "tp"]].drop_duplicates()
ci["col_id"] = ci["chem"] + " S" + ci["conc"].astype(str) + "·T" + ci["tp"].astype(str)
col_order = ci["col_id"].tolist()
df_sorted["col_id"] = df_sorted["file"].map(dict(zip(ci["file"], ci["col_id"])))

def build_matrices(term_filter=None):
    ps = df_sorted.pivot_table(index="Description", columns="col_id", values="score", aggfunc="first").reindex(columns=col_order).fillna(0.0)
    pg = df_sorted.pivot_table(index="Description", columns="col_id", values="sig", aggfunc="first").reindex(columns=col_order).fillna(False).astype(bool)
    if term_filter is not None:
        ps = ps.loc[ps.index.isin(term_filter)]; pg = pg.loc[pg.index.isin(term_filter)]
    mask = pg.any(axis=1)
    return ps.loc[mask], pg.loc[mask]

def cluster_order(score_df):
    if score_df.shape[0] < 3:
        return list(score_df.index)
    dist = np.nan_to_num(pdist(score_df.values, metric="euclidean"), nan=0.0)
    return [score_df.index[i] for i in leaves_list(linkage(dist, method="ward"))]

# ── Draw ──────────────────────────────────────────────────────────────────────
def draw(score_df, sig_df, out_path, figsize, row_fs, col_fs, legend_fs, dot_size, anno_bar_h=0.6):
    nrow, ncol = score_df.shape
    labels = [shorten_term(d) for d in score_df.index]

    fig, ax = plt.subplots(figsize=figsize, dpi=DPI)
    vmax = max(np.nanmax(np.abs(score_df.values)), 0.5)
    im = ax.imshow(score_df.values, aspect="auto", cmap="RdBu_r",
                   norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), interpolation="nearest")

    S = sig_df.values
    for i in range(nrow):
        for j in range(ncol):
            if S[i, j]:
                ax.plot(j, i, "o", color="black", markersize=dot_size, markeredgewidth=0)

    # chemical bar flush on top
    bar_top = -0.5 - anno_bar_h
    for j, colid in enumerate(col_order):
        ax.add_patch(mpatches.Rectangle((j - 0.5, bar_top), 1, anno_bar_h,
                     facecolor=CHEM_COLORS.get(colid.split(" ")[0], "#666"), edgecolor="none", clip_on=False))
    ax.set_xlim(-0.5, ncol - 0.5)
    ax.set_ylim(nrow - 0.5, bar_top)

    # vertical separators between chemicals
    prev = None
    for j, colid in enumerate(col_order):
        ch = colid.split(" ")[0]
        if prev is not None and ch != prev:
            ax.axvline(j - 0.5, color="white", linewidth=1.3)
        prev = ch

    short_labels = [c.split(" ", 1)[1] if " " in c else c for c in col_order]
    ax.set_xticks(range(ncol)); ax.set_xticklabels(short_labels, rotation=90, fontsize=col_fs)
    ax.set_yticks(range(nrow)); ax.set_yticklabels(labels, fontsize=row_fs)
    ax.tick_params(axis="both", which="both", length=0)

    # compact colorbar + chemical legend tight against the heatmap right edge
    cax = ax.inset_axes([1.006, 0.0, 0.015, 0.55])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("sign(NES) · −log₁₀(FDR)", fontsize=legend_fs)
    cbar.ax.tick_params(labelsize=legend_fs - 1)

    chems_present = sorted(set(c.split(" ")[0] for c in col_order))
    patches = [mpatches.Patch(facecolor=CHEM_COLORS[ch], label=ch) for ch in chems_present]
    patches.append(plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="black", markersize=4, label="FDR < 0.05"))
    ax.legend(handles=patches, loc="upper left", bbox_to_anchor=(1.006, 1.0),
              fontsize=legend_fs, framealpha=0, edgecolor="none",
              title="Chemical", title_fontsize=legend_fs + 1,
              handlelength=1.0, handletextpad=0.4, labelspacing=0.25, borderpad=0.2)

    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white", pad_inches=0.02)
    if out_path.endswith(".png"):
        fig.savefig(out_path.replace(".png", ".pdf"), dpi=DPI, bbox_inches="tight", facecolor="white", pad_inches=0.02)
    plt.close(fig)
    print("  saved:", os.path.basename(out_path), f"({nrow} rows x {ncol} cols)")

# ── 1) A5 top-25 ───────────────────────────────────────────────────────────────
sig_only = df[df["sig"]]
term_stats = sig_only.groupby("Description").agg(n_sig=("sig", "sum"),
              mean_abs_nes=("NES", lambda x: np.mean(np.abs(x)))).sort_values(["n_sig", "mean_abs_nes"], ascending=False)
top_terms = term_stats.index[:TOP_N].tolist()
ps, pg = build_matrices(term_filter=top_terms)
order = cluster_order(ps); ps, pg = ps.loc[order], pg.loc[order]
print("A5 top-25:")
draw(ps, pg, os.path.join(OUT_DIR, "THESIS_A5_top25_heatmap.png"),
     figsize=(8.0, 5.9), row_fs=5.0, col_fs=3.4, legend_fs=5.0, dot_size=1.0, anno_bar_h=0.5)

# ── 2) Article: all sig terms ──────────────────────────────────────────────────
ps2, pg2 = build_matrices()
order2 = cluster_order(ps2); ps2, pg2 = ps2.loc[order2], pg2.loc[order2]
nrow2, ncol2 = ps2.shape
art_w = max(14, 2.0 + 0.22 * ncol2); art_h = max(9, 1.4 + 0.20 * nrow2)
print("Article (all sig terms):")
draw(ps2, pg2, os.path.join(OUT_DIR, "ARTICLE_heatmap_no_dendro.png"),
     figsize=(art_w, art_h), row_fs=5.5, col_fs=4.6, legend_fs=8.0, dot_size=1.6, anno_bar_h=0.6)

print("Done. Outputs in:", OUT_DIR)
