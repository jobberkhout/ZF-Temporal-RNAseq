#!/usr/bin/env python3
"""
Control GSVA heatmaps (NTC focus set) — same style as the big GSEA heatmap,
with the timepoint colour bar using the PCA palette.
Rows = GO terms (clustered, no dendrogram); columns = timepoints T1-T5.
Outputs: combined DMSO|E3 heatmaps, in Output/Figures/.
"""
import pandas as pd, numpy as np, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

from pathlib import Path
ROOT = str(Path(__file__).resolve().parents[2])   # repository root (this file lives in Scripts/figures/)
OUT_DIR = os.path.join(ROOT, "Output", "Figures")
DMSO_CSV = os.path.join(ROOT, "Output/GSVA_NTC/Group/DMSO/GSVA_scores_timepoint_means.csv")
E3_CSV   = os.path.join(ROOT, "Output/GSVA_NTC/Group/E3/GSVA_scores_timepoint_means.csv")
TPS = ["T1", "T2", "T3", "T4", "T5"]
# Selection mode: "range" = all terms changing across stages (>= MIN_RANGE, ~141); "core" = NTC-specific + top non-core.
MIN_RANGE       = 0.05   # ALL set: keep terms whose ssGSEA changes >= this across T1-T5 (either vehicle) -> ~141
CORE_MIN_RANGE  = 0.02   # NTC-core figure: drop only dead-flat NTC terms -> ~51
OTHER_N         = 50     # Other figure: this many most-dynamic non-NTC terms
NTC_CORE_RX = (r"neural tube|neurulation|neural plate|neural fold|neuropore|neural rod|neural crest|floor plate|"
               r"roof plate|hinge|notochord|convergent extension|planar cell polarity|\bPCP\b|apical constriction|"
               r"midline|dorsoventral|dorsal.?ventral|somitogen|\bsomite|axis elongation|body axis|gastrulation|"
               r"epiboly|mesenchymal transition|delamination|hedgehog|smoothened|\bshh\b|\bbmp\b|\bwnt\b|"
               r"retinoic acid|retinoid|noggin|chordin")
# PCA palette (Create_control_PCAs_poster.R) mapped T1..T5
TP_COLORS = {"T1": "#C1A3CD", "T2": "#604A7B", "T3": "#B0E298", "T4": "#FC98A9", "T5": "#8596CC"}
MAXLEN = 42
DPI = 300
os.makedirs(OUT_DIR, exist_ok=True)

def shorten(s, maxlen=MAXLEN):
    s = str(s).strip()
    for a, b in [("positive regulation of", "pos. reg."), ("negative regulation of", "neg. reg."),
                 ("regulation of", "reg."), ("involved in", "in"), ("via plasma membrane", "via PM"),
                 ("cell adhesion molecules", "CAMs"), ("signaling pathway", "signaling"),
                 ("morphogenesis", "morph."), ("development", "dev."), ("differentiation", "diff."),
                 ("organization", "org."), ("specification", "spec."), ("extracellular matrix", "ECM"),
                 ("epithelial to mesenchymal transition", "EMT"), ("biosynthetic process", "biosynth.")]:
        s = s.replace(a, b)
    return s if len(s) <= maxlen else s[:maxlen-1].rstrip() + "…"

dmso = pd.read_csv(DMSO_CSV, index_col=0)
e3   = pd.read_csv(E3_CSV, index_col=0)
import re as _re
common0 = [t for t in dmso.index if t in set(e3.index)]
Dd0 = dmso.loc[common0, TPS].values
De0 = e3.loc[common0, TPS].values
rngD = Dd0.max(1) - Dd0.min(1); rngE = De0.max(1) - De0.min(1)
interest = np.maximum(rngD, rngE)
is_core = np.array([bool(_re.compile(NTC_CORE_RX, _re.I).search(c)) for c in common0])

# --- subsets ---
all_idx   = np.where(interest >= MIN_RANGE)[0]                       # full "changes across stages" set (~141)
core_idx  = np.where(is_core & (interest >= CORE_MIN_RANGE))[0]      # NTC-specific terms (~51)  -> A5 figure 1
noncore   = np.where(~is_core)[0]
other_idx = noncore[np.argsort(interest[noncore])[::-1][:OTHER_N]]   # most-dynamic non-NTC terms -> A5 figure 2
# shared colour scale so all figures are directly comparable
_u = np.unique(np.concatenate([all_idx, core_idx, other_idx]))
vmax_shared = max(float(np.nanmax(np.abs(np.column_stack([Dd0, De0])[_u]))), 0.2)

def draw(mat, col_meta, labels, vmax, out_path, figsize, row_fs, col_fs, legend_fs, vehicle_groups=True, anno_h=0.6):
    nrow, ncol = mat.shape
    fig, ax = plt.subplots(figsize=figsize, dpi=DPI)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r",
                   norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), interpolation="nearest")
    # timepoint colour bar (PCA colours) flush on top
    bar_top = -0.5 - anno_h
    for j, (veh, tp) in enumerate(col_meta):
        ax.add_patch(mpatches.Rectangle((j-0.5, bar_top), 1, anno_h, facecolor=TP_COLORS[tp], edgecolor="none", clip_on=False))
    head = max(2.5, 0.035 * nrow) if vehicle_groups else 0.0   # headroom for DMSO/E3 labels
    ax.set_xlim(-0.5, ncol-0.5)
    ax.set_ylim(nrow-0.5, bar_top - head)
    # vehicle separators + labels
    if vehicle_groups:
        prev = None; starts = {}
        for j, (veh, tp) in enumerate(col_meta):
            if veh != prev and prev is not None:
                ax.axvline(j-0.5, color="black", linewidth=1.4)
            starts.setdefault(veh, []).append(j); prev = veh
        for veh, js in starts.items():
            ax.text(np.mean(js), bar_top - 0.55, veh, ha="center", va="bottom",
                    fontsize=legend_fs+2, fontweight="bold", clip_on=False)
    ax.set_xticks(range(ncol)); ax.set_xticklabels([tp for _, tp in col_meta], fontsize=col_fs)
    ax.set_yticks(range(nrow)); ax.set_yticklabels(labels, fontsize=row_fs)
    ax.tick_params(axis="both", which="both", length=0)

    # colorbar + timepoint legend (compact, right)
    cax = ax.inset_axes([1.02, 0.0, 0.02, 0.5])
    cbar = fig.colorbar(im, cax=cax); cbar.set_label("ssGSEA score", fontsize=legend_fs)
    cbar.ax.tick_params(labelsize=legend_fs-1)
    patches = [mpatches.Patch(facecolor=TP_COLORS[t], label=t) for t in TPS]
    leg = ax.legend(handles=patches, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                    fontsize=legend_fs, framealpha=0, edgecolor="none",
                    title="Timepoint", title_fontsize=legend_fs+1,
                    handlelength=1.0, handletextpad=0.4, labelspacing=0.25, borderpad=0.2)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white", pad_inches=0.02, bbox_extra_artists=(leg,))
    if out_path.endswith(".png"):
        fig.savefig(out_path.replace(".png", ".pdf"), dpi=DPI, bbox_inches="tight", facecolor="white", pad_inches=0.02, bbox_extra_artists=(leg,))
    plt.close(fig); print("  saved:", os.path.basename(out_path), f"({nrow}x{ncol})")

col_meta_c = [("DMSO", t) for t in TPS] + [("E3", t) for t in TPS]

def render_set(idx, out_name):
    cc = [common0[i] for i in idx]; dd = Dd0[idx]; ee = De0[idx]
    M = np.column_stack([dd, ee])
    order = leaves_list(linkage(np.nan_to_num(pdist(M, metric="euclidean")), method="ward"))
    cc = [cc[i] for i in order]; dd = dd[order]; ee = ee[order]
    labels = [shorten(t) for t in cc]; n = len(cc)
    rfs = 6.5 if n <= 60 else (5.0 if n <= 80 else 4.2)
    Hc = max(4.5, 0.105 * n); Wc = max(7.0, Hc / 1.45)   # wide aspect -> fits at \textwidth
    print(f"{out_name}: {n} terms")
    draw(np.column_stack([dd, ee]), col_meta_c, labels, vmax_shared,
         os.path.join(OUT_DIR, out_name + ".png"),
         figsize=(Wc, Hc), row_fs=rfs, col_fs=7, legend_fs=8, vehicle_groups=True)

render_set(all_idx,   "GSVA_controls_DMSO_E3_heatmap")          # full set (the big one, ~141)
render_set(core_idx,  "GSVA_controls_NTCcore_DMSO_E3_heatmap")  # A5 figure 1: NTC-specific
render_set(other_idx, "GSVA_controls_other_DMSO_E3_heatmap")    # A5 figure 2: most-dynamic non-NTC
print("Done.")
