"""Figure 5 v5 (rev2): mechanistic hypothesis schematic, programmatic matplotlib drawing.
Fixes from media check: box-2 text overflow (3-line layout), bottom note realigned under
main column, wider boxes, shortened labels.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

fig, ax = plt.subplots(figsize=(10.2, 7.9), dpi=300)
ax.set_xlim(0, 10.2)
ax.set_ylim(0, 7.9)
ax.axis("off")

C_MAIN = "#DCE9F7"
C_EDGE = "#1F4E79"
C_OUT = "#FDEBD3"
C_OEDGE = "#B26A00"
C_GREY = "#EFEDEA"
C_GEDGE = "#7A7A7A"


def box(x, y, w, h, text, fc, ec, fs=9.2, lw=1.6):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.06,rounding_size=0.12",
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs,
            color="#111111", zorder=3, linespacing=1.45)


def arrow(x1, y1, x2, y2, color=C_EDGE, lw=2.0, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=16, lw=lw, color=color,
                                 linestyle=ls, zorder=1))


X = 3.55
W = 5.35

box(X, 7.02, W, 0.78,
    "rs11635906-G allele (15q15.1)\nITPKA cis-eQTL index variant (EAF 0.26–0.30)",
    C_MAIN, C_EDGE, fs=9.6)
box(X, 5.68, W, 1.06,
    "Up-regulation of blood ITPKA in circulating T-lineage cells\n"
    "(eQTLGen β=+0.150, p=7.4×10⁻²⁹;\n"
    "colocalization with eczema GWAS PP.H4=0.962)",
    C_MAIN, C_EDGE, fs=9.0)
box(X, 4.24, W, 1.02,
    "ITPKA converts Ins(1,4,5)P₃ to Ins(1,3,4,5)P₄,\n"
    "modulating T-cell-receptor–evoked Ca²⁺/NFAT signalling\n"
    "(IP3-kinase dosage sensitivity; direction not committed)",
    C_MAIN, C_EDGE, fs=9.0)
box(X, 2.94, W, 0.82,
    "Type-2-skewed (Tc2/Th2) activation of\nantigen-experienced memory T cells",
    C_MAIN, C_EDGE, fs=9.2)
box(X, 1.74, W, 0.72,
    "Skin homing of antigen-experienced memory T cells",
    C_MAIN, C_EDGE, fs=9.2)
box(X, 0.74, W, 0.78,
    "Cutaneous type-2 inflammation — increased dermatitis/eczema risk\n"
    "(OR 1.27, 95% CI 1.17–1.38, delta-method p=6.9×10⁻⁷)",
    C_OUT, C_OEDGE, fs=9.2)

for y1, y2 in [(6.63, 6.21), (5.15, 4.75), (3.73, 3.35), (2.53, 2.10), (1.38, 1.13)]:
    arrow(X, y1, X, y2)

# grey evidence box (right)
gx, gw, gh = 8.55, 3.08, 3.75
ax.add_patch(FancyBboxPatch((gx - gw / 2, 3.85 - gh / 2), gw, gh,
                            boxstyle="round,pad=0.06,rounding_size=0.12",
                            fc=C_GREY, ec=C_GEDGE, lw=1.4, zorder=2))
ax.text(gx, 5.28, "Uniformly negative\nskin-intrinsic evidence",
        ha="center", va="center", fontsize=8.8, fontweight="bold", color="#333333", linespacing=1.35)
ax.text(gx, 3.55,
        "• no skin cis-eQTL (GTEx v8/v10,\n   two skin sites)\n"
        "• no lesional-skin differential\n   expression (5/5 contrasts)\n"
        "• no detectable skin protein (HPA)\n"
        "• ITPKA⁺ cells rare in skin (≤0.66%)\n\n"
        "Consistent with the effect entering\nskin via migrating T cells rather than\nacting in skin-resident cells",
        ha="center", va="center", fontsize=8.2, color="#333333", linespacing=1.4)
arrow(gx - gw / 2 - 0.05, 2.75, X + W / 2 + 0.08, 1.74, color=C_GEDGE, lw=1.6, ls="--")

ax.text(X, 0.06,
        "Hypothesis with falsifiable predictions (Section 4.3) — not a demonstrated mechanism.",
        ha="center", va="bottom", fontsize=8.2, style="italic", color="#444444")

fig.savefig("/workspace/P2-02_Figure5_mechanism_hypothesis_v6.png",
            bbox_inches="tight", facecolor="white")
fig.savefig("/workspace/P2-02_Figure5_mechanism_hypothesis_v6.svg",
            bbox_inches="tight", facecolor="white")
print("saved rev2")
