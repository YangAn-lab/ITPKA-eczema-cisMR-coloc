fig = plt.figure(figsize=(13.5, 5.4))
gs = GridSpec(1, 2, width_ratios=[1.05, 1.35], wspace=0.30)

ax = fig.add_subplot(gs[0])
ys = list(range(len(rows), 0, -1))
for (lab, fc, p, c), y in zip(rows, ys):
    ax.scatter(fc, y, s=64, color=c, zorder=3)
    ax.annotate(f"p = {p:.3f}", (0.265, y), va='center', ha='left',
                fontsize=8.5, color='#333333')  # fixed label column, clear of zero line
ax.axvline(0, color='black', lw=0.9, ls='--', zorder=1, ymin=0.0, ymax=0.88)
ax.set_yticks(ys)
ax.set_yticklabels([r[0] for r in rows], fontsize=8.8)
ax.set_xlim(-0.42, 0.42)
ax.set_xticks([-0.3, -0.15, 0, 0.15, 0.3])
ax.set_xlabel("ITPKA log2 fold-change", fontsize=10)
ax.set_ylim(0.3, 6.4)
ax.text(0.0, 6.08, "GSE121212 paired Wilcoxon cross-check: p = 0.471 (n.s.)",
        ha='center', fontsize=8.2, color=GREY, style='italic')
ax.spines[['top', 'right']].set_visible(False)
ax.set_title("A", loc='left', fontsize=14, fontweight='bold', x=-0.28)

axb = fig.add_subplot(gs[1])
axb.axis('off')
tbl = axb.table(cellText=[["Evidence layer", "Resource", "Result"]] + tbl_rows,
                colWidths=[0.32, 0.34, 0.46], loc='center', cellLoc='left')
tbl.auto_set_font_size(False)
tbl.set_fontsize(8.8)
tbl.scale(1, 1.75)
for j in range(3):
    c = tbl[0, j]
    c.set_facecolor('#ECE9E2'); c.set_text_props(fontweight='bold')
for i in range(1, 7):
    for j in range(3):
        tbl[i, j].set_facecolor('#FAF9F3' if i % 2 else 'white')
axb.set_title("B", loc='left', fontsize=14, fontweight='bold', x=-0.06)

fig.savefig('/mnt/results/P2-02_fig4_组织定位证据汇总.png', dpi=300, bbox_inches='tight')
fig.savefig('/mnt/results/P2-02_fig4_组织定位证据汇总.svg', bbox_inches='tight')
plt.show()
print("saved v3")