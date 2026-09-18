# -*- coding: utf-8 -*-
"""English-labeled versions of the 4 manuscript figures, for the English review package."""
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

OUT_DIR = "figures_en/"
os.makedirs(OUT_DIR, exist_ok=True)


def draw_box(ax, xy, w, h, text, fc="#eaf2fb", ec="#2c6fbb", fontsize=8.6):
    box = FancyBboxPatch((xy[0] - w / 2, xy[1] - h / 2), w, h,
                          boxstyle="round,pad=0.02,rounding_size=0.04",
                          linewidth=1.3, edgecolor=ec, facecolor=fc, zorder=2)
    ax.add_patch(box)
    ax.text(xy[0], xy[1], text, ha="center", va="center", fontsize=fontsize,
             zorder=3, linespacing=1.35)


def draw_arrow(ax, xy1, xy2, color="#555"):
    arr = FancyArrowPatch(xy1, xy2, arrowstyle="-|>", mutation_scale=14,
                           linewidth=1.3, color=color, zorder=1)
    ax.add_patch(arr)


def fig_workflow():
    fig, ax = plt.subplots(figsize=(14.2, 8.6))
    ax.set_xlim(0, 16.0)
    ax.set_ylim(2.2, 15)
    ax.axis("off")

    main_w, main_h = 8.0, 0.95
    qc_w, qc_h = 3.4, 1.05
    xc = 8.0

    stages = [
        (14.2, "Topic word-list construction (23+25=48 seed/expanded topics)\n+ systematic PubMed E-utilities search (4,022+1,243 hits)"),
        (12.9, "PMC Cloud Service programmatic download of candidate article packages\n(1,662 final candidates)"),
        (11.6, "Rule-based scanning + keyword scoring/ranking of tables\n(distinguishing candidate pool / intermediate product / final model coefficient table)"),
        (10.3, "Manual semantic review, one by one (11 batches, 445 reviewed)\n→ 108 independent models with complete gene panels, extractable coefficients"),
        (9.0, "Single-cell reference atlas matched by cancer type\n(21 datasets, TISCH2/GEO public resources)"),
        (7.7, "z-score cell-type attribution\n(log-normalize expression → mean by cell type → z-score across types → argmax)"),
        (6.4, "Exact binomial test (H0: equal-probability attribution across K types) + Benjamini-Hochberg FDR\n(all 108 models corrected under one natural family)"),
        (5.1, "Cross-topic meta-analysis: 3-category grouping (Malignant/Stromal_Vascular/\nImmune_Hematopoietic/Epithelial_mixed) + H1a/H1b split"),
        (3.8, "Cross-topic heterogeneity permutation test (20,000 reshuffles, 9 topic categories)\n+ dataset-coverage / rare-cell-type sensitivity analyses"),
    ]
    for y, text in stages:
        draw_box(ax, (xc, y), main_w, main_h, text)
    for i in range(len(stages) - 1):
        y1 = stages[i][0] - main_h / 2
        y2 = stages[i + 1][0] + main_h / 2
        draw_arrow(ax, (xc, y1), (xc, y2))

    qc_items = [
        (10.3, "Extraction-consistency sampling check\n(2 rounds, independent blind re-review;\nfound & corrected 2 genuine omissions)"),
        (6.4, "Systematic candidate-pool-tautology screening\n(red-flag regex + manual review;\nconfirmed 5 cases, 4 construction templates)"),
        (5.1, "Full-scale pairwise duplicate-publication screening\n(Jaccard similarity; found & corrected\n1 mislabeled cancer type)"),
    ]
    qc_x = xc + main_w / 2 + qc_w / 2 + 0.4
    for y, text in qc_items:
        draw_box(ax, (qc_x, y), qc_w, qc_h, text, fc="#fdecea", ec="#c0392b", fontsize=7.2)
        draw_arrow(ax, (qc_x - qc_w / 2, y), (xc + main_w / 2 + 0.05, y), color="#c0392b")

    va_x = xc - main_w / 2 - qc_w / 2 - 0.4
    draw_box(ax, (va_x, 7.7), qc_w, 1.3,
             "Validation anchor, layer 1: cross-comparison\nagainst 13 gold-standard models from the\nprior angiogenesis study, 119 gene×dataset\npairs 100% match",
             fc="#eafaf1", ec="#27ae60", fontsize=7.0)
    draw_arrow(ax, (va_x + qc_w / 2, 7.7), (xc - main_w / 2 - 0.05, 7.7), color="#27ae60")

    draw_box(ax, (va_x, 9.0), qc_w, 1.3,
             "Validation anchor, supplement: author-\nreported single-cell observations in 7 models,\n16 gene-level comparisons, 9 (56%) exact matches",
             fc="#eafaf1", ec="#27ae60", fontsize=7.0)
    draw_arrow(ax, (va_x + qc_w / 2, 9.0), (xc - main_w / 2 - 0.05, 9.0), color="#27ae60")

    ax.text(xc, 14.85, "Figure 1. Study workflow diagram", ha="center", fontsize=13, fontweight="bold")
    ax.text(xc, 2.9,
            "Main pipeline (blue, top to bottom) links all method steps in §2.1–2.5;\n"
            "red boxes are quality-control procedures running throughout; green boxes are validation-anchor comparisons independent of the main pipeline.",
            ha="center", fontsize=8.0, color="#333")

    fig.tight_layout()
    fig.savefig(OUT_DIR + "fig0_study_workflow.pdf")
    fig.savefig(OUT_DIR + "fig0_study_workflow.png", dpi=200)
    plt.close(fig)
    print("wrote fig0:", OUT_DIR + "fig0_study_workflow.png/.pdf")


def fig_h1a_h1b_distribution():
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))

    all_cats = ["Immune/\nhematopoietic", "Stromal/\nvascular", "Malignant", "Mixed\nepithelial"]
    all_counts = [36, 34, 30, 8]
    colors_all = ["#e67e22", "#2c6fbb", "#c0392b", "#8e44ad"]
    ax = axes[0]
    bars = ax.bar(all_cats, all_counts, color=colors_all, width=0.6)
    for b, c in zip(bars, all_counts):
        ax.annotate(str(c), (b.get_x() + b.get_width() / 2, b.get_height()),
                    textcoords="offset points", xytext=(0, 4), ha="center", fontsize=10)
    ax.set_ylabel("Number of models")
    ax.set_title("Baseline distribution,\nall 108 models (regardless of significance)", fontsize=9.6)
    ax.set_ylim(0, 42)

    h1_cats = ["H1a", "H1b", "Other"]
    h1_counts = [16, 6, 1]
    colors_h1 = ["#2c6fbb", "#c0392b", "#8e44ad"]
    ax2 = axes[1]
    bars2 = ax2.bar(h1_cats, h1_counts, color=colors_h1, width=0.55)
    for b, c in zip(bars2, h1_counts):
        ax2.annotate(str(c), (b.get_x() + b.get_width() / 2, b.get_height()),
                     textcoords="offset points", xytext=(0, 4), ha="center", fontsize=10)
    ax2.set_ylabel("Number of models")
    ax2.set_title("H1a/H1b split, 23 genuinely significant\nmodels after removing 2 tautology cases\n(full 108-model scale)", fontsize=9.2)
    ax2.set_ylim(0, 20)

    h1_counts_exp = [6, 5, 1]
    ax3 = axes[2]
    bars3 = ax3.bar(h1_cats, h1_counts_exp, color=colors_h1, width=0.55)
    for b, c in zip(bars3, h1_counts_exp):
        ax3.annotate(str(c), (b.get_x() + b.get_width() / 2, b.get_height()),
                     textcoords="offset points", xytext=(0, 4), ha="center", fontsize=10)
    ax3.set_ylabel("Number of models")
    ax3.set_title("After excluding 4 datasets with\n1/4–2/4 lineage coverage (83 models):\nH1a/H1b has narrowed to near 1:1", fontsize=9.2)
    ax3.set_ylim(0, 20)

    legend_handles = [
        mpatches.Patch(color="#2c6fbb", label="H1a: TME-cell concentration (genuine disconnect evidence)"),
        mpatches.Patch(color="#c0392b", label="H1b: Malignant-cell-self concentration (interpretation uncertain)"),
        mpatches.Patch(color="#8e44ad", label="Other: mixed-epithelial category (not entered into the binary split)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=1, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.14))

    fig.suptitle("Figure 4. Mode-cell-type broad-category distribution: full scale vs. after excluding\nlineage-coverage-incomplete datasets", fontsize=12, y=1.08)
    fig.tight_layout()
    fig.savefig(OUT_DIR + "fig3_h1a_h1b_distribution.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR + "fig3_h1a_h1b_distribution.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote fig3:", OUT_DIR + "fig3_h1a_h1b_distribution.png/.pdf")


def fig_trajectories():
    n_models_a = [19, 25, 27, 37, 48, 61, 76, 77, 93, 102, 108]
    rates_a = [47.4, 40.0, 33.3, 24.3, 20.8, 13.1, 19.7, 20.8, 20.4, 19.6, 23.1]

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.plot(n_models_a, rates_a, marker="o", color="#2c6fbb", linewidth=1.8, markersize=6)
    # Per-point label offsets (points, in xytext units) chosen to keep labels clear of the
    # connecting lines; the two close-together checkpoints at x=76/77 are split above/below.
    label_offsets = {
        25: (-4, 16),
        27: (10, 12),
        76: (-6, -16),
        77: (10, 14),
    }
    for x, y in zip(n_models_a, rates_a):
        dx, dy = label_offsets.get(x, (0, 9))
        ax.annotate(f"{y:.1f}%", (x, y), textcoords="offset points", xytext=(dx, dy),
                    ha="center", fontsize=8, color="#333")
    ax.set_xlabel("Cumulative number of models included (N)")
    ax.set_ylabel("Proportion significant at FDR<0.05 (%)")
    ax.set_title("Figure 2. Full trajectory of the significant-concentration proportion\nacross 11 sample-size checkpoints")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 57)
    fig.tight_layout()
    fig.savefig(OUT_DIR + "fig_significant_rate_trajectory.pdf")
    fig.savefig(OUT_DIR + "fig_significant_rate_trajectory.png", dpi=200)
    plt.close(fig)
    print("wrote trajectory A:", OUT_DIR + "fig_significant_rate_trajectory.png/.pdf")

    n_models_b = [77, 93, 102, 108]
    pvals_b = [0.20, 0.0787, 0.1547, 0.3220]

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.plot(n_models_b, pvals_b, marker="s", color="#c0392b", linewidth=1.8, markersize=7)
    ax.axhline(0.05, color="gray", linestyle="--", linewidth=1, label="Significance threshold p=0.05")
    # Per-point label offsets: 102 is nudged up-left so the ascending segment to 108 doesn't
    # cross it, and 108 is nudged left so its label doesn't clip the right/top plot edge.
    label_offsets = {
        102: (-20, 8),
        108: (-18, 10),
    }
    for x, y in zip(n_models_b, pvals_b):
        dx, dy = label_offsets.get(x, (0, 9))
        ax.annotate(f"{y:.3f}", (x, y), textcoords="offset points", xytext=(dx, dy),
                    ha="center", fontsize=8, color="#333")
    ax.set_xlabel("Cumulative number of models included (N)")
    ax.set_ylabel("Permutation-test p-value (cross-topic heterogeneity)")
    ax.set_title("Figure 3. Non-monotonic trajectory of the cross-topic heterogeneity\npermutation-test p-value across 4 checkpoints")
    ax.set_ylim(0, 0.43)
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR + "fig_heterogeneity_pvalue_trajectory.pdf")
    fig.savefig(OUT_DIR + "fig_heterogeneity_pvalue_trajectory.png", dpi=200)
    plt.close(fig)
    print("wrote trajectory B:", OUT_DIR + "fig_heterogeneity_pvalue_trajectory.png/.pdf")


if __name__ == "__main__":
    fig_workflow()
    fig_h1a_h1b_distribution()
    fig_trajectories()
