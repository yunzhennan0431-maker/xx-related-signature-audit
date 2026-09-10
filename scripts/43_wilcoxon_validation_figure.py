# -*- coding: utf-8 -*-
"""
生成图5:z-score / 未下采样Wilcoxon / 下采样Wilcoxon 三方对比图。
左图:108模型显著率(%)对比;右图:H1a(TME细胞集中)与H1b(肿瘤细胞自身集中)计数对比。
数据来自41_full_wilcoxon_robustness_check.py与42_full_wilcoxon_downsampled_check.py的输出,
用于3.3节独立算法验证结果的可视化呈现。
"""
import matplotlib.pyplot as plt
import matplotlib
import os

matplotlib.rcParams["font.sans-serif"] = ["Arial Unicode MS", "PingFang SC", "Heiti SC", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

OUT_DIR = "analysis_output/figures/"
os.makedirs(OUT_DIR, exist_ok=True)

METHODS = ["z-score\n(现有结论)", "Wilcoxon\n(未下采样)", "Wilcoxon\n(下采样均衡)"]
PCT_SIG = [23.1, 58.3, 37.0]
N_SIG = [25, 63, 40]
H1A = [16, 19, 18]
H1B = [6, 29, 16]
COLORS = ["#2c6fbb", "#c0392b", "#e67e22"]


def fig_wilcoxon_validation():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.6))

    bars = ax1.bar(METHODS, PCT_SIG, color=COLORS, width=0.55)
    for b, pct, n in zip(bars, PCT_SIG, N_SIG):
        ax1.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.2,
                  f"{n}/108\n({pct:.1f}%)", ha="center", fontsize=9.5)
    ax1.set_ylabel("FDR<0.05显著模型占比 (%)", fontsize=10.5)
    ax1.set_ylim(0, 70)
    ax1.set_title("(A) 显著率对比", fontsize=11)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    x = range(len(METHODS))
    w = 0.32
    b1 = ax2.bar([i - w / 2 for i in x], H1A, width=w, color="#2c6fbb", label="H1a(TME细胞集中)")
    b2 = ax2.bar([i + w / 2 for i in x], H1B, width=w, color="#c0392b", label="H1b(肿瘤细胞自身集中)")
    for b, v in zip(b1, H1A):
        ax2.text(b.get_x() + b.get_width() / 2, v + 0.6, str(v), ha="center", fontsize=9.5)
    for b, v in zip(b2, H1B):
        ax2.text(b.get_x() + b.get_width() / 2, v + 0.6, str(v), ha="center", fontsize=9.5)
    ax2.set_xticks(list(x))
    ax2.set_xticklabels(METHODS, fontsize=9.5)
    ax2.set_ylabel("显著模型数", fontsize=10.5)
    ax2.set_ylim(0, 34)
    ax2.set_title("(B) H1a : H1b 分布对比", fontsize=11)
    ax2.legend(fontsize=8.8, loc="upper left", frameon=False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.suptitle("独立算法(Wilcoxon)对归因结果稳健性的验证", fontsize=12.5, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR + "fig5_wilcoxon_validation.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR + "fig5_wilcoxon_validation.pdf", bbox_inches="tight")
    plt.close(fig)
    print("saved fig5_wilcoxon_validation.png/.pdf")


if __name__ == "__main__":
    fig_wilcoxon_validation()
