# -*- coding: utf-8 -*-
"""
生成两张补充图:
图0(方法流程图):2.1-2.5节文字描述的完整流水线,从主题词表构建到跨主题meta分析,
                 并标注贯穿其中的三类质量控制环节(候选池同义重复核查、抽取一致性抽样检验、
                 重复发表排查)。
图3(H1a/H1b三大类分布对比图):左图为全部108模型的众数细胞类型三大类分布(基线构成),
                              右图为剔除2个同义重复案例后23个真正显著模型的H1a/H1b/其他
                              三分类分布,直观呈现"TME细胞集中多于肿瘤细胞自身集中"这一
                              核心发现,避免读者只能从大段文字里拼数字。
"""
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

matplotlib.rcParams["font.sans-serif"] = ["Arial Unicode MS", "PingFang SC", "Heiti SC", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

OUT_DIR = "analysis_output/figures/"
os.makedirs(OUT_DIR, exist_ok=True)


def draw_box(ax, xy, w, h, text, fc="#eaf2fb", ec="#2c6fbb", fontsize=9.5):
    box = FancyBboxPatch((xy[0] - w / 2, xy[1] - h / 2), w, h,
                          boxstyle="round,pad=0.02,rounding_size=0.04",
                          linewidth=1.3, edgecolor=ec, facecolor=fc, zorder=2)
    ax.add_patch(box)
    ax.text(xy[0], xy[1], text, ha="center", va="center", fontsize=fontsize,
             zorder=3, linespacing=1.4)


def draw_arrow(ax, xy1, xy2, color="#555"):
    arr = FancyArrowPatch(xy1, xy2, arrowstyle="-|>", mutation_scale=14,
                           linewidth=1.3, color=color, zorder=1)
    ax.add_patch(arr)


def fig_workflow():
    fig, ax = plt.subplots(figsize=(8.5, 8.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(2.2, 15)
    ax.axis("off")

    main_w, main_h = 7.6, 0.95
    qc_w, qc_h = 2.6, 0.85
    xc = 4.6

    stages = [
        (14.2, "主题词表构建(23+25=48个种子/扩充主题)\n+ PubMed E-utilities系统检索(4,022+1,243条命中)"),
        (12.9, "PMC Cloud Service程序化下载候选文章包\n(1,662篇最终候选池)"),
        (11.6, "规则化脚本扫描表格+关键词打分排序\n(区分候选池/中间产物/最终模型系数表)"),
        (10.3, "人工语义复核逐篇判断(11轮批次,445篇已复核)\n→ 108个基因面板完整、系数可提取的独立模型"),
        (9.0, "按癌种匹配单细胞参考图谱\n(21个数据集,TISCH2/GEO公开资源)"),
        (7.7, "z-score细胞类型归因\n(log标准化表达→按细胞类型求均值→跨类型z-score→取最大值)"),
        (6.4, "精确二项检验(H0:K类型等概率归因)+ Benjamini-Hochberg FDR校正\n(全部108模型同一自然家族统一校正)"),
        (5.1, "跨主题meta分析:三大类归类(Malignant/Stromal_Vascular/\nImmune_Hematopoietic/Epithelial_mixed)+ H1a/H1b分野"),
        (3.8, "跨主题异质性置换检验(2万次重排,9个主题类别)\n+ 数据集覆盖/稀有细胞类型敏感性分析"),
    ]
    for y, text in stages:
        draw_box(ax, (xc, y), main_w, main_h, text)
    for i in range(len(stages) - 1):
        y1 = stages[i][0] - main_h / 2
        y2 = stages[i + 1][0] + main_h / 2
        draw_arrow(ax, (xc, y1), (xc, y2))

    qc_items = [
        (10.3, "抽取一致性抽样检验\n(2轮独立盲复核,\n发现并纠正2处真实遗漏)"),
        (6.4, "候选池同义重复系统性核查\n(红旗正则+人工核查,\n确认4例、3种构造模板)"),
        (5.1, "全量模型两两重复发表排查\n(Jaccard相似度,\n发现并改正1处误标癌种)"),
    ]
    qc_x = xc + main_w / 2 + qc_w / 2 + 0.35
    for y, text in qc_items:
        draw_box(ax, (qc_x, y), qc_w, qc_h, text, fc="#fdecea", ec="#c0392b", fontsize=8.3)
        draw_arrow(ax, (qc_x - qc_w / 2, y), (xc + main_w / 2 + 0.05, y), color="#c0392b")

    va_x = xc - main_w / 2 - qc_w / 2 - 0.35
    draw_box(ax, (va_x, 7.7), qc_w, 1.15,
             "验证锚点第一层:\n与前序血管生成课题13个\n模型金标准比对,\n119组基因×数据集100%一致",
             fc="#eafaf1", ec="#27ae60", fontsize=8.0)
    draw_arrow(ax, (va_x + qc_w / 2, 7.7), (xc - main_w / 2 - 0.05, 7.7), color="#27ae60")

    draw_box(ax, (va_x, 9.0), qc_w, 1.15,
             "验证锚点补充:\n7个模型作者自述单细胞\n观察,16组基因级比对\n9组(56%)精确匹配",
             fc="#eafaf1", ec="#27ae60", fontsize=8.0)
    draw_arrow(ax, (va_x + qc_w / 2, 9.0), (xc - main_w / 2 - 0.05, 9.0), color="#27ae60")

    ax.text(xc, 14.85, "图1. 研究方法流程图", ha="center", fontsize=13, fontweight="bold")
    ax.text(xc, 2.9,
            "主流程(蓝色,自上而下)串联2.1–2.5节全部方法步骤;\n"
            "红色为贯穿全流程的质量控制环节,绿色为独立于主流程的验证锚点比对。",
            ha="center", fontsize=8.8, color="#333")

    fig.tight_layout()
    fig.savefig(OUT_DIR + "fig0_study_workflow.pdf")
    fig.savefig(OUT_DIR + "fig0_study_workflow.png", dpi=200)
    plt.close(fig)
    print("已生成图0:", OUT_DIR + "fig0_study_workflow.png/.pdf")


def fig_h1a_h1b_distribution():
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))

    all_cats = ["免疫/造血\n细胞", "基质血管\n细胞", "肿瘤/恶性\n细胞", "上皮细胞\n混合"]
    all_counts = [36, 34, 30, 8]
    colors_all = ["#e67e22", "#2c6fbb", "#c0392b", "#8e44ad"]
    ax = axes[0]
    bars = ax.bar(all_cats, all_counts, color=colors_all, width=0.6)
    for b, c in zip(bars, all_counts):
        ax.annotate(str(c), (b.get_x() + b.get_width() / 2, b.get_height()),
                    textcoords="offset points", xytext=(0, 4), ha="center", fontsize=10)
    ax.set_ylabel("模型数")
    ax.set_title("全部108个模型\n(不论显著与否)的基线分布", fontsize=10.5)
    ax.set_ylim(0, 42)

    h1_cats = ["H1a\nTME细胞集中\n(真实脱节证据)", "H1b\n肿瘤细胞自身集中\n(解读存疑)", "其他\n(上皮混合类别)"]
    h1_counts = [16, 6, 1]
    colors_h1 = ["#2c6fbb", "#c0392b", "#8e44ad"]
    ax2 = axes[1]
    bars2 = ax2.bar(h1_cats, h1_counts, color=colors_h1, width=0.55)
    for b, c in zip(bars2, h1_counts):
        ax2.annotate(str(c), (b.get_x() + b.get_width() / 2, b.get_height()),
                     textcoords="offset points", xytext=(0, 4), ha="center", fontsize=10)
    ax2.set_ylabel("模型数")
    ax2.set_title("剔除2例候选池同义重复后\n23个真正显著模型的H1a/H1b分野\n(全量108模型)", fontsize=10.5)
    ax2.set_ylim(0, 20)

    h1_counts_exp = [6, 5, 1]
    ax3 = axes[2]
    bars3 = ax3.bar(h1_cats, h1_counts_exp, color=colors_h1, width=0.55)
    for b, c in zip(bars3, h1_counts_exp):
        ax3.annotate(str(c), (b.get_x() + b.get_width() / 2, b.get_height()),
                     textcoords="offset points", xytext=(0, 4), ha="center", fontsize=10)
    ax3.set_ylabel("模型数")
    ax3.set_title("排除4个谱系覆盖度1/4~2/4\n的数据集后(83模型),\nH1a/H1b已收窄至接近1:1", fontsize=10.5)
    ax3.set_ylim(0, 20)

    fig.suptitle("图4. 众数细胞类型三大类分布:全部模型基线 vs H1a/H1b分野(全量 vs 排除覆盖不完整数据集后)",
                 fontsize=12.5, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR + "fig3_h1a_h1b_distribution.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR + "fig3_h1a_h1b_distribution.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("已生成图3:", OUT_DIR + "fig3_h1a_h1b_distribution.png/.pdf")


if __name__ == "__main__":
    fig_workflow()
    fig_h1a_h1b_distribution()
