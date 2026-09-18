# -*- coding: utf-8 -*-
"""
生成PRISMA风格的文献检索/筛选流程图。

重要说明:本研究不是一项正式的系统综述/meta分析(不做疗效合成),而是一项方法学审计,
候选文献的优先级排序是基于规则化打分脚本而非人工逐篇双人独立筛选;因此这张图是"参照
PRISMA 2020流程图结构改编"而非严格意义上的PRISMA流程图,图题与脚注均明确注明这一点,
避免读者误认为本研究经过了正式的PRISMA双人独立筛选流程。

输出: analysis_output/figures/fig6_prisma_flow.png/.pdf
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib

matplotlib.rcParams["font.sans-serif"] = ["Arial Unicode MS", "PingFang SC", "Heiti SC", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

import os
OUT_DIR = "analysis_output/figures/"
os.makedirs(OUT_DIR, exist_ok=True)

fig, ax = plt.subplots(figsize=(9.4, 12))
ax.set_xlim(-0.5, 10.5)
ax.set_ylim(0, 34)
ax.axis("off")


def box(cx, cy, w, h, text, fc="#eef3fb", ec="#2c6fbb", fontsize=9):
    rect = mpatches.FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.08,rounding_size=0.15",
        linewidth=1.3, edgecolor=ec, facecolor=fc,
    )
    ax.add_patch(rect)
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fontsize, linespacing=1.35)


def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color="#444", lw=1.2))


# Identification
box(5, 32.5, 8.6, 2.4,
    "标识(Identification)\nPubMed检索,23个初始种子主题词 + 25个后续扩充主题词\n共48个过程/通路关键词,3轮检索\n原始命中(未去重):4,022 + 1,243 = 5,265条", fontsize=9)
arrow(5, 31.3, 5, 30.1)

box(5, 29.3, 8.6, 2.0,
    "经idconv批量PMID→PMCID转换、PMC Cloud Service可得性核查后\n三轮共获得1,662篇完整文章包(正文+补充材料),构成候选文献池", fontsize=9)
arrow(5, 28.3, 5, 27.1)

# Screening
box(5, 26.3, 8.6, 1.6,
    "筛选(Screening)\n规则化脚本对候选表格打分排序(正向/负向关键词规则)\n1,662篇候选按打分优先级排队", fontsize=9)
arrow(5, 25.5, 5, 24.3)

arrow(3.4, 23.5, 1.9, 22.3)
arrow(6.6, 23.5, 8.1, 22.3)
box(5, 23.8, 3.2, 0.9, "决定人工复核", fc="none", ec="none", fontsize=8.5)

box(1.9, 21.3, 3.6, 2.2,
    "已人工复核\n445篇 (26.8%)\n[单一研究者逐篇阅读判定]", fc="#eef3fb", fontsize=9)
box(8.1, 21.3, 3.6, 2.2,
    "未人工复核\n1,217篇 (73.2%)\n主动决定暂不推进\n(成本收益评估,详见5节)", fc="#fbeeee", ec="#c0392b", fontsize=8.5)

arrow(1.9, 20.2, 1.9, 19.0)
box(1.9, 18.0, 4.0, 2.4,
    "445篇复核结果\n初判命中(hit):131篇 (=445-314)\n判定为miss/未分类:314篇", fontsize=8.7)

arrow(1.9, 16.8, 1.9, 15.6)
box(1.9, 14.6, 4.0, 2.0,
    "第三轮验证(独立人类复核者\n对60篇抽样重判,kappa=0.106)\n发现miss判定流程系统性缺陷", fc="#fff8e1", ec="#b8860b", fontsize=8.5)

arrow(1.9, 13.6, 1.9, 12.4)
box(1.9, 11.4, 4.0, 2.4,
    "系统性复查全部314篇miss/未分类候选\n(先扫全文,再下载补充材料扫描)\n187篇(59.6%)确认可提取最终模型", fontsize=8.5)

arrow(1.85, 10.2, 1.85, 9.0)
arrow(3.9, 9.9, 6.3, 8.7)

box(1.85, 8.0, 3.0, 2.0, "127篇(40.4%)\n未解决\n(系数仅图片/分子类型\n超范围/确系候选池)", fc="#fbeeee", ec="#c0392b", fontsize=7.6)

# Included
box(6.5, 7.5, 5.8, 2.6,
    "纳入(Included)\n初判108个模型 + 复查新增173个模型\n= 281个独立签名模型\n(34+过程主题,21癌种,21个单细胞参考数据集)", fc="#eafaf1", ec="#1e8449", fontsize=9.0)

arrow(6.5, 6.2, 6.5, 5.0)
box(6.5, 4.0, 5.8, 2.2,
    "定量分析集: 281个模型\n(2个模型因目标癌种本地无参考单细胞数据未纳入定量归因,\n仅作方法学证据留痕,不在281之列)", fontsize=8.3)

fig.suptitle("图6. 文献检索与候选池筛选流程(参照PRISMA 2020结构改编)", fontsize=12, y=0.985)
fig.text(0.5, 0.008,
         "注:本研究为方法学审计而非正式系统综述,候选优先级基于规则化打分而非双人独立筛选;\n"
         "此图为按PRISMA 2020流程图结构改编的示意,非严格意义上的PRISMA流程图。",
         ha="center", fontsize=7.5, color="#555")

fig.tight_layout(rect=[0, 0.02, 1, 0.97])
fig.savefig(OUT_DIR + "fig6_prisma_flow.pdf")
fig.savefig(OUT_DIR + "fig6_prisma_flow.png", dpi=200)
plt.close(fig)
print("已生成:", OUT_DIR + "fig6_prisma_flow.png/.pdf")
