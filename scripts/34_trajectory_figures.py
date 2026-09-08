# -*- coding: utf-8 -*-
"""
生成两张论文正文引用的趋势图:
图A. 显著集中比例随样本规模变化的11个checkpoint轨迹(3.5/4/5/6节反复用文字描述的非单调曲线)。
图B. 跨主题异质性置换检验p值随样本规模变化的4个checkpoint轨迹。
两图均保存为矢量PDF+PNG,供正式排版时选用;当前中文稿以Markdown形式呈现,先产出图片文件,
正文中以图注文字说明其存在并给出文件路径。
"""
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams["font.sans-serif"] = ["Arial Unicode MS", "PingFang SC", "Heiti SC", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

OUT_DIR = "analysis_output/figures/"
import os
os.makedirs(OUT_DIR, exist_ok=True)

# 图A: 显著集中比例轨迹(11个checkpoint)
n_models_a = [19, 25, 27, 37, 48, 61, 76, 77, 93, 102, 108]
rates_a = [47.4, 40.0, 33.3, 24.3, 20.8, 13.1, 19.7, 20.8, 20.4, 19.6, 23.1]

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(n_models_a, rates_a, marker="o", color="#2c6fbb", linewidth=1.8, markersize=6)
for x, y in zip(n_models_a, rates_a):
    ax.annotate(f"{y:.1f}%", (x, y), textcoords="offset points", xytext=(0, 8),
                ha="center", fontsize=8, color="#333")
ax.set_xlabel("累计纳入模型数(N)")
ax.set_ylabel("FDR<0.05 显著集中比例 (%)")
ax.set_title("图2. 显著集中比例随样本规模变化的完整轨迹(11个checkpoint)")
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 55)
fig.tight_layout()
fig.savefig(OUT_DIR + "fig_significant_rate_trajectory.pdf")
fig.savefig(OUT_DIR + "fig_significant_rate_trajectory.png", dpi=200)
plt.close(fig)
print("已生成图A:", OUT_DIR + "fig_significant_rate_trajectory.png/.pdf")

# 图B: 异质性置换检验p值轨迹(4个checkpoint)
n_models_b = [77, 93, 102, 108]
# 77模型取论文报告区间(0.19-0.21,两种主题归类粒度)的中点0.20作代表值;
# 93/102/108模型为脚本实测精确值(见25_heterogeneity_test.py历次运行记录)。
pvals_b = [0.20, 0.0787, 0.1547, 0.3220]

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(n_models_b, pvals_b, marker="s", color="#c0392b", linewidth=1.8, markersize=7)
ax.axhline(0.05, color="gray", linestyle="--", linewidth=1, label="显著性阈值 p=0.05")
for x, y in zip(n_models_b, pvals_b):
    ax.annotate(f"{y:.3f}", (x, y), textcoords="offset points", xytext=(0, 8),
                ha="center", fontsize=8, color="#333")
ax.set_xlabel("累计纳入模型数(N)")
ax.set_ylabel("置换检验 p 值(跨主题异质性)")
ax.set_title("图3. 跨主题异质性置换检验p值的非单调轨迹(4个checkpoint)")
ax.set_ylim(0, 0.4)
ax.legend(loc="upper left", fontsize=8)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(OUT_DIR + "fig_heterogeneity_pvalue_trajectory.pdf")
fig.savefig(OUT_DIR + "fig_heterogeneity_pvalue_trajectory.png", dpi=200)
plt.close(fig)
print("已生成图B:", OUT_DIR + "fig_heterogeneity_pvalue_trajectory.png/.pdf")
