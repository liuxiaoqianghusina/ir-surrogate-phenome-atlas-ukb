"""mr04：图 7 —— 孟德尔随机化。a 单变量 IVW（TG、BMI、WHRadjBMI）；b 多变量 MR（三者同时）。每 1 SD 遗传预测暴露的 OR（95% CI），FinnGen R12。"""
import pandas as pd, numpy as np, io, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
RES, FIG = "../results", "../figures"; MM = 1 / 25.4
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "axes.titlesize": 6.8, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6, "axes.linewidth": 0.5,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
def save(fig, name):
    fig.savefig(f"{FIG}/{name}.pdf", bbox_inches="tight"); buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=600, bbox_inches="tight", facecolor="white"); buf.seek(0)
    im = Image.open(buf).convert("RGB"); im.save(f"{FIG}/{name}.tiff", format="TIFF", compression="tiff_lzw", dpi=(600, 600)); im.resize((im.size[0] // 3, im.size[1] // 3)).save(f"{FIG}/{name}_preview.png"); plt.close(fig)
ORD = ["I9_CHD", "T2D", "M13_GOUT", "K11_CHOLELITH", "K11_ACUTPANC", "NAFLD", "N14_CHRONKIDNEYDIS", "I9_HEARTFAIL", "I9_AF", "I9_STR"]
LBL = {"I9_CHD": "Coronary heart disease", "T2D": "Type 2 diabetes", "M13_GOUT": "Gout", "K11_CHOLELITH": "Gallstones", "K11_ACUTPANC": "Acute pancreatitis",
       "NAFLD": "Fatty liver disease", "N14_CHRONKIDNEYDIS": "Chronic kidney disease", "I9_HEARTFAIL": "Heart failure", "I9_AF": "Atrial fibrillation", "I9_STR": "Stroke"}
EX = [("TG", "Triglycerides", "#B2182B"), ("BMI", "BMI", "#2166AC"), ("WHRadjBMI", "WHR adjusted for BMI", "#542788")]
U = pd.read_csv(f"{RES}/mr03_univariable.csv"); U = U[U.method == "IVW"]; MV = pd.read_csv(f"{RES}/mr03_mvmr.csv"); MV = MV[MV.model == "TG+BMI+WHRadjBMI"]
fig, axs = plt.subplots(1, 2, figsize=(174 * MM, 92 * MM), sharey=True, gridspec_kw={"wspace": 0.12})
for ax, D, ttl, lab in [(axs[0], U, "Univariable MR (IVW)", "a"), (axs[1], MV, "Multivariable MR (triglycerides, BMI and WHR adjusted for BMI)", "b")]:
    for j, o in enumerate(ORD):
        for k, (e, el, c) in enumerate(EX):
            r = D[(D.endpoint == o) & (D.exposure == e)]
            if r.empty: continue
            r = r.iloc[0]; y = -j - (k - 1) * 0.24; sig = r.p < 0.05
            ax.hlines(y, r.lo, r.hi, color=c, lw=0.7); ax.plot(r.OR, y, "o", ms=2.8, mfc=c if sig else "white", mec=c, mew=0.7, label=el if j == 0 else None)
    ax.axvline(1, color="#7F7F7F", ls="--", lw=0.4); ax.axhline(-6.5, color="#BDBDBD", lw=0.4); ax.set_xscale("log")
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.5, 0.7, 1, 1.5, 2, 3, 5])); ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.set_xlabel("OR per SD of genetically predicted exposure (95% CI)")
    ax.set_title(ttl, loc="left", fontsize=6.4, pad=4); ax.text(-0.05 if lab == "b" else -0.36, 1.045, lab, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="right" if lab == "b" else "left")
axs[0].set_yticks([-j for j in range(len(ORD))]); axs[0].set_yticklabels([LBL[o] for o in ORD]); axs[0].set_ylim(-len(ORD) + 0.4, 0.6)
from matplotlib.lines import Line2D
fig.legend([Line2D([], [], marker="o", color=c, mfc=c, ms=3, lw=0.7) for _, _, c in EX], [el for _, el, _ in EX], frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(0.55, -0.04))
fig.text(0.99, -0.03, "Filled: $P$<0.05; open: $P$≥0.05", ha="right", fontsize=5.6)
save(fig, "Fig6"); print("Fig6 saved")
