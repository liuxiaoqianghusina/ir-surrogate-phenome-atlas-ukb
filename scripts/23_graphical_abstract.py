"""ga01：图形摘要（Cardiovascular Diabetology：约 920×300 像素，图注 ≤30 词，须在图片内、位于图下方）。
四栏：研究设计 → 校准前后显著疾病数 → TyG 特异与不特异的结局 → 孟德尔随机化。输出 920×300 和 3 倍分辨率 PNG。"""
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
RES, FIG = "../results", "../figures"
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})
RED, BLUE, GREY, DARK, GREEN = "#B2182B", "#2166AC", "#9E9E9E", "#262626", "#1B7837"
S = pd.read_csv(f"{RES}/u04_index_summary.csv"); S = S[S.model == "M1"].set_index("index")
U = pd.read_csv(f"{RES}/mr03_univariable.csv"); U = U[U.method == "IVW"].set_index(["exposure", "endpoint"])
fig = plt.figure(figsize=(9.2001, 3.0001), dpi=100); fig.patch.set_facecolor("white")
W = fig.add_axes([0, 0, 1, 1]); W.set_xlim(0, 920); W.set_ylim(0, 300); W.axis("off")
CAP = ("Most disease associations of insulin resistance surrogates reflect systematic error and adiposity; TyG is specific for metabolic, "
       "hepatobiliary and coronary disease, not heart failure, atrial fibrillation, stroke or death.")
assert len(CAP.split()) <= 30
W.add_patch(FancyBboxPatch((0, 0), 920, 46, boxstyle="square,pad=0", fc="#F2F2F2", ec="none"))
W.text(460, 23, CAP, ha="center", va="center", fontsize=8.0, color=DARK, wrap=True)
col = [(8, 205), (225, 230), (467, 225), (704, 208)]   # 各栏 x 起点与宽度
TOP, BOT = 292, 56
for k, (x, w) in enumerate(col):
    W.add_patch(FancyBboxPatch((x, BOT), w, TOP - BOT, boxstyle="round,pad=0,rounding_size=8", fc="white", ec="#BDBDBD", lw=0.8))
    if k < 3:
        W.add_patch(FancyArrowPatch((x + w + 1, (TOP + BOT) / 2), (col[k + 1][0] - 1, (TOP + BOT) / 2), arrowstyle="-|>", mutation_scale=9, color="#757575", lw=1.0))
heads = ["Design", "Broad but not specific", "Specific profile of TyG", "Genetic evidence (MR)"]
for (x, w), h in zip(col, heads): W.text(x + w / 2, TOP - 15, h, ha="center", va="center", fontsize=9.5, fontweight="bold", color=DARK)
# ---- 1 研究设计
x0 = col[0][0] + 12; y = TOP - 42
for big, small in [("374,814", "UK Biobank adults\nwithout diabetes"), ("12", "insulin resistance\nsurrogate indices"), ("540", "incident diseases;\n60 negative controls")]:
    W.text(x0, y, big, fontsize=12, fontweight="bold", color=BLUE, va="center")
    W.text(x0 + 66, y, small, fontsize=7.4, color=DARK, va="center", linespacing=1.15); y -= 50
W.text(x0, BOT + 14, "+ NHANES (HOMA-IR)\n+ FinnGen (Mendelian randomisation)", fontsize=6.8, color="#616161", va="center", linespacing=1.2)
# ---- 2 校准前后
ax = fig.add_axes([0.315, 0.27, 0.165, 0.54]); rows = [("TyG", "TyG"), ("TyG_WC", "TyG-WC"), ("METS_IR", "METS-IR"), ("TyG_BMI", "TyG-BMI")]
yy = np.arange(len(rows))[::-1]
for yv, (ix, lb) in zip(yy, rows):
    raw, cal = S.loc[ix, "sig_raw_FDR"], S.loc[ix, "sig_calibrated_FDR"]
    ax.barh(yv + 0.18, raw, height=0.34, color="#D9D9D9"); ax.barh(yv - 0.18, cal, height=0.34, color=RED)
    ax.text(raw + 6, yv + 0.18, f"{int(raw)}", va="center", fontsize=6.6, color="#616161"); ax.text(cal + 6, yv - 0.18, f"{int(cal)}", va="center", fontsize=6.6, color=RED, fontweight="bold")
ax.set_yticks(yy); ax.set_yticklabels([r[1] for r in rows], fontsize=7.6); ax.set_xlim(0, 450); ax.set_xticks([])
for sp in ["top", "right", "bottom"]: ax.spines[sp].set_visible(False)
ax.spines["left"].set_linewidth(0.6); ax.tick_params(length=0); ax.patch.set_alpha(0)
ax.set_title("Diseases with FDR <0.05", fontsize=7.2, color=DARK, pad=3)
W.add_patch(FancyBboxPatch((col[1][0] + 28, BOT + 9), 10, 8, boxstyle="square,pad=0", fc="#D9D9D9", ec="none")); W.text(col[1][0] + 42, BOT + 13, "before", fontsize=7.0, va="center")
W.add_patch(FancyBboxPatch((col[1][0] + 100, BOT + 9), 10, 8, boxstyle="square,pad=0", fc=RED, ec="none")); W.text(col[1][0] + 114, BOT + 13, "after calibration", fontsize=7.0, va="center")
# ---- 3 TyG 特异谱系
x0, w = col[2]
W.add_patch(FancyBboxPatch((x0 + 10, 148), w - 20, 112, boxstyle="round,pad=0,rounding_size=5", fc="#EAF4EC", ec="none"))
W.text(x0 + 20, 248, "Specific", fontsize=8, fontweight="bold", color=GREEN, va="center")
for i, t in enumerate(["Type 2 diabetes", "Fatty liver, gallstones, pancreatitis", "Gout", "Myocardial infarction, chronic IHD", "Chronic kidney disease"]):
    W.text(x0 + 24, 230 - i * 16, "•  " + t, fontsize=7.4, color=DARK, va="center")
W.add_patch(FancyBboxPatch((x0 + 10, 64), w - 20, 76, boxstyle="round,pad=0,rounding_size=5", fc="#F5F5F5", ec="none"))
W.text(x0 + 20, 128, "Not specific", fontsize=8, fontweight="bold", color="#616161", va="center")
W.text(x0 + 24, 111, "•  Heart failure, atrial fibrillation", fontsize=7.4, color=DARK, va="center")
W.text(x0 + 24, 95, "•  Stroke, all-cause death", fontsize=7.4, color=DARK, va="center")
W.text(x0 + 24, 76, "HF: HR 1.05 → 0.94 after BMI adjustment", fontsize=6.9, color=RED, va="center")
# ---- 4 孟德尔随机化
ax = fig.add_axes([0.80, 0.30, 0.17, 0.48])
MR = [("TG", "I9_CHD", "Triglycerides → CHD", RED), ("TG", "I9_HEARTFAIL", "Triglycerides → HF", RED), ("BMI", "I9_HEARTFAIL", "BMI → HF", BLUE)]
for j, (e, o, lb, c) in enumerate(MR):
    r = U.loc[(e, o)]; yv = len(MR) - 1 - j
    ax.hlines(yv, r.lo, r.hi, color=c, lw=1.4); ax.plot(r.OR, yv, "o", ms=5, color=c, mec="white", mew=0.6)
    ax.text(0.62, yv + 0.33, lb, fontsize=7.2, color=DARK, va="bottom"); ax.text(2.05, yv, f"{r.OR:.2f}", fontsize=7, color=c, va="center", ha="left", fontweight="bold")
ax.axvline(1, color=GREY, ls="--", lw=0.7); ax.set_xscale("log"); ax.set_xlim(0.6, 2.4); ax.set_ylim(-0.6, len(MR) - 0.1)
ax.set_xticks([0.7, 1, 1.5, 2]); ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.tick_params(labelsize=6.2, length=2); ax.set_yticks([])
for sp in ["top", "right", "left"]: ax.spines[sp].set_visible(False)
ax.set_xlabel("OR per SD (FinnGen)", fontsize=6.4, labelpad=1)
fig.savefig(f"{FIG}/Graphical_abstract.png", dpi=100, facecolor="white")
fig.savefig(f"{FIG}/Graphical_abstract_hires.png", dpi=300, facecolor="white")
fig.savefig(f"{FIG}/Graphical_abstract.svg", facecolor="white")
from PIL import Image
im = Image.open(f"{FIG}/Graphical_abstract.png"); print("size", im.size, "caption words", len(CAP.split()))
