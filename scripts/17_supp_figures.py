"""u12：补充图 —— S2 剂量反应（TyG、TyG-WC、eGDR、eGDR 去高血压项）；S3 按种族亚组；S4 TyG–ApoB 不一致；S5 M1 vs M2（调整 BMI）特异性对比。"""
import pandas as pd, numpy as np, io, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.transforms
from PIL import Image
RES, FIG = "../results", "../figures"; MM = 1 / 25.4
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "axes.titlesize": 6.8, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 5.8,
                     "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
def save(fig, name):
    fig.savefig(f"{FIG}/{name}.pdf", bbox_inches="tight"); buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=600, bbox_inches="tight", facecolor="white"); buf.seek(0)
    im = Image.open(buf).convert("RGB"); im.save(f"{FIG}/{name}.tiff", format="TIFF", compression="tiff_lzw", dpi=(600, 600)); im.resize((im.size[0] // 3, im.size[1] // 3)).save(f"{FIG}/{name}_preview.png"); plt.close(fig)
LAB = {"c_lnGlu": "Glucose", "c_lnTG": "TG", "c_negHDL": "HDL-C (inv.)", "c_HbA1c": "HbA1c", "ref_ApoB": "ApoB", "AIP": "AIP", "TG_HDL": "TG/HDL", "CTI": "CTI", "TyG": "TyG", "TyG_fs": "TyG (fasting-std.)", "TyG_WC": "TyG-WC", "neg_eGDR": "eGDR", "neg_eGDR_noHTN": "eGDR without HTN", "METS_IR": "METS-IR", "TG_HDL": "TG/HDL", "CHG": "CHG"}
COL = {"TyG": "#B2182B", "TyG_WC": "#E08214", "neg_eGDR": "#2166AC", "neg_eGDR_noHTN": "#542788"}
R = pd.read_csv(f"{RES}/u07_rcs.csv"); OUT = ["T2D", "MASLD", "Gout", "Gallstones", "MI", "HF", "AF", "Stroke", "CKD", "Death"]
fig, axs = plt.subplots(2, 5, figsize=(174 * MM, 80 * MM), sharex=True)
for ax, o in zip(axs.flat, OUT):
    for ix, c in COL.items():
        d = R[(R.outcome == o) & (R["index"] == ix)]; ax.plot(d.x, d.HR, color=c, lw=0.9, label=LAB[ix]); ax.fill_between(d.x, d.lo, d.hi, color=c, alpha=0.12, lw=0)
    ax.axhline(1, color="#7F7F7F", lw=0.4, ls="--"); ax.set_yscale("log"); ax.set_title(o)
    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10])); ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
for ax in axs[-1]: ax.set_xlabel("Index (SD)")
for ax in axs[:, 0]: ax.set_ylabel("HR (ref = mean)")
fig.tight_layout(rect=(0, 0.06, 1, 1)); h_, l_ = axs[0, 0].get_legend_handles_labels(); fig.legend(h_, l_, frameon=False, loc="lower center", ncol=4); save(fig, "SuppFig2_dose_response")
S = pd.read_csv(f"{RES}/u07_subgroups.csv"); e = S[(S.subgroup == "ethgrp") & S["index"].isin(["TyG", "neg_eGDR", "TyG_WC"])]
fig, axs = plt.subplots(1, 3, figsize=(174 * MM, 85 * MM), sharey=True)
for ax, ix in zip(axs, ["TyG", "TyG_WC", "neg_eGDR"]):
    d = e[e["index"] == ix]; y0 = 0; yt, yl = [], []
    for o in OUT:
        dd = d[d.outcome == o]
        for k, (lv, c) in enumerate([("White", "#7F7F7F"), ("South Asian", "#E08214"), ("Black", "#2166AC")]):
            r = dd[dd.level == lv]
            if len(r): r = r.iloc[0]; ax.hlines(-y0 - k * 0.25, r.lo, r.hi, color=c, lw=0.7); ax.plot(r.HR, -y0 - k * 0.25, "s", color=c, ms=2, label=lv if o == OUT[0] else None)
        pint = dd.p_int.iloc[0] if len(dd) else np.nan; yt.append(-y0 - 0.25); yl.append(o)
        ax.text(1.0, -y0 - 0.25, f"{pint:.2g}", transform=matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData), fontsize=5.4, ha="right", va="center"); y0 += 1
    ax.text(1.0, 0.55, "$P$ for interaction", transform=matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData), fontsize=5.4, ha="right", va="center")
    ax.axvline(1, color="#7F7F7F", ls="--", lw=0.4); ax.set_xscale("log"); ax.set_title(LAB[ix], pad=10); ax.set_xlabel("HR per SD"); ax.set_xlim(0.5, 5.2); ax.set_xticks([0.5, 1, 2, 4])
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axs[0].set_yticks(yt); axs[0].set_yticklabels(yl); fig.tight_layout(rect=(0, 0.05, 1, 1))
h_, l_ = axs[0].get_legend_handles_labels(); fig.legend(h_, l_, frameon=False, loc="lower center", ncol=3, fontsize=5.8); save(fig, "SuppFig3_ethnicity")
A = pd.read_csv(f"{RES}/u09_apob_discordance.csv"); A = A[A.term.str.contains("/")]
fig, ax = plt.subplots(figsize=(110 * MM, 85 * MM)); outs = A.outcome.unique(); y = 0; yt, yl = [], []
for o in outs:
    d = A[A.outcome == o]
    for k, (t, c) in enumerate(zip(["TyG high / ApoB low", "TyG low / ApoB high", "TyG high / ApoB high"], ["#B2182B", "#2166AC", "#542788"])):
        r = d[d.term == t].iloc[0]; ax.hlines(-y - k * 0.25, r.lo, r.hi, color=c, lw=0.7); ax.plot(r.HR, -y - k * 0.25, "s", color=c, ms=2, label=t if o == outs[0] else None)
    yt.append(-y - 0.25); yl.append(o); y += 1
ax.axvline(1, color="#7F7F7F", ls="--", lw=0.4); ax.set_xscale("log"); ax.set_yticks(yt); ax.set_yticklabels(yl); ax.set_xlabel("HR vs TyG low / ApoB low (median split)")
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.legend(frameon=False, loc="lower right")
save(fig, "SuppFig4_apob_discordance")
Sm = pd.read_csv(f"{RES}/u04_index_summary.csv"); p = Sm.pivot(index="index", columns="model", values="sig_calibrated_FDR").dropna()
fig, ax = plt.subplots(figsize=(90 * MM, 70 * MM)); y = np.arange(len(p)); p = p.sort_values("M2")
ax.barh(y + 0.2, p.M1, height=0.38, color="#BDBDBD", label="Model 1"); ax.barh(y - 0.2, p.M2, height=0.38, color="#B2182B", label="Model 2 (+ BMI)")
ax.set_yticks(y); ax.set_yticklabels([LAB.get(i, i.replace("c_", "").replace("ref_", "")) for i in p.index]); ax.set_xlabel("Diseases significant after calibration"); ax.legend(frameon=False)
save(fig, "SuppFig5_bmi_adjusted")
# S6：阴性对照结局的估计值（M1）与拟合的系统误差分布
U = pd.read_csv(f"{RES}/u03_phewas.csv"); NCs = pd.read_csv(f"{RES}/u04_negative_controls.csv").outcome.tolist(); meta = pd.read_csv("../cache/phe_meta.csv").set_index("col")
ST = ["H52", "H61", "H00", "H01", "H02", "H04", "H10", "H11", "H16", "H18", "H43", "H65", "H66", "H93", "J01", "J02", "J03", "J30", "J31", "J33", "J34", "J35",
      "K02", "K04", "K08", "K12", "K13", "K35", "L20", "L57", "L60", "L70", "L72", "L80", "L82", "N60", "N43", "N47", "D22", "D23", "G43", "G44"]
fig, axs = plt.subplots(1, 4, figsize=(174 * MM, 150 * MM), sharey=True)
u0 = U[(U.model == "M1") & (U["index"] == "TyG") & U.outcome.isin(NCs)].sort_values("logHR"); order = u0.outcome.tolist(); yy = {o: i for i, o in enumerate(order)}
for ax, ix in zip(axs, ["TyG", "TyG_WC", "neg_eGDR", "c_lnGlu"]):
    d = U[(U.model == "M1") & (U["index"] == ix) & U.outcome.isin(NCs)].copy(); d["y"] = d.outcome.map(yy)
    mu, tau = Sm[(Sm["index"] == ix) & (Sm.model == "M1")][["NC_mean_HR", "NC_tau"]].iloc[0]
    ax.axvspan(np.exp(np.log(mu) - 1.96 * tau), np.exp(np.log(mu) + 1.96 * tau), color="#FDDBC7", lw=0); ax.axvline(mu, color="#B2182B", lw=0.6)
    for _, r in d.iterrows():
        strict = r.outcome[3:] in ST; c = "black" if strict else "#7F7F7F"
        ax.hlines(r.y, np.exp(r.logHR - 1.96 * r.se), np.exp(r.logHR + 1.96 * r.se), color=c, lw=0.5); ax.plot(np.exp(r.logHR), r.y, "o", ms=1.8, mfc=c if strict else "white", mec=c, mew=0.5)
    ax.axvline(1, color="#7F7F7F", ls="--", lw=0.4); ax.set_xscale("log"); ax.set_title(f"{LAB.get(ix, 'Glucose')}\nbackground HR {mu:.3f}", fontsize=6.3)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xticks([0.8, 1, 1.25]); ax.set_xlabel("HR per SD")
axs[0].set_yticks(range(len(order))); axs[0].set_yticklabels([(str(meta.loc[o, "title"])[:1].upper() + str(meta.loc[o, "title"])[1:34]) + f" ({o[3:]})" for o in order], fontsize=4.6)
fig.tight_layout(); save(fig, "SuppFig6_negative_controls"); print("supp figs saved")
