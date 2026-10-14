"""u10：主图 1–6 与补充图（Arial，PDF + 600 dpi TIFF + 预览 PNG）。"""
import pandas as pd, numpy as np, io, os, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import FancyBboxPatch
from PIL import Image
RES, FIG = "../results", "../figures"; os.makedirs(FIG, exist_ok=True)
MM = 1 / 25.4
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
                     "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "lines.linewidth": 0.9})
RED, BLUE, GREY, LG = "#B2182B", "#2166AC", "#7F7F7F", "#D9D9D9"
LAB = {"TyG": "TyG", "TyG_fs": "TyG (fasting-std.)", "TyG_BMI": "TyG-BMI", "TyG_WC": "TyG-WC", "TyG_WHtR": "TyG-WHtR", "TyG_ABSI": "TyG-ABSI", "TG_HDL": "TG/HDL",
       "AIP": "AIP", "METS_IR": "METS-IR", "neg_eGDR": "eGDR", "neg_eGDR_noHTN": "eGDR without HTN", "CHG": "CHG", "CTI": "CTI", "c_lnTG": "TG", "c_lnGlu": "Glucose",
       "c_negHDL": "HDL-C (inv.)", "c_BMI": "BMI", "c_WC": "Waist", "c_HbA1c": "HbA1c", "ref_ApoB": "ApoB", "HOMA_rin": "HOMA-IR"}
COMPOSITE = ["TyG", "TyG_BMI", "TyG_WC", "TyG_WHtR", "TyG_ABSI", "TG_HDL", "AIP", "METS_IR", "neg_eGDR", "neg_eGDR_noHTN", "CHG", "CTI"]
COMPS = ["c_lnTG", "c_lnGlu", "c_negHDL", "c_BMI", "c_WC", "c_HbA1c", "ref_ApoB"]
TGT = {"fo_E11": "Type 2 diabetes", "fo_K76": "Fatty liver (K76)", "fo_M10": "Gout", "fo_K80": "Gallstones", "fo_K85": "Acute pancreatitis", "fo_G47": "Sleep disorders",
       "fo_I10": "Hypertension", "fo_I21": "Myocardial infarction", "fo_I25": "Chronic IHD", "fo_I70": "Atherosclerosis", "fo_I73": "PVD", "fo_I63": "Ischaemic stroke",
       "fo_I50": "Heart failure", "fo_I48": "Atrial fibrillation", "fo_I35": "Aortic valve", "fo_N18": "CKD", "fo_I26": "Pulm. embolism", "fo_I80": "DVT",
       "fo_F03": "Dementia", "ca_C18": "Colon cancer", "ca_C22": "Liver cancer", "ca_C25": "Pancreatic cancer", "ca_C64": "Kidney cancer", "cvdeath": "CV death", "death": "All-cause death"}
def save(fig, name):
    fig.savefig(f"{FIG}/{name}.pdf", bbox_inches="tight"); buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=600, bbox_inches="tight", facecolor="white"); buf.seek(0)
    im = Image.open(buf).convert("RGB"); im.save(f"{FIG}/{name}.tiff", format="TIFF", compression="tiff_lzw", dpi=(600, 600)); im.resize((im.size[0] // 3, im.size[1] // 3)).save(f"{FIG}/{name}_preview.png"); plt.close(fig)
def plab(ax, t, x=-0.12, y=1.04): ax.text(x, y, t, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
C = pd.read_csv(f"{RES}/u04_calibrated.csv"); S = pd.read_csv(f"{RES}/u04_index_summary.csv")
B = pd.read_csv("../cache/u_base.csv", usecols=[f"z_{i}" for i in COMPOSITE + COMPS])
# ================= 补充图 S7：指数间相关矩阵（原图 1a）；研究设计另存为补充图 S1 =================
H = pd.read_csv(f"{RES}/v01_nhanes_homa.csv")
fig, ax = plt.subplots(figsize=(110 * MM, 92 * MM)); order = COMPOSITE + ["c_lnTG", "c_lnGlu", "c_WC", "c_HbA1c"]
Rm = B[[f"z_{i}" for i in order]].corr().values; im = ax.imshow(Rm, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(order))); ax.set_xticklabels([LAB[i] for i in order], rotation=90); ax.set_yticks(range(len(order))); ax.set_yticklabels([LAB[i] for i in order]); ax.spines[:].set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.ax.tick_params(labelsize=5); ax.set_title("Pearson correlation between indices (UK Biobank, n = 374,814)", loc="left", fontsize=6.5)
save(fig, "SuppFig7_correlation")
fig, ax = plt.subplots(figsize=(90 * MM, 95 * MM)); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
def box(y, h, t, fc="#F2F2F2"):
    ax.add_patch(FancyBboxPatch((3, y), 94, h, boxstyle="round,pad=0.8,rounding_size=2", fc=fc, ec="#595959", lw=0.6)); ax.text(50, y + h / 2, t, ha="center", va="center", fontsize=5.6, linespacing=1.3)
FL = pd.read_csv(f"{RES}/u14_flow.csv").set_index("step").n
box(82, 15, f"UK Biobank participants, n = {FL.iloc[0]:,}\nExcluded: diabetes at baseline (n = {FL.iloc[1]:,});\nmissing index components or covariates (n = {FL.iloc[2]:,})\n"
    f"Analytic cohort without diabetes, n = {FL.iloc[3]:,}; follow-up to 31 Dec 2022", "#DCE9F5")
box(62, 16, "12 insulin-resistance surrogate indices\n(TyG family, TG/HDL, AIP, METS-IR, eGDR, CHG, CTI)\n+ 6 components, ApoB, fasting-standardised TyG")
box(40, 16, "540 incident diseases (ICD-10 A–N + cancer registry),\nall-cause and CV death;\n60 pre-specified negative-control outcomes")
box(3, 31, "Empirical calibration (specificity)\nComponent decomposition\nIncremental C-index; dose–response; subgroups\nMRI organ fat; TyG–ApoB discordance\nNHANES 2005–2018: HOMA-IR anchoring and mortality", "#E5F1E0")
for y1, y2 in [(82, 79.5), (62, 57.5), (40, 35.5)]: ax.annotate("", xy=(50, y2), xytext=(50, y1), arrowprops=dict(arrowstyle="-|>", color="#595959", lw=0.6, mutation_scale=6))
save(fig, "SuppFig1_design")
# ================= 图 1：全表型（3 个指数） =================
SHORT = {"fo_E11": "Type 2 diabetes", "fo_E14": "Unspecified diabetes", "fo_E10": "Insulin-treated diabetes", "fo_E13": "Other diabetes", "fo_H36": "Diabetic retinopathy",
         "fo_G63": "Secondary polyneuropathy", "fo_E66": "Obesity", "fo_E78": "Dyslipidaemia", "fo_M10": "Gout", "fo_K76": "Fatty liver (K76)", "fo_K80": "Gallstones",
         "fo_K85": "Acute pancreatitis", "fo_G47": "Sleep disorders", "fo_I10": "Hypertension", "fo_I11": "Hypertensive heart dis.", "fo_I12": "Hypertensive renal dis.",
         "fo_I20": "Angina", "fo_I21": "Myocardial infarction", "fo_I25": "Chronic IHD", "fo_N18": "CKD", "fo_L97": "Leg ulcer", "fo_K42": "Umbilical hernia",
         "fo_I50": "Heart failure", "fo_K74": "Liver fibrosis/cirrhosis", "fo_K41": "Femoral hernia", "fo_M81": "Osteoporosis", "fo_I89": "Lymphoedema",
         "fo_E29": "Testicular dysfunction", "fo_N03": "Chronic nephritic syndrome", "fo_N17": "Acute kidney injury", "cvdeath": "CV death", "ca_C54": "Uterine cancer"}
def short(o, t): return SHORT.get(o, (str(t)[:1].upper() + str(t)[1:24]).rstrip(" ,"))
def place_labels(ax, pts, fs=5.0):
    """标注带：所有标签排在数据上方同一行；用保序回归（PAVA）在不重叠的前提下让每个标签尽量靠近其圆点，
    再通过两两交换标签顺序消除引线交叉。引线从标签底部连到圆点边缘，不穿过文字。"""
    fig = ax.figure; fig.canvas.draw(); rr = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    y_top = max(y for _, y, _ in pts); y_lab = y_top + 3.2; x_lo, x_hi = ax.get_xlim(); pad = 5
    def width(s):
        t = ax.text(0, 0, s, fontsize=fs); bb = t.get_window_extent(rr); t.remove()
        (x0, _), (x1, _) = inv.transform([(0, 0), (bb.width, 0)]); return x1 - x0
    W = {s: width(s) for _, _, s in pts}
    def layout(P):
        w = [W[s] for _, _, s in P]; off = [0.0]
        for i in range(1, len(P)): off.append(off[-1] + (w[i - 1] + w[i]) / 2 + pad)
        blocks = []
        for v in [x - o for (x, _, _), o in zip(P, off)]:          # 保序回归：q 非降，最小化 Σ(q − t)²
            blocks.append([v, 1])
            while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
                a = blocks.pop(); blocks[-1][0] += a[0]; blocks[-1][1] += a[1]
        cx = [qi + o for qi, o in zip(sum([[b[0] / b[1]] * b[1] for b in blocks], []), off)]
        sh = max(0, (x_lo + w[0] / 2 + 2) - cx[0]) - max(0, cx[-1] - (x_hi - w[-1] / 2 - 2))
        return [c + sh for c in cx]
    def crossings(P, cx):
        def inter(a, b, c, d):
            o = lambda p, q, r: (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
            return o(a, b, c) * o(a, b, d) < 0 and o(c, d, a) * o(c, d, b) < 0
        seg = [((c, y_lab), (x, y)) for (x, y, _), c in zip(P, cx)]
        return sum(inter(*seg[i], *seg[j]) for i in range(len(seg)) for j in range(i + 1, len(seg)))
    P = sorted(pts); cx = layout(P); best = crossings(P, cx); improved = True
    while improved and best > 0:                                   # 交换顺序，直到没有交叉或无法再减少
        improved = False
        for i in range(len(P) - 1):
            for j in range(i + 1, len(P)):                         # 任意两两交换
                Q = P.copy(); Q[i], Q[j] = Q[j], Q[i]; cq = layout(Q); k = crossings(Q, cq)
                if k < best: P, cx, best, improved = Q, cq, k, True
    if best > 0:                                                   # 仍有交叉：在“每个标签最多偏离原位 3 个位置”的排列中搜索
        base = sorted(pts); n = len(base); cand = []
        def rec(prefix, used):
            if len(prefix) == n: cand.append(prefix); return
            i = len(prefix)
            for j in range(max(0, i - 3), min(n, i + 4)):
                if j not in used: rec(prefix + [j], used | {j})
        rec([], frozenset())
        for perm in cand:
            Q = [base[j] for j in perm]; cq = layout(Q); k = crossings(Q, cq)
            if k < best: P, cx, best = Q, cq, k
            if best == 0: break
    for (x, y, s), c in zip(P, cx):
        ax.annotate(s, xy=(x, y), xytext=(c, y_lab), textcoords="data", fontsize=fs, ha="center", va="bottom", zorder=4,
                    arrowprops=dict(arrowstyle="-", color="#404040", lw=0.5, shrinkA=0.8, shrinkB=1.9))
    ax.set_ylim(top=y_lab + 2.4)
    ax.set_yticks([v for v in ax.get_yticks() if v <= y_top + 0.5 and v >= ax.get_ylim()[0]])   # 刻度只覆盖数据区
    return best
chap_order = [c for c in sorted(C.chapter.dropna().unique()) if c != "Death"] + ["Death"]
cmap = plt.get_cmap("tab20"); col = {c: cmap(i % 20) for i, c in enumerate(chap_order)}
fig, axs = plt.subplots(3, 1, figsize=(174 * MM, 150 * MM), sharex=True)
base = C[(C.model == "M1") & (C["index"] == "TyG")].copy(); base["ck"] = base.chapter.map({c: i for i, c in enumerate(chap_order)})
outs = base.sort_values(["ck", "outcome"]).outcome.tolist(); pos = {o: i for i, o in enumerate(outs)}
for ax, ix, lab in zip(axs, ["TyG", "TyG_WC", "neg_eGDR"], "abc"):
    d = C[(C.model == "M1") & (C["index"] == ix)].copy(); d["x"] = d.outcome.map(pos)
    sig = d.fdr_cal < 0.05
    ax.scatter(d.x[~sig], d.z_cal[~sig], s=3, c=[col[c] for c in d.chapter[~sig]], alpha=0.5, lw=0)
    ax.scatter(d.x[sig], d.z_cal[sig], s=9, c=[col[c] for c in d.chapter[sig]], edgecolor="black", lw=0.3)
    ax.axhline(0, color=GREY, lw=0.4); ax.set_ylim(min(-6, d.z_cal.min() - 1), d.z_cal.max() + 2)
    top = d[sig].nlargest(9, "z_cal"); nx = place_labels(ax, [(r.x, r.z_cal, short(r.outcome, r.title)) for _, r in top.iterrows()]); print(ix, "leader-line crossings:", nx)
    ax.set_ylabel("Calibrated z")
    ax.set_title(f"{LAB[ix]}: {int(sig.sum())} of {len(d)} outcomes significant after calibration ({int((d.fdr_raw < 0.05).sum())} before calibration)", loc="left", fontsize=6.5)
    plab(ax, lab, x=-0.07)
cc = base.assign(x=lambda d: d.outcome.map(pos)).groupby("chapter").x.mean().reindex(chap_order)
axs[-1].set_xticks(cc.values); axs[-1].set_xticklabels(cc.index, fontsize=6); axs[-1].set_xlabel("ICD-10 chapter (C and D include cancer registry; Death = all-cause and CV death)")
fig.tight_layout(); save(fig, "Fig1")
# ================= 图 2：阴性对照揭示偏倚（a 校准前后疾病数；b 系统误差背景；c 复现常规分析） =================
s1 = S[S.model == "M1"].set_index("index").loc[COMPOSITE + COMPS]; s1 = s1.sort_values("sig_calibrated_FDR")
SNc = pd.read_csv(f"{RES}/u04b_strict_nc.csv"); SNc = SNc[SNc.model == "M1"].set_index("index").reindex(s1.index)
fig = plt.figure(figsize=(174 * MM, 145 * MM)); gs = gridspec.GridSpec(2, 2, height_ratios=[1.15, 1], width_ratios=[1.3, 1], hspace=0.72, wspace=0.55)
y = np.arange(len(s1))
ax = fig.add_subplot(gs[0, 0]); ax.hlines(y, s1.sig_calibrated_FDR, s1.sig_raw_FDR, color=LG, lw=1.6)
ax.plot(s1.sig_raw_FDR, y, "o", color=GREY, ms=3, label="Before calibration")
ax.plot(s1.sig_calibrated_FDR, y, "o", color=RED, ms=3, label="Calibrated (60 negative controls)")
ax.plot(SNc.sig_cal, y, "o", mfc="none", mec=RED, ms=3.4, mew=0.6, label="Calibrated (40 strict negative controls)")
ax.set_yticks(y); ax.set_yticklabels([LAB[i] for i in s1.index]); ax.set_xlabel("Diseases with FDR<0.05 (of 482 non-control outcomes)")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.2), ncol=2, fontsize=5.5, columnspacing=1.2, handletextpad=0.3); plab(ax, "a", x=-0.3)
ax = fig.add_subplot(gs[0, 1])
ax.errorbar(s1.NC_mean_HR, y, xerr=[s1.NC_mean_HR - np.exp(np.log(s1.NC_mean_HR) - 1.96 * s1.NC_tau), np.exp(np.log(s1.NC_mean_HR) + 1.96 * s1.NC_tau) - s1.NC_mean_HR],
            fmt="s", color="black", ms=2.5, elinewidth=0.6, capsize=0); ax.axvline(1, color=GREY, ls="--", lw=0.5)
ax.set_yticks(y); ax.set_yticklabels([]); ax.set_xlabel("Negative-control background HR per SD\n(mean and 95% range of systematic error)"); plab(ax, "b", x=-0.08)
CV = pd.read_csv(f"{RES}/u15_conventional.csv")
GO = [("HF", "Heart failure"), ("Stroke", "Ischaemic stroke"), ("Death", "All-cause death"), ("AF", "Atrial fibrillation"),
      ("Acute tonsillitis (NC)", "Acute tonsillitis"), ("Dental caries (NC)", "Dental caries"), ("MI", "Myocardial infarction")]
MD = [("A", "#7F7F7F", "Demographics and lifestyle"), ("B", "#E08214", "+ blood pressure, LDL-C, drugs"), ("C", RED, "+ BMI")]
ax = fig.add_subplot(gs[1, :]); ax.axvspan(3.5, 5.5, color="#F2F2F2", lw=0, zorder=0); ax.axvspan(5.5, 6.5, color="#EAF1F8", lw=0, zorder=0)
for j, (o, _) in enumerate(GO):
    for k, (m, c, lb) in enumerate(MD):
        r = CV[(CV.outcome == o) & (CV.model == m)].iloc[0]; xj = j + (k - 1) * 0.22; filled = r.p < 0.05
        ax.vlines(xj, r.lo, r.hi, color=c, lw=0.8); ax.plot(xj, r.HR_Q4, "o", ms=3.4, mfc=c if filled else "white", mec=c, mew=0.8, label=lb if j == 0 else None)
ax.axhline(1, color=GREY, ls="--", lw=0.5); ax.set_yscale("log"); ax.set_xticks(range(len(GO))); ax.set_xticklabels([g[1] for g in GO], rotation=25, ha="right")
ax.set_yticks([0.7, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 1.75, 2]); ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_ylim(0.68, 2.25); ax.set_xlim(-0.6, len(GO) - 0.4); ax.set_ylabel("TyG Q4 vs Q1: HR (95% CI)")
ax.text(4.5, 2.05, "Negative controls", ha="center", fontsize=5.8, color="#595959"); ax.text(6, 2.05, "Positive control", ha="center", fontsize=5.8, color="#2166AC")
ax.legend(frameon=False, loc="upper left", fontsize=5.6)
ax.set_title("Conventional quartile analysis: filled symbols $P$<0.05", loc="left", fontsize=6); plab(ax, "c", x=-0.07)
save(fig, "Fig2")
# ================= 图 3：TyG 的特异谱系（a 关键结局校准 HR 热图；b TyG 调整 BMI 前后） =================
T = C[(C.model == "M1") & C.outcome.isin(TGT)].copy(); rows = COMPOSITE + COMPS
W = T.pivot(index="outcome", columns="index", values="HR_cal").reindex(index=list(TGT), columns=rows); Q = T.pivot(index="outcome", columns="index", values="fdr_cal").reindex(index=list(TGT), columns=rows)
Qr = T.pivot(index="outcome", columns="index", values="fdr_raw").reindex(index=list(TGT), columns=rows)
fig = plt.figure(figsize=(174 * MM, 185 * MM)); gs = gridspec.GridSpec(2, 1, height_ratios=[1.55, 1], hspace=0.62)
ax = fig.add_subplot(gs[0]); im = ax.imshow(np.log(W.values.T), cmap="RdBu_r", vmin=np.log(0.6), vmax=np.log(1 / 0.6), aspect="auto")
for i in range(W.shape[1]):
    for j in range(W.shape[0]):
        if Q.values[j, i] < 0.05: ax.text(j, i, "●", ha="center", va="center", fontsize=4.5, color="black")
        elif Qr.values[j, i] < 0.05: ax.text(j, i, "○", ha="center", va="center", fontsize=4.5, color="black")
ax.set_xticks(range(len(TGT))); ax.set_xticklabels(list(TGT.values()), rotation=60, ha="right"); ax.set_yticks(range(len(rows))); ax.set_yticklabels([LAB[r] for r in rows]); ax.axhline(len(COMPOSITE) - 0.5, color="black", lw=0.6)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, ticks=np.log([0.6, 0.8, 1, 1.25, 1.67])); cb.ax.set_yticklabels(["0.6", "0.8", "1.0", "1.25", "1.67"]); cb.set_label("Calibrated HR per SD", fontsize=6)
ax.set_title("● calibrated FDR<0.05   ○ significant before calibration only", loc="left", fontsize=6); plab(ax, "a", x=-0.12, y=1.03)
KO = ["fo_E11", "fo_K76", "fo_M10", "fo_K80", "fo_K85", "fo_I21", "fo_I25", "fo_N18", "fo_I63", "fo_I50", "fo_I48", "cvdeath", "death"]
U = pd.read_csv(f"{RES}/u03_phewas.csv"); U = U[(U["index"] == "TyG") & U.outcome.isin(KO)].set_index(["outcome", "model"])
Cc = C[(C["index"] == "TyG") & C.outcome.isin(KO)].set_index(["outcome", "model"])
ax = fig.add_subplot(gs[1]); xx = np.arange(len(KO))
for k, (mdl, c, lb) in enumerate([("M1", GREY, "Model 1 (demographics, lifestyle)"), ("M2", RED, "Model 2 (+ BMI)")]):
    for j, o in enumerate(KO):
        r = U.loc[(o, mdl)]; hr, lo, hi = np.exp(r.logHR), np.exp(r.logHR - 1.96 * r.se), np.exp(r.logHR + 1.96 * r.se)
        xj = j + (-0.15 if k == 0 else 0.15); filled = Cc.loc[(o, mdl), "fdr_cal"] < 0.05
        ax.vlines(xj, lo, hi, color=c, lw=0.8); ax.plot(xj, hr, "o", ms=3.4, mfc=c if filled else "white", mec=c, mew=0.8, label=lb if j == 0 else None)
ax.axhline(1, color=GREY, ls="--", lw=0.5); ax.set_yscale("log"); ax.set_xticks(xx); ax.set_xticklabels([TGT[o] for o in KO], rotation=35, ha="right")
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_yticks([0.8, 0.9, 1, 1.25, 1.5, 2])
ax.set_ylabel("TyG: HR per SD (95% CI)"); ax.legend(frameon=False, loc="upper right", fontsize=5.6)
ax.set_title("Filled symbols: significant after negative-control calibration (FDR<0.05); open symbols: not specific", loc="left", fontsize=6); plab(ax, "b", x=-0.07)
save(fig, "Fig3")
# ================= 图 4：成分与 MRI =================
Dc = pd.read_csv(f"{RES}/u05_decompose.csv"); Dc = Dc[Dc.outcome.isin(TGT)]; Dc["z"] = Dc.logHR / Dc.se
cm = ["lnTG", "lnGlu", "negHDL", "WC", "HbA1c", "ApoB"]; cml = ["TG", "Glucose", "HDL-C (inv.)", "Waist", "HbA1c", "ApoB"]
Wz = Dc.pivot(index="outcome", columns="component", values="logHR").reindex(index=list(TGT), columns=cm)
Mr = pd.read_csv(f"{RES}/u08_mri.csv"); tr = ["Liver PDFF", "Pancreas PDFF", "VAT", "ASAT", "Muscle fat infiltration", "LV mass index", "LVEF"]
Wm = Mr.pivot(index="index", columns="trait", values="beta").reindex(index=COMPOSITE + ["c_lnTG", "c_lnGlu", "c_WC", "c_HbA1c", "ref_ApoB"], columns=tr)
fig, axs = plt.subplots(1, 2, figsize=(174 * MM, 92 * MM), gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.55})
ax = axs[0]; im = ax.imshow(np.exp(Wz.values), cmap="RdBu_r", vmin=0.7, vmax=1 / 0.7, aspect="auto")
for i in range(Wz.shape[0]):
    for j in range(Wz.shape[1]): v = np.exp(Wz.values[i, j]); ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.4, color="white" if (v > 1.36 or v < 0.78) else "black")
ax.set_yticks(range(len(TGT))); ax.set_yticklabels(list(TGT.values())); ax.set_xticks(range(len(cm))); ax.set_xticklabels(cml, rotation=40, ha="right"); plab(ax, "a", x=-0.45)
ax.set_title("Mutually adjusted HR per SD of each component", loc="left", fontsize=6.3)
ax = axs[1]; im2 = ax.imshow(Wm.values, cmap="Reds", vmin=0, vmax=0.8, aspect="auto")
for i in range(Wm.shape[0]):
    for j in range(Wm.shape[1]): ax.text(j, i, f"{Wm.values[i, j]:.2f}", ha="center", va="center", fontsize=4.2, color="white" if Wm.values[i, j] > 0.5 else "black")
ax.set_yticks(range(Wm.shape[0])); ax.set_yticklabels([LAB[i] for i in Wm.index]); ax.set_xticks(range(len(tr))); ax.set_xticklabels(tr, rotation=40, ha="right"); plab(ax, "b", x=-0.6)
ax.set_title("MRI 9–11 y later: β (SD) per SD of index", loc="left", fontsize=6.3); save(fig, "Fig4")
# ================= 图 5：构念效度与预测（a 与 HOMA-IR 的相关；b 临床增量 C；c NHANES 心血管死亡） =================
P6 = pd.read_csv(f"{RES}/u06_prediction.csv"); P6 = P6[P6.base == "B"]
oc = ["T2D", "MASLD", "Gout", "Gallstones", "MI", "MACE", "HF", "AF", "Stroke", "CKD", "Death"]; ix6 = ["TyG", "TyG_WC", "METS_IR", "TG_HDL", "CHG", "CTI", "neg_eGDR", "neg_eGDR_noHTN", "c_HbA1c", "ref_ApoB"]
Wd = P6.pivot(index="index", columns="outcome", values="dC").reindex(index=ix6, columns=oc) * 1000
Mo = pd.read_csv(f"{RES}/v01_nhanes_mortality.csv")
fig = plt.figure(figsize=(174 * MM, 150 * MM)); gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1.55], height_ratios=[1.3, 1], hspace=0.62, wspace=0.62)
ax = fig.add_subplot(gs[0, 0]); plab(ax, "a", x=-0.5)
h = H[H.group == "All"].set_index("index").loc[[i for i in COMPOSITE + ["c_lnTG", "c_lnGlu", "c_WC", "c_HbA1c"] if i in set(H["index"])]].sort_values("spearman_r_HOMA")
y = np.arange(len(h)); ax.barh(y, h.spearman_r_HOMA, color=LG, edgecolor=GREY, lw=0.4)
for grp, c, mk in [("NH White", BLUE, "o"), ("NH Black", RED, "s")]:
    g = H[H.group == grp].set_index("index").reindex(h.index); ax.plot(g.spearman_r_HOMA, y, mk, color=c, ms=2.2, label=grp)
ax.set_yticks(y); ax.set_yticklabels([LAB[i] for i in h.index]); ax.set_xlabel("Spearman r with HOMA-IR (NHANES)"); ax.legend(frameon=False, loc="lower right", fontsize=5.4); ax.set_xlim(0, 0.75)
ax.set_title("Construct validity against HOMA-IR", loc="left", fontsize=6.3)
ax = fig.add_subplot(gs[0, 1]); im = ax.imshow(Wd.values, cmap="Purples", vmin=0, vmax=25, aspect="auto")
for i in range(Wd.shape[0]):
    for j in range(Wd.shape[1]): ax.text(j, i, f"{Wd.values[i, j]:.1f}", ha="center", va="center", fontsize=4.6, color="white" if Wd.values[i, j] > 15 else "black")
ax.set_xticks(range(len(oc))); ax.set_xticklabels([{"MASLD": "Fatty liver"}.get(o, o) for o in oc], rotation=40, ha="right"); ax.set_yticks(range(len(ix6))); ax.set_yticklabels([LAB[i] for i in ix6]); plab(ax, "b", x=-0.3)
ax.set_title("ΔC-index ×1000 over clinical model (UK Biobank)", loc="left", fontsize=6.3)
ax = fig.add_subplot(gs[1, :]); mi = ["TyG", "TyG_WC", "METS_IR", "neg_eGDR", "neg_eGDR_noHTN", "CHG", "c_HbA1c", "HOMA_rin"]
m = Mo[(Mo.outcome == "cvd")].copy(); xx = np.arange(len(mi))
for k, (md, c) in enumerate([("base", GREY), ("+HOMA", RED)]):
    g = m[m.model == md].set_index("index").reindex(mi); xk = xx + (-0.13 if k == 0 else 0.13)
    ax.vlines(xk, g.lo, g.hi, color=c, lw=0.8); ax.plot(xk, g.HR, "s", color=c, ms=2.8, label={"base": "Adjusted for age, sex, race and ethnicity", "+HOMA": "+ HOMA-IR"}[md])
ax.axhline(1, color=GREY, ls="--", lw=0.5); ax.set_xticks(xx); ax.set_xticklabels([LAB[i] for i in mi]); ax.set_yscale("log")
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_yticks([0.8, 0.9, 1, 1.25, 1.5])
ax.set_ylim(0.8, 1.85); ax.set_ylabel("CV death: HR per SD (95% CI)"); ax.legend(frameon=False, loc="upper right", fontsize=5.6, ncol=2)
ax.set_title("NHANES 2005–2018: 10,714 participants, 276 cardiovascular deaths", loc="left", fontsize=6.3); plab(ax, "c", x=-0.07)
save(fig, "Fig5")
print("figures saved")
