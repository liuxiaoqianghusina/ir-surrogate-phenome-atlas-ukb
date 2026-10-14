"""u04：阴性对照经验校准与特异性指标（全表型）。
系统误差模型：对每个指数×模型，用阴性对照估计 logHR_i ~ N(mu, tau² + se_i²)（极大似然），校准 z = (b − mu)/sqrt(se² + tau²)。
输出：u04_calibrated.csv（逐结局）、u04_index_summary.csv（各指数特异性）、u04_best_index.csv（每个结局最佳指数）、u04_targets.csv（预设目标结局矩阵）。"""
import pandas as pd, numpy as np
from scipy import stats, optimize
R = pd.read_csv("../results/u03_phewas.csv"); meta = pd.read_csv("../cache/phe_meta.csv")
meta["outcome"] = meta.col
TITLE = dict(zip(meta.col, meta.title)); CHAP = dict(zip(meta.col, meta.chapter))
TITLE.update({"death": "All-cause death", "cvdeath": "Cardiovascular death"}); CHAP.update({"death": "Death", "cvdeath": "Death"})
# 预设阴性对照：与胰岛素抵抗无既定生物学联系的眼、耳、口腔、鼻咽、皮肤、肛门、乳腺/外生殖器良性病变及局部机械性疾病（分析前确定）
NC_CAND = ["H52", "H61", "H00", "H01", "H02", "H04", "H10", "H11", "H16", "H18", "H33", "H43", "H53", "H65", "H66", "H81", "H83", "H93",
           "J01", "J02", "J03", "J30", "J31", "J32", "J33", "J34", "J35", "K02", "K04", "K08", "K12", "K13", "K35", "K40", "K60", "K61", "K62", "K64",
           "L20", "L23", "L30", "L57", "L60", "L70", "L72", "L80", "L82", "L90", "M20", "M21", "M22", "M24", "N60", "N62", "N84", "N87", "N43", "N47", "G43", "G44", "D22", "D23"]
NC = [c for c in R.outcome.unique() if c[3:] in NC_CAND]
print("negative controls available:", len(NC))
pd.Series(NC).to_csv("../results/u04_negative_controls.csv", index=False, header=["outcome"])
def fit_null(b, se):
    nll = lambda th: -np.sum(stats.norm.logpdf(b, th[0], np.sqrt(se ** 2 + np.exp(th[1]) ** 2)))
    r = optimize.minimize(nll, [np.mean(b), np.log(max(np.std(b), 1e-3))], method="Nelder-Mead"); return r.x[0], np.exp(r.x[1])
def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = np.empty_like(p); n = len(p); c = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; q[o] = np.minimum(c, 1); return q
out, summ = [], []
for (ix, mdl), g in R.groupby(["index", "model"]):
    nc = g[g.outcome.isin(NC)]; mu, tau = fit_null(nc.logHR.values, nc.se.values)
    t = g[~g.outcome.isin(NC)].copy()
    t["z_raw"] = t.logHR / t.se; t["z_cal"] = (t.logHR - mu) / np.sqrt(t.se ** 2 + tau ** 2)
    t["p_cal"] = 2 * stats.norm.sf(t.z_cal.abs()); t["fdr_raw"] = bh(t.p.values); t["fdr_cal"] = bh(t.p_cal.values)
    t["HR"] = np.exp(t.logHR); t["HR_cal"] = np.exp(t.logHR - mu)
    t["cal_lo"] = np.exp(t.logHR - mu - 1.96 * np.sqrt(t.se ** 2 + tau ** 2)); t["cal_hi"] = np.exp(t.logHR - mu + 1.96 * np.sqrt(t.se ** 2 + tau ** 2))
    t["mu"], t["tau"] = mu, tau; out.append(t)
    nraw = int((t.fdr_raw < 0.05).sum()); ncal = int((t.fdr_cal < 0.05).sum()); ncal_pos = int(((t.fdr_cal < 0.05) & (t.z_cal > 0)).sum())
    summ.append([ix, mdl, len(nc), np.exp(mu), tau, int((nc.p < 0.05).sum()), len(t), nraw, ncal, ncal_pos, ncal / max(nraw, 1), t.z_cal.abs().median(),
                 t.loc[t.fdr_cal < 0.05, "outcome"].map(CHAP).value_counts().to_dict()])
C = pd.concat(out); C["title"] = C.outcome.map(TITLE); C["chapter"] = C.outcome.map(CHAP)
C.to_csv("../results/u04_calibrated.csv", index=False)
S = pd.DataFrame(summ, columns=["index", "model", "n_NC", "NC_mean_HR", "NC_tau", "NC_p<0.05", "n_outcomes", "sig_raw_FDR", "sig_calibrated_FDR", "sig_cal_positive",
                                "retained_fraction", "median_abs_z_cal", "calibrated_sig_by_chapter"])
S = S.sort_values(["model", "sig_calibrated_FDR"], ascending=[True, False]); S.to_csv("../results/u04_index_summary.csv", index=False)
pd.set_option("display.width", 250); print(S.drop(columns="calibrated_sig_by_chapter").round(3).to_string(index=False))
COMPOSITE = ["TyG", "TyG_BMI", "TyG_WC", "TyG_WHtR", "TyG_ABSI", "TG_HDL", "AIP", "METS_IR", "neg_eGDR", "neg_eGDR_noHTN", "CHG", "CTI"]
B = []
for o, g in C[(C.model == "M1")].groupby("outcome"):
    gi = g[g["index"].isin(COMPOSITE)]; top = gi.loc[gi.z_cal.idxmax()]
    sig = gi[(gi.fdr_cal < 0.05) & (gi.z_cal > 0)]["index"].tolist()
    B.append([o, TITLE.get(o, ""), CHAP.get(o, ""), int(g.events.max()), top["index"], top.HR, top.z_cal, top.fdr_cal, len(sig), ";".join(sig)])
B = pd.DataFrame(B, columns=["outcome", "title", "chapter", "events", "best_index", "best_HR", "best_z_cal", "best_fdr_cal", "n_indices_sig_cal", "indices_sig_cal"])
B.to_csv("../results/u04_best_index.csv", index=False)
print("\nbest index counts among outcomes with ≥1 calibrated-significant index:\n", B[B.n_indices_sig_cal > 0].best_index.value_counts().to_string())
print("outcomes with any calibrated-significant composite index:", int((B.n_indices_sig_cal > 0).sum()), "of", len(B))
TG = ["fo_E11", "fo_I10", "fo_I20", "fo_I21", "fo_I25", "fo_I50", "fo_I48", "fo_I63", "fo_I70", "fo_I73", "fo_I35", "fo_N18", "fo_K76", "fo_K80", "fo_M10", "fo_K85", "fo_G47",
      "fo_E66", "fo_E78", "fo_F03", "fo_G30", "fo_I26", "fo_I80", "ca_C18", "ca_C22", "ca_C25", "ca_C64", "ca_C54", "death", "cvdeath"]
T = C[(C.model == "M1") & C.outcome.isin(TG)].copy()
T["cell"] = T.apply(lambda r: f"{r.HR:.2f}" + ("‡" if r.fdr_cal < 0.05 else ("†" if r.fdr_raw < 0.05 else "")), axis=1)
W = T.pivot(index="outcome", columns="index", values="cell"); W.insert(0, "title", W.index.map(TITLE))
W.to_csv("../results/u04_targets.csv"); print("\n", W[["title", "TyG", "TyG_WC", "METS_IR", "neg_eGDR", "neg_eGDR_noHTN", "CHG", "c_lnGlu", "c_HbA1c", "ref_ApoB"]].to_string())
