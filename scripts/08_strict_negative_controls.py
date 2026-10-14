"""u04b：严格阴性对照敏感性分析 —— 去掉与肥胖/机械负荷/便秘/皮肤皱褶等可能相关的对照（疝、痔、肛门病、膝/足畸形、皮炎、子宫/宫颈病变、眩晕、鼻窦炎），只保留眼、耳、口腔牙齿、上呼吸道感染/过敏、痤疮、色素/良性皮肤病变、乳腺良性、阑尾炎、头痛等。"""
import pandas as pd, numpy as np
from scipy import stats, optimize
R = pd.read_csv("../results/u03_phewas.csv")
STRICT = ["H52", "H61", "H00", "H01", "H02", "H04", "H10", "H11", "H16", "H18", "H43", "H65", "H66", "H93", "J01", "J02", "J03", "J30", "J31", "J33", "J34", "J35",
          "K02", "K04", "K08", "K12", "K13", "K35", "L20", "L57", "L60", "L70", "L72", "L80", "L82", "N60", "N43", "N47", "D22", "D23", "G43", "G44"]
ALL_NC = pd.read_csv("../results/u04_negative_controls.csv").outcome.tolist()
NC = [c for c in R.outcome.unique() if c[3:] in STRICT]
print("strict NCs:", len(NC))
def fit_null(b, se):
    nll = lambda th: -np.sum(stats.norm.logpdf(b, th[0], np.sqrt(se ** 2 + np.exp(th[1]) ** 2)))
    r = optimize.minimize(nll, [np.mean(b), np.log(max(np.std(b), 1e-3))], method="Nelder-Mead"); return r.x[0], np.exp(r.x[1])
def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = np.empty_like(p); n = len(p); c = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; q[o] = np.minimum(c, 1); return q
rows = []; tg = ["fo_E11", "fo_K76", "fo_M10", "fo_K80", "fo_I21", "fo_I25", "fo_I50", "fo_I48", "fo_I63", "fo_N18", "cvdeath", "death"]
for (ix, m), g in R.groupby(["index", "model"]):
    nc = g[g.outcome.isin(NC)]; mu, tau = fit_null(nc.logHR.values, nc.se.values)
    t = g[~g.outcome.isin(ALL_NC)].copy(); t["z"] = (t.logHR - mu) / np.sqrt(t.se ** 2 + tau ** 2); t["fdr"] = bh(2 * stats.norm.sf(t.z.abs()))
    sig = t[t.fdr < 0.05]
    rows.append([ix, m, np.exp(mu), tau, int((nc.p < 0.05).sum()), len(nc), int((t.p < 0.05).sum()), int(len(sig)), int((sig.z > 0).sum())] +
                [("‡" if (t.set_index("outcome").loc[o, "fdr"] < 0.05) else "·") if o in set(t.outcome) else "" for o in tg])
S = pd.DataFrame(rows, columns=["index", "model", "NC_mean_HR", "tau", "NC_p<0.05", "n_NC", "raw_p<0.05", "sig_cal", "sig_cal_pos"] + [o.replace("fo_", "") for o in tg])
S = S.sort_values(["model", "sig_cal"], ascending=[True, False]); S.to_csv("../results/u04b_strict_nc.csv", index=False)
pd.set_option("display.width", 250); print(S.round(3).to_string(index=False))
