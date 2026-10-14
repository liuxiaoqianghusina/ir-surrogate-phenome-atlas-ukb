"""v01：NHANES 2005–2018 外部锚定 —— 空腹亚样本，成年人（≥20 岁）无糖尿病者。
(1) 各指数与 HOMA-IR 的相关（HOMA-IR 在各周期内做秩逆正态转换，以消除胰岛素检测方法变化），总体及按种族/性别；
(2) 识别 HOMA-IR 最高四分位（周期内）的 AUC；
(3) 与全因及心血管死亡的关联（链接死亡数据至 2019 年），以及调整 HOMA-IR 后的残余关联。
用 python3 -P 运行（读取下载文件）。"""
import glob, os, numpy as np, pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
D = "/path/to/project/external/nhanes"
OUT = "/path/to/project/results"
CY = {"D": 2005, "E": 2007, "F": 2009, "G": 2011, "H": 2013, "I": 2015, "J": 2017}
rd = lambda f: pd.read_sas(f"{D}/{f}.xpt", format="xport", encoding="latin1") if os.path.exists(f"{D}/{f}.xpt") else None
rows = []
for s, y in CY.items():
    demo = rd(f"DEMO_{s}"); glu = rd(f"GLU_{s}")
    if demo is None or glu is None: print("skip cycle", s); continue
    X = demo[["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1"] + (["RIDRETH3"] if "RIDRETH3" in demo else [])].merge(glu, on="SEQN")
    ins = rd(f"INS_{s}")
    if ins is not None: X = X.merge(ins[["SEQN", "LBXIN"]], on="SEQN", how="left", suffixes=("_g", ""))
    if "LBXIN" not in X and "LBXIN_g" in X: X["LBXIN"] = X["LBXIN_g"]
    for f, cols in [("TRIGLY", ["LBXTR"]), ("HDL", ["LBDHDD"]), ("TCHOL", ["LBXTC"]), ("GHB", ["LBXGH"]), ("BMX", ["BMXBMI", "BMXWAIST", "BMXHT"]),
                    ("BPQ", ["BPQ020", "BPQ040A"]), ("DIQ", ["DIQ010", "DIQ050", "DIQ070"])]:
        t = rd(f"{f}_{s}")
        if t is not None: X = X.merge(t[["SEQN"] + [c for c in cols if c in t]], on="SEQN", how="left")
    bp = rd(f"BPX_{s}")
    if bp is not None:
        X = X.merge(pd.DataFrame({"SEQN": bp.SEQN, "sbp": bp[[c for c in bp if c.startswith("BPXSY")]].mean(axis=1), "dbp": bp[[c for c in bp if c.startswith("BPXDI")]].replace(0, np.nan).mean(axis=1)}), on="SEQN", how="left")
    crp = rd(f"CRP_{s}"); hs = rd(f"HSCRP_{s}")
    if crp is not None: X = X.merge(pd.DataFrame({"SEQN": crp.SEQN, "crp": crp.LBXCRP * 10}), on="SEQN", how="left")
    elif hs is not None: X = X.merge(pd.DataFrame({"SEQN": hs.SEQN, "crp": hs.LBXHSCRP}), on="SEQN", how="left")
    mf = glob.glob(f"{D}/NHANES_{y}_{y + 1}_MORT_2019_PUBLIC.dat")
    if mf:
        m = pd.read_fwf(mf[0], colspecs=[(0, 6), (14, 15), (15, 16), (16, 19), (42, 45), (45, 48)], names=["SEQN", "elig", "mortstat", "ucod", "permth_int", "permth_exm"], dtype=str)
        for c in ["SEQN", "elig", "mortstat", "permth_exm"]: m[c] = pd.to_numeric(m[c], errors="coerce")
        X = X.merge(m, on="SEQN", how="left")
    X["cycle"] = s; rows.append(X); print(s, len(X), flush=True)
N = pd.concat(rows, ignore_index=True)
N = N[(N.RIDAGEYR >= 20) & N.LBXIN.notna() & N.LBXGLU.notna() & N.LBXTR.notna() & N.LBDHDD.notna()]
dm = (N.DIQ010 == 1) | (N.DIQ050 == 1) | (N.DIQ070 == 1) | (N.LBXGH >= 6.5) | (N.LBXGLU >= 126)
N = N[~dm].copy()
N["HOMA"] = N.LBXGLU * N.LBXIN / 405
N["HOMA_rin"] = N.groupby("cycle").HOMA.transform(lambda s: stats.norm.ppf((s.rank() - 0.5) / s.notna().sum()))
N["HOMA_q4"] = N.groupby("cycle").HOMA.transform(lambda s: (s >= s.quantile(0.75)).astype(int))
htn = (N.sbp >= 140) | (N.dbp >= 90) | (N.BPQ040A == 1) | (N.BPQ020 == 1)
I = pd.DataFrame(index=N.index)
I["TyG"] = np.log(N.LBXTR * N.LBXGLU / 2)
I["TyG_BMI"] = I.TyG * N.BMXBMI; I["TyG_WC"] = I.TyG * N.BMXWAIST; I["TyG_WHtR"] = I.TyG * N.BMXWAIST / N.BMXHT
I["TyG_ABSI"] = I.TyG * ((N.BMXWAIST / 100) / (N.BMXBMI ** (2 / 3) * (N.BMXHT / 100) ** 0.5))
I["TG_HDL"] = N.LBXTR / N.LBDHDD; I["AIP"] = np.log10((N.LBXTR / 88.57) / (N.LBDHDD / 38.67))
I["METS_IR"] = np.log(2 * N.LBXGLU + N.LBXTR) * N.BMXBMI / np.log(N.LBDHDD)
I["neg_eGDR"] = -(21.158 - 0.09 * N.BMXWAIST - 3.407 * htn.astype(float) - 0.551 * N.LBXGH)
I["neg_eGDR_noHTN"] = -(21.158 - 0.09 * N.BMXWAIST - 0.551 * N.LBXGH)
I["CHG"] = np.log(N.LBXTC * N.LBXGLU / (2 * N.LBDHDD))
I["CTI"] = 0.412 * np.log(N.crp.clip(lower=0.1)) + I.TyG if "crp" in N else np.nan
I["c_lnTG"] = np.log(N.LBXTR); I["c_lnGlu"] = np.log(N.LBXGLU); I["c_negHDL"] = -N.LBDHDD; I["c_BMI"] = N.BMXBMI; I["c_WC"] = N.BMXWAIST; I["c_HbA1c"] = N.LBXGH
eth = N.RIDRETH1.map({1: "Mexican American", 2: "Other Hispanic", 3: "NH White", 4: "NH Black", 5: "Other/Asian"})
print("NHANES analytic n =", len(N), eth.value_counts().to_dict())
res = []
for ix in I.columns:
    ok = I[ix].notna() & N.HOMA_rin.notna()
    for grp, m in [("All", ok), ("NH White", ok & (eth == "NH White")), ("NH Black", ok & (eth == "NH Black")), ("Mexican American", ok & (eth == "Mexican American")),
                   ("Men", ok & (N.RIAGENDR == 1)), ("Women", ok & (N.RIAGENDR == 2))]:
        if m.sum() < 200: continue
        r = stats.spearmanr(I.loc[m, ix], N.loc[m, "HOMA_rin"]).statistic
        y = N.loc[m, "HOMA_q4"].values; sc = I.loc[m, ix].values
        auc = stats.mannwhitneyu(sc[y == 1], sc[y == 0]).statistic / ((y == 1).sum() * (y == 0).sum())
        res.append([ix, grp, int(m.sum()), r, auc])
R = pd.DataFrame(res, columns=["index", "group", "n", "spearman_r_HOMA", "AUC_HOMA_Q4"]); R.to_csv(f"{OUT}/v01_nhanes_homa.csv", index=False)
print(R[R.group == "All"].sort_values("spearman_r_HOMA", ascending=False).round(3).to_string(index=False))
print(R.pivot(index="index", columns="group", values="spearman_r_HOMA").round(2).to_string())
# 死亡
if "mortstat" in N:
    M_ = N[(N.elig == 1) & N.permth_exm.notna()].copy(); M_["t"] = M_.permth_exm / 12; M_["death"] = (M_.mortstat == 1).astype(int)
    M_["cvd"] = ((M_.mortstat == 1) & M_.ucod.astype(str).str.strip().isin(["001", "005"])).astype(int)
    pd.DataFrame([["n", len(M_)], ["median_followup_y", round(M_.t.median(), 1)], ["q1", round(M_.t.quantile(.25), 1)], ["q3", round(M_.t.quantile(.75), 1)], ["deaths", int(M_.death.sum())], ["cv_deaths", int(M_.cvd.sum())]],
                 columns=["item", "value"]).to_csv(f"{OUT}/v01_nhanes_followup.csv", index=False)
    M_["male"] = (M_.RIAGENDR == 1).astype(int); M_["black"] = (M_.RIDRETH1 == 4).astype(int); M_["hisp"] = M_.RIDRETH1.isin([1, 2]).astype(int)
    out = []
    for ix in list(I.columns) + ["HOMA_rin"]:
        z = (I[ix] if ix != "HOMA_rin" else N.HOMA_rin).loc[M_.index]; M_["z"] = (z - z.mean()) / z.std()
        for ev in ["death", "cvd"]:
            for adj in ["base", "+HOMA"]:
                if adj == "+HOMA" and ix == "HOMA_rin": continue
                cols = ["t", ev, "z", "RIDAGEYR", "male", "black", "hisp"] + (["HOMA_rin"] if adj == "+HOMA" else [])
                d = M_[cols].dropna()
                s = CoxPHFitter(penalizer=1e-6).fit(d, "t", ev).summary.loc["z"]
                out.append([ix, ev, adj, len(d), int(d[ev].sum()), s["exp(coef)"], s["exp(coef) lower 95%"], s["exp(coef) upper 95%"], s["p"]])
    O = pd.DataFrame(out, columns=["index", "outcome", "model", "n", "events", "HR", "lo", "hi", "p"]); O.to_csv(f"{OUT}/v01_nhanes_mortality.csv", index=False)
    print(O.pivot_table(index="index", columns=["outcome", "model"], values="HR").round(2).to_string())
