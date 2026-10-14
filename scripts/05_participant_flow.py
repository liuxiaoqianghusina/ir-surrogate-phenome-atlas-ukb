"""u14：参与者流程（复用 u02 的数据合并与定义，只统计各步排除人数，不改动分析基表）。"""

import pandas as pd, numpy as np
import statsmodels.formula.api as smf
IDEA = "/path/to/upstream"
SPECIAL = {1900, 1901, 1902, 1903, 1909, 2037}
b = pd.read_pickle(f"{IDEA}/ddw_ppdm/cache/cohort_full.pkl")[["eid", "d0", "p31", "birth", "bmi", "smoke", "tdi", "ethnic"]]
bs = pd.read_pickle(f"{IDEA}/cd_screening/cache/bs_lipgluc.pkl"); bs["eid"] = bs.eid.astype(int)
bc = pd.read_pickle(f"{IDEA}/diet_barrier/cache/biochem.pkl")[["eid", "p30710_i0", "p30700_i0", "p30720_i0"]]; bc["eid"] = bc.eid.astype(int)
ac = pd.read_pickle(f"{IDEA}/gut_hormone/cache/ac_cov.pkl"); ac["eid"] = ac.eid.astype(int)
ex = pd.read_pickle(f"{IDEA}/idea2_feasibility/cache/extra_cov.pkl")[["eid", "p6177_i0", "p6153_i0"]]; ex["eid"] = ex.eid.astype(int)
ae = pd.read_pickle("../cache/ac_extra.pkl").drop(columns=["p2443_i0"]); ae["eid"] = ae.eid.astype(int)
fo = pd.read_pickle("../cache/fo.pkl"); fo["eid"] = fo.eid.astype(int)
X = b.merge(bs, on="eid", how="left").merge(bc, on="eid", how="left").merge(ac, on="eid", how="left").merge(ex, on="eid", how="left").merge(ae, on="eid", how="left").merge(fo, on="eid", how="left")
num = lambda c: pd.to_numeric(X[c], errors="coerce")
X["age"] = (X.d0 - X.birth).dt.days / 365.25; X["sex"] = X.p31
eth = num("p21000_i0")
X["eth"] = np.select([eth.isin([1, 1001, 1002, 1003]), eth.isin([3001, 3002, 3003]), eth.isin([4001, 4002, 4003, 4]), eth == 5], ["White", "South Asian", "Black", "Chinese"], "Other")
X["white"] = (X.eth == "White").astype(int)
X["smk"] = X.smoke.astype(str); alc = num("p1558_i0")
X["alc"] = np.select([alc == 1, alc.isin([2, 3]), alc.isin([4, 5]), alc == 6], ["daily", "weekly", "occasional", "never"], "missing")
X["centre"] = X.p54_i0.astype(str); X["fast"] = num("p74_i0")
dt = pd.to_datetime(X.p3166_i0_a0, errors="coerce"); X["hour"] = dt.dt.hour + dt.dt.minute / 60; X["month"] = num("p55_i0")
tg, glu, hdl, tc, hb, ldl = num("p30870_i0"), num("p30740_i0"), num("p30760_i0"), num("p30690_i0"), num("p30750_i0"), num("p30780_i0")
wc = num("p48_i0"); ht = num("p50_i0"); sbp = pd.concat([num("p4080_i0_a0"), num("p4080_i0_a1")], axis=1).mean(axis=1)
dbp = pd.concat([num("p4079_i0_a0"), num("p4079_i0_a1")], axis=1).mean(axis=1).fillna(pd.concat([num("p94_i0_a0"), num("p94_i0_a1")], axis=1).mean(axis=1))
has = lambda s, c: s.fillna("").astype(str).str.contains(rf"(?:^|\|){c}(?:\||$)")
bpmed = has(X.p6177_i0, 2) | has(X.p6153_i0, 2); lipmed = has(X.p6177_i0, 1) | has(X.p6153_i0, 1); dmmed = has(X.p6177_i0, 3) | has(X.p6153_i0, 3)
sr = X[[c for c in X.columns if c.startswith("p20002_i0_a")]].astype(str)
htn_sr = sr.isin(["1065", "1072"]).any(axis=1) | X.p6150_i0.fillna("").astype(str).str.contains(r"(?:^|\|)4(?:\||$)")
fod = lambda c: pd.to_datetime(X["fo_" + c], errors="coerce")
htn = (sbp >= 140) | (dbp >= 90) | bpmed | htn_sr | (fod("I10") <= X.d0)
cvd_prev = pd.concat([fod(c) for c in ["I21", "I25", "I63", "I50"]], axis=1).min(axis=1) <= X.d0
TGmg, GLUmg, HDLmg, TCmg = tg * 88.57, glu * 18.0, hdl * 38.67, tc * 38.67
hb_pct = hb / 10.929 + 2.15
I = pd.DataFrame(index=X.index)
I["TyG"] = np.log(TGmg * GLUmg / 2)
I["TyG_BMI"] = I.TyG * X.bmi; I["TyG_WC"] = I.TyG * wc; I["TyG_WHtR"] = I.TyG * wc / ht
absi = (wc / 100) / (X.bmi ** (2 / 3) * (ht / 100) ** 0.5); I["TyG_ABSI"] = I.TyG * absi
I["TG_HDL"] = tg / hdl; I["AIP"] = np.log10(tg / hdl)
I["METS_IR"] = np.log(2 * GLUmg + TGmg) * X.bmi / np.log(HDLmg)
I["neg_eGDR"] = -(21.158 - 0.09 * wc - 3.407 * htn.astype(float) - 0.551 * hb_pct)
I["neg_eGDR_noHTN"] = -(21.158 - 0.09 * wc - 0.551 * hb_pct)
I["CHG"] = np.log(TCmg * GLUmg / (2 * HDLmg))
I["CTI"] = 0.412 * np.log(num("p30710_i0").clip(lower=0.1)) + I.TyG
I["c_lnTG"] = np.log(tg); I["c_lnGlu"] = np.log(glu); I["c_negHDL"] = -hdl; I["c_BMI"] = X.bmi; I["c_WC"] = wc; I["c_HbA1c"] = hb; I["ref_ApoB"] = num("p30640_i0")
dm_prev = pd.concat([fod(c) for c in ["E10", "E11", "E12", "E13", "E14"]], axis=1).min(axis=1) <= X.d0
dm_bad = pd.concat([pd.to_datetime(X["fo_" + c], errors="coerce").dt.year.isin(SPECIAL) for c in ["E10", "E11", "E12", "E13", "E14"]], axis=1).any(axis=1)

n0 = len(X)
dm = dm_prev | dm_bad | (num("p2443_i0") == 1) | (hb >= 48) | dmmed
n_dm = int(dm.sum())
core = ["TyG", "TyG_BMI", "TyG_WC", "TyG_WHtR", "TyG_ABSI", "TG_HDL", "AIP", "METS_IR", "neg_eGDR", "CHG", "CTI", "c_lnTG", "c_lnGlu", "c_negHDL", "c_BMI", "c_WC", "c_HbA1c", "ref_ApoB"]
miss = ~(I[core].notna().all(axis=1) & X[["age", "sex", "tdi"]].notna().all(axis=1))
n_miss = int((~dm & miss).sum()); n_final = int((~dm & ~miss).sum())
# 缺失来源（在无糖尿病者中）
nd = ~dm
src = {"甘油三酯/血糖/HDL/总胆固醇": int((nd & (tg.isna() | glu.isna() | hdl.isna() | tc.isna())).sum()), "HbA1c": int((nd & hb.isna()).sum()),
       "CRP": int((nd & num("p30710_i0").isna()).sum()), "ApoB": int((nd & num("p30640_i0").isna()).sum()), "腰围/身高/BMI": int((nd & (wc.isna() | ht.isna() | X.bmi.isna())).sum()),
       "TDI": int((nd & X.tdi.isna()).sum())}
import pandas as pd
pd.DataFrame([["UK Biobank 参与者", n0], ["基线糖尿病（排除）", n_dm], ["无糖尿病但指数成分或核心协变量缺失（排除）", n_miss], ["分析人群", n_final]] +
             [[f"缺失来源（无糖尿病者中，可重叠）：{k}", v] for k, v in src.items()], columns=["step", "n"]).to_csv("../results/u14_flow.csv", index=False)
print(open("../results/u14_flow.csv", encoding="utf-8").read())
