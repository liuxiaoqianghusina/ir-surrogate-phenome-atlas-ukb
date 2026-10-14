"""t02：构建分析表（基线无糖尿病者）—— 11 种 IR 替代指数、6 个成分、ApoB；各结局的随访时间、事件与基线患病标记。输出 cache/ir_atlas.csv（供 R）。"""
import pandas as pd, numpy as np
IDEA = "/path/to/upstream"
END = pd.Timestamp("2022-12-31"); SPECIAL = {1900, 1901, 1902, 1903, 1909, 2037}
b = pd.read_pickle(f"{IDEA}/ddw_ppdm/cache/cohort_full.pkl")[["eid", "d0", "p31", "birth", "bmi", "smoke", "tdi", "ethnic"]]
bs = pd.read_pickle(f"{IDEA}/cd_screening/cache/bs_lipgluc.pkl"); bs["eid"] = bs.eid.astype(int)
bc = pd.read_pickle(f"{IDEA}/diet_barrier/cache/biochem.pkl")[["eid", "p30710_i0"]]; bc["eid"] = bc.eid.astype(int)
ac = pd.read_pickle(f"{IDEA}/gut_hormone/cache/ac_cov.pkl"); ac["eid"] = ac.eid.astype(int)
ex = pd.read_pickle(f"{IDEA}/idea2_feasibility/cache/extra_cov.pkl")[["eid", "p6177_i0", "p6153_i0"]]; ex["eid"] = ex.eid.astype(int)
ae = pd.read_pickle("../cache/ac_extra.pkl").drop(columns=["p2443_i0"]); ae["eid"] = ae.eid.astype(int)
fo = pd.read_pickle("../cache/fo.pkl"); fo["eid"] = fo.eid.astype(int)
X = b.merge(bs, on="eid", how="left").merge(bc, on="eid", how="left").merge(ac, on="eid", how="left").merge(ex, on="eid", how="left").merge(ae, on="eid", how="left").merge(fo, on="eid", how="left")
num = lambda c: pd.to_numeric(X[c], errors="coerce")
X["age"] = (X.d0 - X.birth).dt.days / 365.25; X["sex"] = X.p31
X["white"] = X.ethnic.astype(str).str.contains("White|1").astype(int); X["smk"] = X.smoke.astype(str)
tg, glu, hdl, tc, hb = num("p30870_i0"), num("p30740_i0"), num("p30760_i0"), num("p30690_i0"), num("p30750_i0")
wc = num("p48_i0"); ht = num("p50_i0"); sbp = pd.concat([num("p4080_i0_a0"), num("p4080_i0_a1")], axis=1).mean(axis=1)
dbp = pd.concat([num("p4079_i0_a0"), num("p4079_i0_a1")], axis=1).mean(axis=1).fillna(pd.concat([num("p94_i0_a0"), num("p94_i0_a1")], axis=1).mean(axis=1))
has = lambda s, c: s.fillna("").astype(str).str.contains(rf"(?:^|\|){c}(?:\||$)")
bpmed = has(X.p6177_i0, 2) | has(X.p6153_i0, 2); dmmed = has(X.p6177_i0, 3) | has(X.p6153_i0, 3)
sr = X[[c for c in X.columns if c.startswith("p20002_i0_a")]].astype(str)
htn_sr = sr.isin(["1065", "1072"]).any(axis=1) | X.p6150_i0.fillna("").astype(str).str.contains(r"(?:^|\|)4(?:\||$)")
fod = lambda c: pd.to_datetime(X["fo_" + c], errors="coerce")
htn = (sbp >= 140) | (dbp >= 90) | bpmed | htn_sr | (fod("I10") <= X.d0)
TGmg, GLUmg, HDLmg, TCmg = tg * 88.57, glu * 18.0, hdl * 38.67, tc * 38.67
hb_pct = hb / 10.929 + 2.15
X["TyG"] = np.log(TGmg * GLUmg / 2)
X["TyG_BMI"] = X.TyG * X.bmi; X["TyG_WC"] = X.TyG * wc; X["TyG_WHtR"] = X.TyG * wc / ht
absi = (wc / 100) / (X.bmi ** (2 / 3) * (ht / 100) ** 0.5); X["TyG_ABSI"] = X.TyG * absi
X["TG_HDL"] = tg / hdl; X["AIP"] = np.log10(tg / hdl)
X["METS_IR"] = np.log(2 * GLUmg + TGmg) * X.bmi / np.log(HDLmg)
X["neg_eGDR"] = -(21.158 - 0.09 * wc - 3.407 * htn.astype(float) - 0.551 * hb_pct)      # 取负，使数值越高 = 越胰岛素抵抗
X["CHG"] = np.log(TCmg * GLUmg / (2 * HDLmg))
X["CTI"] = 0.412 * np.log(num("p30710_i0").clip(lower=0.1)) + X.TyG
X["c_lnTG"] = np.log(tg); X["c_lnGlu"] = np.log(glu); X["c_negHDL"] = -hdl; X["c_BMI"] = X.bmi; X["c_WC"] = wc; X["c_HbA1c"] = hb; X["ref_ApoB"] = num("p30640_i0")
IDX = ["TyG", "TyG_BMI", "TyG_WC", "TyG_WHtR", "TyG_ABSI", "TG_HDL", "AIP", "METS_IR", "neg_eGDR", "CHG", "CTI", "c_lnTG", "c_lnGlu", "c_negHDL", "c_BMI", "c_WC", "c_HbA1c", "ref_ApoB"]
dm_prev = pd.concat([fod(c) for c in ["E10", "E11", "E12", "E13", "E14"]], axis=1).min(axis=1) <= X.d0
dm_bad = pd.concat([pd.to_datetime(X["fo_" + c], errors="coerce").dt.year.isin(SPECIAL) for c in ["E10", "E11", "E12", "E13", "E14"]], axis=1).any(axis=1)
X = X[~(dm_prev | dm_bad | (num("p2443_i0") == 1) | (hb >= 48) | dmmed)]
X = X.dropna(subset=IDX + ["age", "sex", "tdi"]).copy()
for c in IDX: X["z_" + c] = (X[c] - X[c].mean()) / X[c].std()
death = pd.to_datetime(X.p40000_i0, errors="coerce"); stop = pd.concat([death, pd.Series(END, index=X.index)], axis=1).min(axis=1)
OUT = [c[3:] for c in X.columns if c.startswith("fo_")]
cols = {"eid": X.eid, "age": X.age, "sex": X.sex, "tdi": X.tdi, "white": X.white, "smk": X.smk}
for c in IDX: cols["z_" + c] = X["z_" + c]
for o in OUT:
    d = pd.to_datetime(X["fo_" + o], errors="coerce"); bad = d.dt.year.isin(SPECIAL); d = d.mask(bad)
    ev = d.notna() & (d > X.d0) & (d <= stop)
    cols["t_" + o] = np.where(ev, (d - X.d0).dt.days, (stop - X.d0).dt.days) / 365.25
    cols["e_" + o] = ev.astype(int); cols["p_" + o] = ((d.notna() & (d <= X.d0)) | bad).astype(int)
cvd = X.p40001_i0.fillna("").str.match(r"I(0\d|1\d|2\d|3\d|4\d|5[0-2]|6\d|7\d)")
cols["t_death"] = (stop - X.d0).dt.days / 365.25; cols["e_death"] = (death.notna() & (death <= END)).astype(int); cols["p_death"] = 0
cols["t_cvdeath"] = cols["t_death"]; cols["e_cvdeath"] = (cols["e_death"].astype(bool) & cvd).astype(int); cols["p_cvdeath"] = 0
D = pd.DataFrame(cols); D.to_csv("../cache/ir_atlas.csv", index=False)
print("n =", len(D), "| indices:", len(IDX), "| outcomes:", len(OUT) + 2)
print(D[["z_" + c for c in IDX[:11]]].corr().round(2).iloc[:, :6].to_string())
print({o: int(D["e_" + o][D["p_" + o] == 0].sum()) for o in OUT + ["death", "cvdeath"]})
