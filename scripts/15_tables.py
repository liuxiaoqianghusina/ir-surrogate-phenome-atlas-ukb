"""u13：主表 —— 表 1 基线特征（全体及按 TyG 四分位）；表 2 指数评分卡（构念效度、特异性、调整 BMI 后特异性、严格阴性对照、T2D 与心血管增量 C）；
另输出随访与关键结局事件数（u13_followup.csv），供正文引用。"""
import pandas as pd, numpy as np
RES = "../results"
b = pd.read_csv("../cache/u_base.csv"); cen = pd.read_csv("../cache/phe_censor.csv"); b = b.merge(cen, on="eid")
b["fu"] = b.cday / 365.25
# TyG 原值（四分位用）：由 TG 和血糖（mmol/L → mg/dL）重算
b["tyg"] = np.log(b.tg * 88.57 * b.glu * 18.0 / 2); b["q"] = pd.qcut(b.tyg, 4, labels=["Q1", "Q2", "Q3", "Q4"])
G = [("全体", b)] + [(f"TyG {q}", b[b.q == q]) for q in ["Q1", "Q2", "Q3", "Q4"]]
mean_sd = lambda s: f"{s.mean():.1f} ± {s.std():.1f}"
med_iqr = lambda s, d=2: f"{s.median():.{d}f} ({s.quantile(.25):.{d}f}–{s.quantile(.75):.{d}f})"
pct = lambda m: f"{int(m.sum()):,} ({m.mean() * 100:.1f})"
rows = [("人数", lambda d: f"{len(d):,}"),
        ("TyG 范围", lambda d: f"{d.tyg.min():.2f}–{d.tyg.max():.2f}"),
        ("年龄，岁", lambda d: mean_sd(d.age)),
        ("女性，n (%)", lambda d: pct(d.sex == 0)),
        ("种族，n (%)", None),
        ("　白人", lambda d: pct(d.eth == "White")),
        ("　南亚裔", lambda d: pct(d.eth == "South Asian")),
        ("　黑人", lambda d: pct(d.eth == "Black")),
        ("　其他", lambda d: pct(d.eth.isin(["Chinese", "Other"]))),
        ("Townsend 剥夺指数", lambda d: med_iqr(d.tdi, 1)),
        ("吸烟，n (%)", None),
        ("　从不", lambda d: pct(d.smk == 0)),
        ("　既往", lambda d: pct(d.smk == 1)),
        ("　当前", lambda d: pct(d.smk == 2)),
        ("每天饮酒，n (%)", lambda d: pct(d.alc == "daily")),
        ("BMI，kg/m²", lambda d: mean_sd(d.bmi)),
        ("腰围，cm", lambda d: mean_sd(d.wc)),
        ("收缩压，mmHg", lambda d: mean_sd(d.sbp.dropna())),
        ("高血压，n (%)", lambda d: pct(d.htn == 1)),
        ("降压药，n (%)", lambda d: pct(d.bpmed == 1)),
        ("降脂药，n (%)", lambda d: pct(d.lipmed == 1)),
        ("既往心血管病，n (%)", lambda d: pct(d.cvd_prev == 1)),
        ("甘油三酯，mmol/L", lambda d: med_iqr(d.tg)),
        ("血糖，mmol/L", lambda d: med_iqr(d.glu)),
        ("HbA1c，mmol/mol", lambda d: mean_sd(d.hba1c)),
        ("总胆固醇，mmol/L", lambda d: mean_sd(d.tc)),
        ("HDL-C，mmol/L", lambda d: mean_sd(d.hdl)),
        ("LDL-C，mmol/L", lambda d: mean_sd(d.ldl.dropna())),
        ("ApoB，g/L", lambda d: mean_sd(d.apob)),
        ("距末次进食时间，h", lambda d: med_iqr(d.fast.dropna(), 0)),
        ("空腹 ≥8 h，n (%)", lambda d: pct(d.fast >= 8)),
        ("随访时间，年", lambda d: med_iqr(d.fu, 1))]
T1 = pd.DataFrame([[lab] + ([f(d) for _, d in G] if f else [""] * len(G)) for lab, f in rows], columns=["特征"] + [g for g, _ in G])
T1.to_csv(f"{RES}/u13_table1.csv", index=False)
# 随访与关键结局事件数（与 u03 一致：排除基线已患）
P = pd.concat([pd.read_csv(f, usecols=lambda c: c in ("eid", "fo_E11", "fo_I21", "fo_I25", "fo_I50", "fo_I48", "fo_I63", "fo_N18", "fo_K76", "fo_M10", "fo_K80", "fo_K85"))
               for f in sorted(__import__("glob").glob("../cache/phe_days_part*.csv"))], axis=1)
P = P.loc[:, ~P.columns.duplicated()].set_index("eid").reindex(b.eid)
fu = [["参与者", len(b)], ["中位随访（年）", round(b.fu.median(), 1)], ["随访四分位间距", f"{b.fu.quantile(.25):.1f}–{b.fu.quantile(.75):.1f}"],
      ["总人年", int(b.fu.sum())], ["全因死亡", int(b.death_day.notna().sum())], ["心血管死亡", int(((b.cvdeath == 1) & b.death_day.notna()).sum())]]
for c in P.columns:
    v = P[c].values; fu.append([c, int(((v > 0) & (v <= b.cday.values)).sum())])
pd.DataFrame(fu, columns=["item", "value"]).to_csv(f"{RES}/u13_followup.csv", index=False)
# ---------------- 表 2：指数评分卡 ----------------
LAB = {"TyG": "TyG", "TyG_BMI": "TyG-BMI", "TyG_WC": "TyG-WC", "TyG_WHtR": "TyG-WHtR", "TyG_ABSI": "TyG-ABSI", "TG_HDL": "TG/HDL-C", "AIP": "AIP", "METS_IR": "METS-IR",
       "neg_eGDR": "eGDR†", "neg_eGDR_noHTN": "eGDR（去高血压项）†", "CHG": "CHG", "CTI": "CTI", "c_lnTG": "甘油三酯", "c_lnGlu": "血糖", "c_negHDL": "HDL-C†",
       "c_BMI": "BMI", "c_WC": "腰围", "c_HbA1c": "HbA1c", "ref_ApoB": "ApoB"}
COMP = {"TyG": "TG、血糖", "TyG_BMI": "TG、血糖、BMI", "TyG_WC": "TG、血糖、腰围", "TyG_WHtR": "TG、血糖、腰围、身高", "TyG_ABSI": "TG、血糖、腰围、BMI、身高",
        "TG_HDL": "TG、HDL-C", "AIP": "TG、HDL-C", "METS_IR": "TG、血糖、HDL-C、BMI", "neg_eGDR": "腰围、高血压、HbA1c", "neg_eGDR_noHTN": "腰围、HbA1c",
        "CHG": "总胆固醇、HDL-C、血糖", "CTI": "TG、血糖、CRP", "c_lnTG": "—", "c_lnGlu": "—", "c_negHDL": "—", "c_BMI": "—", "c_WC": "—", "c_HbA1c": "—", "ref_ApoB": "—"}
S = pd.read_csv(f"{RES}/u04_index_summary.csv"); S1 = S[S.model == "M1"].set_index("index"); S2 = S[S.model == "M2"].set_index("index")
SN = pd.read_csv(f"{RES}/u04b_strict_nc.csv"); SN1 = SN[SN.model == "M1"].set_index("index")
H = pd.read_csv(f"{RES}/v01_nhanes_homa.csv"); H = H[H.group == "All"].set_index("index")
P6 = pd.read_csv(f"{RES}/u06_prediction.csv"); P6 = P6[P6.base == "B"]
dc = lambda ix, o: (lambda r: f"{r.dC.iloc[0] * 1000:.1f}" if len(r) else "—")(P6[(P6["index"] == ix) & (P6.outcome == o)])
order = ["TyG", "TG_HDL", "AIP", "CHG", "CTI", "TyG_BMI", "TyG_WC", "TyG_WHtR", "TyG_ABSI", "METS_IR", "neg_eGDR", "neg_eGDR_noHTN", "c_lnTG", "c_lnGlu", "c_negHDL", "c_BMI", "c_WC", "c_HbA1c", "ref_ApoB"]
T2 = []
for ix in order:
    s1 = S1.loc[ix]
    T2.append([LAB[ix], COMP[ix], f"{H.loc[ix, 'spearman_r_HOMA']:.2f}" if ix in H.index else "—", f"{s1.NC_mean_HR:.3f}", f"{int(s1['NC_p<0.05'])}/60",
               int(s1.sig_raw_FDR), f"{int(s1.sig_calibrated_FDR)} ({int(s1.sig_cal_positive)})", int(SN1.loc[ix, "sig_cal"]),
               int(S2.loc[ix, "sig_calibrated_FDR"]) if ix in S2.index else "—", dc(ix, "T2D"), dc(ix, "MASLD"), dc(ix, "MACE"), dc(ix, "HF")])
T2 = pd.DataFrame(T2, columns=["指数", "组成", "与 HOMA-IR 的 r（NHANES）", "阴性对照背景 HR", "阴性对照 P<0.05", "原始 FDR 显著", "校准后显著（正向）",
                              "校准后显著（严格阴性对照）", "校准后显著（加调 BMI）", "ΔC×1000：T2D", "ΔC×1000：MASLD", "ΔC×1000：MACE", "ΔC×1000：心衰"])
T2.to_csv(f"{RES}/u13_table2.csv", index=False)
pd.set_option("display.width", 250); print(T1.to_string(index=False)); print(T2.to_string(index=False)); print(pd.DataFrame(fu).to_string())
