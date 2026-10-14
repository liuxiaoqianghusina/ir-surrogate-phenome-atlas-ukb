"""u08：机制层 —— 基线指数与约 9 年后腹部/心脏 MRI 表型（肝脏 PDFF、胰腺 PDFF、内脏脂肪、腹部皮下脂肪、肌肉脂肪浸润、左室质量指数、射血分数）。
线性回归，表型取 SD 单位；调整年龄、性别、种族、TDI、吸烟、两次访视间隔。"""
import pandas as pd, numpy as np, re
import statsmodels.formula.api as smf
M = "/path/to/ukb_data"
h = pd.read_csv(f"{M}/ukb_assessment_center.csv", nrows=0).columns
use = ["eid"] + [c for c in ["p21090_i2", "p22435_i2", "p22436_i2", "p21088_i2", "p40061_i2", "p22407_i2", "p22408_i2", "p24105_i2", "p24103_i2", "p53_i2", "p50_i2", "p21002_i2"] if c in h]
A = pd.concat(pd.read_csv(f"{M}/ukb_assessment_center.csv", usecols=use, dtype=str, chunksize=40000), ignore_index=True)
for c in use[1:]:
    if c != "p53_i2": A[c] = pd.to_numeric(A[c], errors="coerce")
A["eid"] = A.eid.astype(int)
A["liver_pdff"] = A.p40061_i2.fillna(A.p22436_i2).fillna(A.p21088_i2); A["panc_pdff"] = A.p21090_i2; A["vat"] = A.p22407_i2; A["asat"] = A.p22408_i2
A["mfi"] = A.p22435_i2; bsa = np.sqrt(A.p50_i2 * A.p21002_i2 / 3600); A["lvmi"] = A.p24105_i2 / bsa; A["lvef"] = A.p24103_i2
B = pd.read_csv("../cache/u_base.csv")
IDEA = "/path/to/upstream"
d0 = pd.read_pickle(f"{IDEA}/ddw_ppdm/cache/cohort_full.pkl")[["eid", "d0"]]
X = B.merge(A[["eid", "p53_i2", "liver_pdff", "panc_pdff", "vat", "asat", "mfi", "lvmi", "lvef"]], on="eid").merge(d0, on="eid")
X["gap"] = (pd.to_datetime(X.p53_i2, errors="coerce") - X.d0).dt.days / 365.25
TR = {"liver_pdff": "Liver PDFF", "panc_pdff": "Pancreas PDFF", "vat": "VAT", "asat": "ASAT", "mfi": "Muscle fat infiltration", "lvmi": "LV mass index", "lvef": "LVEF"}
IDX = [c[2:] for c in B.columns if c.startswith("z_")]
rows = []
for tr in TR:
    d = X[X[tr].notna() & X.gap.notna()].copy()
    d[tr] = np.log(d[tr].clip(lower=0.1)) if tr in ("liver_pdff", "panc_pdff") else d[tr]
    z = (d[tr] - d[tr].mean()) / d[tr].std(); d = d[z.abs() < 5]; d["y"] = (d[tr] - d[tr].mean()) / d[tr].std()
    for ix in IDX:
        m = smf.ols(f"y ~ z_{ix} + age + sex + C(eth) + tdi + C(smk) + gap", d).fit()
        rows.append([TR[tr], ix, int(m.nobs), m.params[f"z_{ix}"], m.bse[f"z_{ix}"], m.pvalues[f"z_{ix}"], m.rsquared])
R = pd.DataFrame(rows, columns=["trait", "index", "n", "beta", "se", "p", "R2_model"]); R.to_csv("../results/u08_mri.csv", index=False)
pd.set_option("display.width", 250)
print(R.pivot(index="index", columns="trait", values="beta").round(3).to_string()); print(R.groupby("trait").n.max().to_dict())
