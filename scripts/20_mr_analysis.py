"""mr03：两样本孟德尔随机化 —— 暴露：甘油三酯（GLGC 2013）、BMI 和腰臀比调整 BMI（GIANT 2015）；结局：FinnGen R12 的 10 个终点。
单变量：IVW（乘法随机效应）、MR-Egger、加权中位数（自助法 SE），Cochran Q 和 I²。多变量 MR（IVW）：TG + BMI；TG + BMI + WHRadjBMI，并给出条件 F 统计量。
所有估计为每 1 SD 遗传预测暴露的 OR。输出 results/mr03_harmonised.csv、mr03_univariable.csv、mr03_mvmr.csv。用 python3 -P 运行。"""
import numpy as np, pandas as pd
from scipy import stats
ROOT = "/path/to/project"
OUT = {"I9_CHD": "Coronary heart disease", "I9_HEARTFAIL": "Heart failure", "I9_AF": "Atrial fibrillation", "I9_STR": "Stroke", "T2D": "Type 2 diabetes",
       "M13_GOUT": "Gout", "K11_CHOLELITH": "Gallstones", "NAFLD": "Non-alcoholic fatty liver disease", "N14_CHRONKIDNEYDIS": "Chronic kidney disease", "K11_ACUTPANC": "Acute pancreatitis"}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
I = pd.read_csv(f"{ROOT}/results/mr01_instruments.csv", dtype={"chr": str})
F = pd.read_csv(f"{ROOT}/results/mr02_finngen_lookup.csv", dtype={"chrom": str, "rsids": str})
FULL = {k: pd.read_pickle(f"{ROOT}/cache/mr_{k}.pkl").drop_duplicates("rsid").set_index("rsid") for k in ["TG", "BMI", "WHRadjBMI"]}


def align(ea, oa, eaf, ref, alt, af_alt):
    """返回 +1（暴露效应等位基因 = FinnGen alt）、−1（= ref）或 None（无法对齐/回文且频率接近 0.5）。"""
    ea, oa = str(ea).upper(), str(oa).upper()
    if len(ea) != 1 or len(oa) != 1 or len(ref) != 1 or len(alt) != 1: return None
    pal = {ea, oa} in ({"A", "T"}, {"C", "G"})
    if pal:
        if pd.isna(eaf) or abs(eaf - 0.5) < 0.08: return None
        return 1 if abs(eaf - af_alt) < abs(eaf - (1 - af_alt)) else -1
    if (ea, oa) == (alt, ref) or (COMP.get(ea), COMP.get(oa)) == (alt, ref): return 1
    if (ea, oa) == (ref, alt) or (COMP.get(ea), COMP.get(oa)) == (ref, alt): return -1
    return None


def outcome_row(ep, rsid, ea, oa):
    d = F[(F.endpoint == ep) & (F.query_rsid == rsid)]
    if d.empty: return None
    hit = d[d.rsids.fillna("").str.split(",").map(lambda s: rsid in s)]
    if hit.empty:   # 按等位基因匹配
        al = {str(ea).upper(), str(oa).upper()}; alc = {COMP.get(x, x) for x in al}
        hit = d[d.apply(lambda r: {r.ref, r.alt} in (al, alc), axis=1)]
    return None if hit.empty else hit.iloc[0]


# ---------------- 对齐 ----------------
H = []
for ep in OUT:
    for r in I.itertuples(index=False):
        o = outcome_row(ep, r.rsid, r.ea, r.oa)
        if o is None: continue
        s = align(r.ea, r.oa, r.eaf, o.ref, o.alt, float(o.af_alt))
        if s is None: continue
        H.append([r.exposure, ep, r.rsid, r.chr, r.pos38, r.ea, r.oa, r.eaf, r.beta, r.se, r.p, s * float(o.beta), float(o.sebeta), float(o.pval), float(o.af_alt)])
H = pd.DataFrame(H, columns=["exposure", "endpoint", "rsid", "chr", "pos38", "ea", "oa", "eaf", "b_exp", "se_exp", "p_exp", "b_out", "se_out", "p_out", "af_alt_finngen"])
H.to_csv(f"{ROOT}/results/mr03_harmonised.csv", index=False)


def ivw(bx, by, sy):
    w = 1 / sy ** 2; b = np.sum(w * bx * by) / np.sum(w * bx ** 2); se_f = np.sqrt(1 / np.sum(w * bx ** 2))
    k = len(bx); Q = np.sum(w * (by - b * bx) ** 2); phi = max(1, Q / (k - 1)) if k > 1 else 1
    return b, se_f * np.sqrt(phi), Q, k


def egger(bx, by, sy):
    s = np.sign(bx); bx, by = bx * s, by * s; w = 1 / sy ** 2
    X = np.column_stack([np.ones_like(bx), bx]); W = np.diag(w)
    XtWX = X.T @ W @ X; beta = np.linalg.solve(XtWX, X.T @ W @ by); res = by - X @ beta
    sig = max(1, np.sqrt(np.sum(w * res ** 2) / (len(bx) - 2))); V = np.linalg.inv(XtWX) * sig ** 2
    return beta[1], np.sqrt(V[1, 1]), beta[0], np.sqrt(V[0, 0])


def wmedian_est(b, w):
    o = np.argsort(b); b, w = b[o], w[o] / w.sum(); cs = np.cumsum(w) - 0.5 * w
    j = np.searchsorted(cs, 0.5); j = min(max(j, 1), len(b) - 1)
    return b[j - 1] + (b[j] - b[j - 1]) * (0.5 - cs[j - 1]) / (cs[j] - cs[j - 1])


def wmedian(bx, by, sx, sy, nboot=1000, seed=1):
    br = by / bx; se_r = np.sqrt(sy ** 2 / bx ** 2 + by ** 2 * sx ** 2 / bx ** 4); w = 1 / se_r ** 2
    est = wmedian_est(br, w); rng = np.random.default_rng(seed); bs = []
    for _ in range(nboot):
        bxb = rng.normal(bx, sx); byb = rng.normal(by, sy); bs.append(wmedian_est(byb / bxb, w))
    return est, np.std(bs)


U = []
for (ex, ep), d in H.groupby(["exposure", "endpoint"]):
    bx, by, sx, sy = d.b_exp.values, d.b_out.values, d.se_exp.values, d.se_out.values
    b, se, Q, k = ivw(bx, by, sy); qp = stats.chi2.sf(Q, k - 1); I2 = max(0, (Q - (k - 1)) / Q) if Q > 0 else 0
    eb, ese, ei, eise = egger(bx, by, sy); mb, mse = wmedian(bx, by, sx, sy)
    for m, bb, ss in [("IVW", b, se), ("MR-Egger", eb, ese), ("Weighted median", mb, mse)]:
        z = bb / ss
        U.append([ex, ep, OUT[ep], m, k, np.exp(bb), np.exp(bb - 1.96 * ss), np.exp(bb + 1.96 * ss), 2 * stats.norm.sf(abs(z)), Q, qp, I2,
                  ei if m == "MR-Egger" else np.nan, 2 * stats.norm.sf(abs(ei / eise)) if m == "MR-Egger" else np.nan])
U = pd.DataFrame(U, columns=["exposure", "endpoint", "outcome", "method", "nsnp", "OR", "lo", "hi", "p", "Q", "Q_p", "I2", "egger_intercept", "egger_intercept_p"])
U.to_csv(f"{ROOT}/results/mr03_univariable.csv", index=False)


# ---------------- 多变量 MR ----------------
def mv_matrix(exps):
    """所有暴露工具变量的并集（去掉 1 Mb 内重复，保留暴露 P 最小者），每个 SNP 取各暴露的效应（对齐到同一等位基因）。"""
    s = I[I.exposure.isin(exps)].sort_values("p"); keep = []; taken = {}
    for r in s.itertuples(index=False):
        if any(abs(r.pos38 - q) < 1_000_000 for q in taken.get(r.chr, [])): continue
        keep.append(r); taken.setdefault(r.chr, []).append(r.pos38)
    rows = []
    for r in keep:
        ea, oa = r.ea.upper(), r.oa.upper(); vals = []
        for e in exps:
            f = FULL[e]
            if r.rsid not in f.index: vals = None; break
            g = f.loc[r.rsid]; a1, a2 = g.ea.upper(), g.oa.upper()
            if (a1, a2) == (ea, oa) or (COMP.get(a1), COMP.get(a2)) == (ea, oa): sg = 1
            elif (a1, a2) == (oa, ea) or (COMP.get(a1), COMP.get(a2)) == (oa, ea): sg = -1
            else: vals = None; break
            vals += [sg * g.beta, g.se]
        if vals is not None: rows.append([r.rsid, r.chr, r.pos38, ea, oa, r.eaf] + vals)
    cols = ["rsid", "chr", "pos38", "ea", "oa", "eaf"] + sum([[f"b_{e}", f"se_{e}"] for e in exps], [])
    return pd.DataFrame(rows, columns=cols)


MV = []
for exps in [["TG", "BMI"], ["TG", "BMI", "WHRadjBMI"]]:
    M = mv_matrix(exps)
    # 条件 F（Sanderson 2021，假设暴露间 SNP 效应协方差为 0）
    condF = {}
    for e in exps:
        oth = [x for x in exps if x != e]; Y = M[f"b_{e}"].values; X = M[[f"b_{x}" for x in oth]].values
        dlt = np.linalg.lstsq(X, Y, rcond=None)[0]; resid = Y - X @ dlt
        var = M[f"se_{e}"].values ** 2 + sum((dlt[i] ** 2) * M[f"se_{x}"].values ** 2 for i, x in enumerate(oth))
        condF[e] = np.sum(resid ** 2 / var) / (len(M) - len(oth))
    for ep in OUT:
        rows = []
        for r in M.itertuples(index=False):
            o = outcome_row(ep, r.rsid, r.ea, r.oa)
            if o is None: continue
            s = align(r.ea, r.oa, r.eaf, o.ref, o.alt, float(o.af_alt))
            if s is None: continue
            rows.append([getattr(r, f"b_{e}") for e in exps] + [s * float(o.beta), float(o.sebeta)])
        A = np.array(rows); Bx, by, sy = A[:, :len(exps)], A[:, -2], A[:, -1]; w = 1 / sy ** 2
        XtWX = Bx.T @ (Bx * w[:, None]); beta = np.linalg.solve(XtWX, Bx.T @ (w * by)); res = by - Bx @ beta
        phi = max(1, np.sum(w * res ** 2) / (len(by) - len(exps))); V = np.linalg.inv(XtWX) * phi
        for i, e in enumerate(exps):
            se = np.sqrt(V[i, i]); MV.append(["+".join(exps), e, ep, OUT[ep], len(by), np.exp(beta[i]), np.exp(beta[i] - 1.96 * se), np.exp(beta[i] + 1.96 * se),
                                               2 * stats.norm.sf(abs(beta[i] / se)), condF[e]])
MV = pd.DataFrame(MV, columns=["model", "exposure", "endpoint", "outcome", "nsnp", "OR", "lo", "hi", "p", "conditional_F"])
MV.to_csv(f"{ROOT}/results/mr03_mvmr.csv", index=False)
pd.set_option("display.width", 250)
u = U[U.method == "IVW"].copy(); u["s"] = u.apply(lambda r: f"{r.OR:.2f} ({r.lo:.2f}-{r.hi:.2f}) p={r.p:.1g}", axis=1)
print(u.pivot(index="outcome", columns="exposure", values="s").to_string())
m = MV[MV.model == "TG+BMI+WHRadjBMI"].copy(); m["s"] = m.apply(lambda r: f"{r.OR:.2f} ({r.lo:.2f}-{r.hi:.2f}) p={r.p:.1g}", axis=1)
print(m.pivot(index="outcome", columns="exposure", values="s").to_string()); print(MV.groupby(["model", "exposure"]).conditional_F.first())
print(U[U.method != "IVW"].assign(s=lambda d: d.OR.round(2).astype(str) + " p=" + d.p.map(lambda v: f"{v:.1g}")).pivot_table(index=["outcome"], columns=["exposure", "method"], values="s", aggfunc="first").to_string())
