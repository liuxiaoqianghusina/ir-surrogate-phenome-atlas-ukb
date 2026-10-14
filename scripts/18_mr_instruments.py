"""mr01：孟德尔随机化工具变量 —— 甘油三酯（GLGC 2013）、BMI（GIANT 2015，欧洲人群）、腰臀比调整 BMI（GIANT 2015，欧洲人群合并）、空腹血糖（MAGIC 2021，欧洲人群，单位 mmol/L）。
取 P<5×10⁻⁸ 的 SNP，经 Ensembl REST 得到 GRCh38 坐标，按距离贪心聚类（保留每 1 Mb 内最显著者），去掉 MHC（chr6:25–34 Mb）。
甘油三酯先用 hg19 坐标预聚类，只对候选主变异查询 Ensembl。输出 results/mr01_instruments.csv 和 cache/mr01_ensembl.csv。用 python3 -P 运行。"""
import json, time, subprocess, os
import pandas as pd, numpy as np
ROOT = "/path/to/project"
G = f"{ROOT}/external/gwas"
PTH, WIN = 5e-8, 1_000_000
tg = pd.read_csv(f"{G}/glgc2013_tg/jointGwasMc_TG.txt.gz", sep="\t", dtype={"rsid": str})
tg = pd.DataFrame({"rsid": tg.rsid, "ea": tg.A1.str.upper(), "oa": tg.A2.str.upper(), "eaf": tg["Freq.A1.1000G.EUR"], "beta": tg.beta, "se": tg.se, "p": tg["P-value"], "n": tg.N,
                   "chr19": tg.SNP_hg19.str.replace("chr", "").str.split(":").str[0], "pos19": pd.to_numeric(tg.SNP_hg19.str.split(":").str[1], errors="coerce")})
bmi = pd.read_csv(f"{G}/giant2015_bmi/SNP_gwas_mc_merge_nogc.tbl.uniq.gz", sep="\t")
bmi = pd.DataFrame({"rsid": bmi.SNP, "ea": bmi.A1.str.upper(), "oa": bmi.A2.str.upper(), "eaf": bmi["Freq1.Hapmap"], "beta": bmi.b, "se": bmi.se, "p": bmi.p, "n": bmi.N})
whr = pd.read_csv(f"{G}/giant2015_whradjbmi/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt", sep="\t")
whr = pd.DataFrame({"rsid": whr.MarkerName, "ea": whr.Allele1.str.upper(), "oa": whr.Allele2.str.upper(), "eaf": whr.FreqAllele1HapMapCEU, "beta": whr.b, "se": whr.se, "p": whr.p, "n": whr.N})
EXP = {"TG": tg, "BMI": bmi, "WHRadjBMI": whr}
for k, d in EXP.items(): d.to_pickle(f"{ROOT}/cache/mr_{k}.pkl")   # 全部汇总统计（供多变量 MR 查询）


def clump(d, chrcol, poscol):
    d = d.sort_values("p"); keep = []; taken = {}
    for r in d.itertuples(index=False):
        c, p = getattr(r, chrcol), getattr(r, poscol)
        if any(abs(p - q) < WIN for q in taken.get(c, [])): continue
        keep.append(r.rsid); taken.setdefault(c, []).append(p)
    return keep


def ensembl(ids):
    """rsID → GRCh38 坐标（主染色体）。"""
    out = {}
    for i in range(0, len(ids), 200):
        body = json.dumps({"ids": ids[i:i + 200]})
        for a in range(5):
            r = subprocess.run(["curl", "-s", "--max-time", "120", "-X", "POST", "https://rest.ensembl.org/variation/homo_sapiens", "-H", "Content-Type: application/json",
                                "-H", "Accept: application/json", "-d", body], capture_output=True, text=True)
            try:
                js = json.loads(r.stdout); break
            except json.JSONDecodeError:
                time.sleep(2 + 3 * a); js = {}
        for rid, v in js.items():
            for m in v.get("mappings", []):
                if m.get("assembly_name") == "GRCh38" and m.get("seq_region_name") in [str(c) for c in range(1, 23)] + ["X"]:
                    out[rid] = (m["seq_region_name"], int(m["start"]), m.get("allele_string", "")); break
        time.sleep(0.15)
    return out


# 候选 SNP：TG 先按 hg19 预聚类（窗口 500 kb，较宽松），BMI 和 WHRadjBMI 取全部显著 SNP
t0 = tg[(tg.p < PTH) & tg.pos19.notna()].sort_values("p"); keep = []; taken = {}
for r in t0.itertuples(index=False):
    if any(abs(r.pos19 - q) < 500_000 for q in taken.get(r.chr19, [])): continue
    keep.append(r.rsid); taken.setdefault(r.chr19, []).append(r.pos19)
cand = {"TG": keep, "BMI": bmi[bmi.p < PTH].rsid.tolist(), "WHRadjBMI": whr[whr.p < PTH].rsid.tolist()}
# 空腹血糖（约 3000 万行）：分块读取，只保留 P<5e-8 的 SNP 和其他暴露的候选 SNP（多变量 MR 要查它们的血糖效应）
need = set(sum(cand.values(), []))
parts = []
for ch in pd.read_csv(f"{G}/magic2021_fg/MAGIC1000G_FG_EUR.tsv.gz", sep="\t", dtype={"variant": str, "chromosome": str}, chunksize=2_000_000):
    parts.append(ch[(ch.p_value < PTH) | ch.variant.isin(need)])
fg = pd.concat(parts)
fg = pd.DataFrame({"rsid": fg.variant, "ea": fg.effect_allele.str.upper(), "oa": fg.other_allele.str.upper(), "eaf": fg.effect_allele_frequency, "beta": fg.beta, "se": fg.standard_error,
                   "p": fg.p_value, "n": fg.sample_size, "chr19": fg.chromosome, "pos19": fg.base_pair_location})
fg = fg[fg.rsid.str.startswith("rs", na=False)]
# 质控：小样本的罕见变异报告的 P 与 β/SE 不符（n≈1,820，β 达 1–2 mmol/L），只保留样本量 ≥ 最大值一半、MAF ≥1%、且 β/SE 也达到全基因组显著的变异作为工具变量
qc = (fg.n >= 0.5 * fg.n.max()) & fg.eaf.between(0.01, 0.99) & ((fg.beta / fg.se) ** 2 >= 29.72)
fg.loc[~qc, "p"] = 1.0
fg.drop(columns=["chr19", "pos19"]).to_pickle(f"{ROOT}/cache/mr_FG.pkl"); EXP["FG"] = fg
t0 = fg[fg.p < PTH].sort_values("p"); keep = []; taken = {}
for r in t0.itertuples(index=False):
    if any(abs(r.pos19 - q) < 500_000 for q in taken.get(r.chr19, [])): continue
    keep.append(r.rsid); taken.setdefault(r.chr19, []).append(r.pos19)
cand["FG"] = keep
ids = sorted(set(sum(cand.values(), [])))
cache = f"{ROOT}/cache/mr01_ensembl.csv"
known = pd.read_csv(cache, dtype={"chr": str}) if os.path.exists(cache) else pd.DataFrame(columns=["rsid", "chr", "pos38", "alleles"])
todo = [i for i in ids if i not in set(known.rsid)]
print({k: len(v) for k, v in cand.items()}, "Ensembl queries:", len(todo), flush=True)
if todo:
    mp = ensembl(todo)
    known = pd.concat([known, pd.DataFrame([[k, v[0], v[1], v[2]] for k, v in mp.items()], columns=known.columns)], ignore_index=True)
    known.to_csv(cache, index=False)
pos = known.drop_duplicates("rsid").set_index("rsid")
rows = []
for k, d in EXP.items():
    s = d[d.rsid.isin(cand[k])].copy(); s = s[s.rsid.isin(pos.index)]
    s["chr"] = s.rsid.map(pos.chr).astype(str); s["pos38"] = s.rsid.map(pos.pos38).astype(int)
    s = s[~((s.chr == "6") & s.pos38.between(25_000_000, 34_000_000))]
    lead = clump(s, "chr", "pos38"); s = s[s.rsid.isin(lead)].copy()
    s["exposure"] = k; s["F"] = (s.beta / s.se) ** 2; rows.append(s)
    print(k, "candidates", len(cand[k]), "mapped", int(d.rsid.isin(cand[k]).sum()), "instruments", len(s), "median F", round(s.F.median(), 1), flush=True)
I = pd.concat(rows)[["exposure", "rsid", "chr", "pos38", "ea", "oa", "eaf", "beta", "se", "p", "n", "F"]]
I.to_csv(f"{ROOT}/results/mr01_instruments.csv", index=False)
print(I.groupby("exposure").agg(n=("rsid", "size"), F_min=("F", "min"), F_med=("F", "median")))
