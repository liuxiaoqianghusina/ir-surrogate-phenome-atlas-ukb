"""u01：全表型结局 —— 全部 First occurrence（3 位 ICD-10）+ 癌症登记（3 位 C/D0 代码），转换为相对基线的天数。
仅保留分析队列（t02 的 374,814 名基线无糖尿病者）。特殊日期（1900/1901/1902/1903/1909/2037）记为 -99999（日期未知，分析时排除）。
输出：cache/phe_days_partK.csv（每份约 120 个结局）、cache/phe_meta.csv（代码、名称、章节、事件数）、cache/phe_censor.csv。"""
import pandas as pd, numpy as np, csv, re, time
M = "/path/to/ukb_data"
IDEA = "/path/to/upstream"
END = pd.Timestamp("2022-12-31"); SPECIAL = {1900, 1901, 1902, 1903, 1909, 2037}
t0 = time.time()
fo, title = {}, {}
for r in csv.DictReader(open(f"{M}/ukb_data_dictionary.csv", encoding="utf-8")):
    m = re.match(r"Date ([A-Z]\d\d) first reported \((.*)\)", r["title"])
    if m and "First occurrences" in r["folder_path"]:
        fo[r["name"]] = m.group(1); title[m.group(1)] = m.group(2)
print("FO fields:", len(fo), flush=True)
# 癌症登记代码的名称取自 UK Biobank coding 19（ICD-10），例如 C18 → malignant neoplasm of colon
ICD = pd.read_csv("../external/coding19.tsv", sep="\t", dtype=str)
ca_title = {r.coding: r.meaning.split(" ", 1)[1].lower() for r in ICD.itertuples() if re.fullmatch(r"(C\d\d|D0\d)", r.coding)}
coh = pd.read_csv("../cache/ir_atlas.csv", usecols=["eid"]).eid
base = pd.read_pickle(f"{IDEA}/ddw_ppdm/cache/cohort_full.pkl")[["eid", "d0"]].set_index("eid").loc[coh]
h = pd.read_csv(f"{M}/ukb_health_outcomes.csv", nrows=0).columns
can = [c for c in h if re.fullmatch(r"p4000[56]_i\d+", c)]
use = ["eid", "p40000_i0", "p40001_i0"] + list(fo) + can
parts = []
for k, ch in enumerate(pd.read_csv(f"{M}/ukb_health_outcomes.csv", usecols=use, dtype=str, chunksize=40000)):
    ch["eid"] = ch.eid.astype(int); ch = ch[ch.eid.isin(base.index)]
    d0 = base.loc[ch.eid, "d0"].values
    out = pd.DataFrame({"eid": ch.eid.values})
    for col, code in fo.items():
        d = pd.to_datetime(ch[col], errors="coerce")
        days = ((d.values - d0) / np.timedelta64(1, "D"))
        days = np.where(d.dt.year.isin(SPECIAL).values, -99999, days)
        out["fo_" + code] = pd.array(np.where(np.isnan(days), np.nan, np.round(days)), dtype="Int32")
    # 癌症登记：每个 3 位代码的最早日期
    cc = {}
    for i in range(22):
        c, dt = f"p40006_i{i}", f"p40005_i{i}"
        if c not in ch: continue
        code = ch[c].fillna("").str[:3]; dd = pd.to_datetime(ch[dt], errors="coerce")
        days = ((dd.values - d0) / np.timedelta64(1, "D"))
        for cd in code[code.str.match(r"^(C\d\d|D0\d)$")].unique():
            m_ = (code == cd).values
            cur = cc.get(cd, np.full(len(ch), np.nan)); new = np.where(m_, days, np.nan)
            cc[cd] = np.fmin(cur, new)
    for cd, v in cc.items(): out["ca_" + cd] = pd.array(np.where(np.isnan(v), np.nan, np.round(v)), dtype="Int32")
    death = pd.to_datetime(ch.p40000_i0, errors="coerce")
    stop = np.minimum(np.where(death.isna(), END.to_datetime64(), death.values), END.to_datetime64())
    out["cday"] = ((stop - d0) / np.timedelta64(1, "D")).round().astype(int)
    out["death_day"] = pd.array(np.where(death.isna() | (death > END), np.nan, ((death.values - d0) / np.timedelta64(1, "D")).round()), dtype="Int32")
    out["cvdeath"] = ch.p40001_i0.fillna("").str.match(r"I(0\d|1\d|2\d|3\d|4\d|5[0-2]|6\d|7\d)").astype(int).values
    parts.append(out); print(f"chunk {k}: {sum(len(p) for p in parts):,} rows, {round(time.time() - t0)} s", flush=True)
D = pd.concat(parts, ignore_index=True).fillna({c: np.nan for c in []})
cols = [c for c in D.columns if c.startswith(("fo_", "ca_"))]
# 元数据与筛选：基线无该病者中新发事件 ≥200
meta = []
for c in cols:
    v = D[c]; prev = (v <= 0) | (v == -99999)
    inc = (v > 0) & (v <= D.cday)
    code = c[3:]
    meta.append([c, code, (ca_title.get(code, "") if c.startswith("ca_") else title.get(code, "")), code[0], int(inc.sum()), int((prev & (v != -99999)).sum())])
meta = pd.DataFrame(meta, columns=["col", "code", "title", "chapter", "n_incident", "n_prevalent"])
excl = meta.chapter.isin(["O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"])
meta["include"] = (meta.n_incident >= 200) & ~excl
meta.to_csv("../cache/phe_meta.csv", index=False)
keep = meta[meta.include].col.tolist()
print("outcomes eligible:", len(keep), "| by chapter:", meta[meta.include].chapter.value_counts().sort_index().to_dict(), flush=True)
D[["eid", "cday", "death_day", "cvdeath"]].to_csv("../cache/phe_censor.csv", index=False)
for k in range(0, len(keep), 120):
    D[["eid"] + keep[k:k + 120]].to_csv(f"../cache/phe_days_part{k // 120}.csv", index=False)
print("done", round(time.time() - t0), "s")
