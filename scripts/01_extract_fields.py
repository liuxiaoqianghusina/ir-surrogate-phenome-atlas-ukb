"""t01：IR 指数图谱预试验 —— 提取结局（First occurrence）与身高、舒张压、高血压诊断。"""
import pandas as pd, csv, re, time
M = "/path/to/ukb_data"
t0 = time.time()
POS = ["E11", "I21", "I25", "I63", "I50", "I48", "N18", "K76", "M10", "K80", "I70", "I73", "I10", "I35", "I26", "I80", "G47", "F03", "G30", "K85"]
NEG = ["K40", "K35", "I83", "H40", "K64", "G43", "J30", "L20", "H90", "K29", "G70", "M54", "J45", "L40", "K57", "M17", "H25", "J18", "N39", "F32",
       "H52", "H61", "K02", "L60", "M20", "N60", "J34", "H00", "L70", "K60"]
EXTRA = ["E10", "E12", "E13", "E14"]
fo = {}
for r in csv.DictReader(open(f"{M}/ukb_data_dictionary.csv", encoding="utf-8")):
    m = re.match(r"Date ([A-Z]\d\d) first reported", r["title"])
    if m and m.group(1) in POS + NEG + EXTRA and "First occurrences" in r["folder_path"]: fo[r["name"]] = "fo_" + m.group(1)
H = pd.concat(pd.read_csv(f"{M}/ukb_health_outcomes.csv", usecols=["eid", "p40000_i0", "p40001_i0"] + list(fo), dtype=str, chunksize=100000), ignore_index=True).rename(columns=fo)
H.to_pickle("../cache/fo.pkl"); print("FO", H.shape, round(time.time() - t0), flush=True)
h = pd.read_csv(f"{M}/ukb_assessment_center.csv", nrows=0).columns
use = ["eid"] + [c for c in h if re.fullmatch(r"p(50_i0|4079_i0_a[01]|94_i0_a[01]|2443_i0|6150_i0|20002_i0_a\d+)", c)]
A = pd.concat(pd.read_csv(f"{M}/ukb_assessment_center.csv", usecols=use, dtype=str, chunksize=25000), ignore_index=True)
A.to_pickle("../cache/ac_extra.pkl"); print("AC", A.shape, round(time.time() - t0))
