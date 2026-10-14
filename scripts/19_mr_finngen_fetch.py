"""mr02：按 GRCh38 位置从 FinnGen R12 公开汇总数据（GCS，bgzip + tabix）做 HTTP 范围读取，只取工具变量所在的行，不下载整个文件。
解析本地 .tbi 的线性索引得到虚拟偏移，再用 Range 请求取相应 BGZF 块并解压（方法同肠道激素轴项目 m02）。多线程。用 python3 -P 运行。
输出 results/mr02_finngen_lookup.csv。"""
import gzip, struct, zlib, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
ROOT = "/path/to/project"
FG = f"{ROOT}/external/finngen"
URL = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_{}.gz"
ENDPOINTS = ["I9_CHD", "I9_HEARTFAIL", "I9_AF", "I9_STR", "T2D", "M13_GOUT", "K11_CHOLELITH", "NAFLD", "N14_CHRONKIDNEYDIS", "K11_ACUTPANC"]
HEAD = ["chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"]


def read_tbi(path):
    b = gzip.open(path).read()
    assert b[:4] == b"TBI\x01"
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", b[4:36])
    names = b[36:36 + l_nm].split(b"\x00")[:-1]
    p = 36 + l_nm; idx = {}
    for r in range(n_ref):
        n_bin = struct.unpack("<i", b[p:p + 4])[0]; p += 4
        for _ in range(n_bin):
            bin_, n_chunk = struct.unpack("<Ii", b[p:p + 8]); p += 8 + 16 * n_chunk
        n_intv = struct.unpack("<i", b[p:p + 4])[0]; p += 4
        ioff = struct.unpack(f"<{n_intv}Q", b[p:p + 8 * n_intv]); p += 8 * n_intv
        idx[names[r].decode()] = ioff
    return idx


def http_range(url, start, length):
    for a in range(4):
        r = subprocess.run(["curl", "-s", "--max-time", "120", "-r", f"{start}-{start + length - 1}", url], capture_output=True)
        if r.returncode == 0 and r.stdout: return r.stdout
        time.sleep(1 + 2 * a)
    return b""


def parse(text, chrom, pos, window):
    rows = []
    for ln in text.split(b"\n"):
        f = ln.decode(errors="ignore").split("\t")
        if len(f) < 10 or f[0] != chrom: continue
        try:
            p = int(f[1])
        except ValueError:
            continue
        if abs(p - pos) <= window: rows.append(f)
    return rows


def fetch(url, ioff, chrom, pos, window=2):
    """返回 [pos-window, pos+window] 区间内的行。"""
    beg = max(pos - window, 1); k = min((beg - 1) >> 14, len(ioff) - 1)
    voff = ioff[k]; coff, uoff = voff >> 16, voff & 0xFFFF
    buf = http_range(url, coff, 120_000); text = b""; first = True
    while True:
        while buf:
            d = zlib.decompressobj(31)
            try:
                out = d.decompress(buf)
            except zlib.error:
                break
            if not d.eof: break                     # 块不完整：需要更多字节
            if first: out = out[uoff:]; first = False
            text += out; coff += len(buf) - len(d.unused_data); buf = d.unused_data
            last = text.rsplit(b"\n", 2)
            if len(last) >= 2:
                try:
                    lp = int(last[-2].split(b"\t")[1]) if last[-1] == b"" else int(last[-1].split(b"\t")[1])
                    if lp > pos + window: return parse(text, chrom, pos, window)
                except (ValueError, IndexError):
                    pass
        more = http_range(url, coff, len(buf) + 120_000)   # 末尾块不完整：从该块起点重新取
        if not more or len(more) <= len(buf): return parse(text, chrom, pos, window)
        buf = more


if __name__ == "__main__":
    I = pd.read_csv(f"{ROOT}/results/mr01_instruments.csv", dtype={"chr": str})
    Q = I.drop_duplicates("rsid")[["rsid", "chr", "pos38"]].reset_index(drop=True)
    eps = sys.argv[1:] or ENDPOINTS; t0 = time.time(); out = []
    for ep in eps:
        ioff = read_tbi(f"{FG}/finngen_R12_{ep}.gz.tbi"); key = next(iter(ioff))
        def one(q):
            ch = q.chr if not key.startswith("chr") else f"chr{q.chr}"
            ch = "23" if ch == "X" else ch
            return [[ep, q.rsid, q.pos38] + f[:len(HEAD)] for f in fetch(URL.format(ep), ioff[ch], ch, int(q.pos38))]
        with ThreadPoolExecutor(8) as ex:
            for res in ex.map(one, [q for q in Q.itertuples(index=False)]): out.extend(res)
        print(ep, "done", len(out), "rows", round(time.time() - t0), "s", flush=True)
    O = pd.DataFrame(out, columns=["endpoint", "query_rsid", "query_pos38"] + HEAD)
    O.to_csv(f"{ROOT}/results/mr02_finngen_lookup.csv", index=False); print(O.shape)
