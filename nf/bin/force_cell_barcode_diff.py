#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path

import pandas as pd


def _read_lines(path: Path):
    if str(path).endswith(".gz"):
        opener = gzip.open
        mode = "rt"
    else:
        opener = open
        mode = "r"
    out = []
    with opener(path, mode) as f:
        for line in f:
            s = line.strip()
            if s:
                out.append(s)
    return out


def _read_whitelist(cbcsv: Path):
    df = pd.read_csv(cbcsv, header=0)
    if "gex_cb" not in df.columns or "m_cb" not in df.columns:
        raise ValueError(f"Invalid whitelist columns in {cbcsv}, expected gex_cb and m_cb")
    df = df[["gex_cb", "m_cb"]].dropna()
    df["gex_cb"] = df["gex_cb"].astype(str)
    df["m_cb"] = df["m_cb"].astype(str)
    gex_to_m = dict(zip(df["gex_cb"], df["m_cb"]))
    m_to_gex = dict(zip(df["m_cb"], df["gex_cb"]))
    return gex_to_m, m_to_gex


def _write_list(path: Path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for v in values:
            f.write(f"{v}\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gex_barcodes", required=True, help="RNA forced barcodes.tsv.gz (one column, no header)")
    p.add_argument("--methy_counts_csv", required=True, help="Previous filtered_barcode_reads_counts.csv (has barcode column)")
    p.add_argument("--chemistry", default=None, help="Chemistry name, e.g. DD-MET3 or DD-MET5")
    p.add_argument("--cbcsv", default=None, help="Whitelist csv with columns gex_cb,m_cb (needed for DD-MET3)")
    p.add_argument("--outdir", default=".", help="Output directory")
    args = p.parse_args()

    outdir = Path(args.outdir).resolve()
    gex_barcodes_path = Path(args.gex_barcodes).resolve()
    methy_counts_path = Path(args.methy_counts_csv).resolve()
    chemistry = (args.chemistry or "").strip()
    cbcsv = Path(args.cbcsv).resolve() if args.cbcsv else None

    gex_barcodes = _read_lines(gex_barcodes_path)
    gex_barcodes_set = set(gex_barcodes)

    methy_df = pd.read_csv(methy_counts_path, header=0)
    if "barcode" not in methy_df.columns:
        raise ValueError(f"Missing barcode column in {methy_counts_path}")
    methy_barcodes = methy_df["barcode"].dropna().astype(str).tolist()
    methy_barcodes_set = set(methy_barcodes)

    use_mapping = (chemistry == "DD-MET3") or (cbcsv is not None and chemistry != "")
    if use_mapping:
        if cbcsv is None:
            raise ValueError("chemistry indicates mapping is needed but --cbcsv not provided")
        gex_to_m, m_to_gex = _read_whitelist(cbcsv)
        mapped = [gex_to_m.get(cb) for cb in gex_barcodes if cb in gex_to_m]
        gex_missing = [cb for cb in gex_barcodes if cb not in gex_to_m]
        a_methy = [m for m in mapped if m]
        a_methy_set = set(a_methy)
    else:
        m_to_gex = None
        gex_missing = []
        a_methy_set = set(gex_barcodes_set)

    add_methy = sorted(a_methy_set - methy_barcodes_set)
    drop_methy = sorted(methy_barcodes_set - a_methy_set)
    target_methy = sorted(a_methy_set)

    if use_mapping:
        add_gex = sorted({m_to_gex[m] for m in add_methy if m in m_to_gex})
    else:
        add_gex = add_methy

    _write_list(outdir / "drop_barcodes.tsv", drop_methy)
    _write_list(outdir / "add_barcodes.tsv", add_methy)
    _write_list(outdir / "add_gex_barcodes.tsv", add_gex)
    _write_list(outdir / "target_methy_barcodes.tsv", target_methy)

    meta_path = outdir / "barcode_diff.meta.tsv"
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    with open(meta_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["gex_barcodes_total", len(gex_barcodes_set)])
        w.writerow(["methy_barcodes_total", len(methy_barcodes_set)])
        w.writerow(["target_methy_total", len(a_methy_set)])
        w.writerow(["add_methy_total", len(add_methy)])
        w.writerow(["drop_methy_total", len(drop_methy)])
        w.writerow(["add_gex_total", len(add_gex)])
        w.writerow(["gex_missing_in_whitelist", len(gex_missing)])


if __name__ == "__main__":
    main()
