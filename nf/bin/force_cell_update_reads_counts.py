#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def _read_barcode_list(path: Path):
    out = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s:
                out.append(s)
    seen = set()
    ordered = []
    for s in out:
        if s not in seen:
            seen.add(s)
            ordered.append(s)
    return ordered


def _load_reads_counts_csv(path: Path):
    df = pd.read_csv(path, header=0)
    cols = list(df.columns)
    if "barcode" in cols and "reads_counts" in cols:
        df = df[["reads_counts", "barcode"]]
    elif len(cols) >= 2:
        df = df.iloc[:, :2]
        df.columns = ["reads_counts", "barcode"]
    else:
        raise ValueError(f"Invalid reads counts file: {path}")
    df["barcode"] = df["barcode"].astype(str)
    df["reads_counts"] = pd.to_numeric(df["reads_counts"], errors="coerce").fillna(0).astype(int)
    df = df.dropna(subset=["barcode"])
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--old_counts_csv", required=True, help="Previous filtered_barcode_reads_counts.csv")
    p.add_argument("--drop_barcodes", required=True, help="drop_barcodes.tsv (one barcode per line)")
    p.add_argument("--add_barcodes", required=True, help="add_barcodes.tsv (one barcode per line)")
    p.add_argument(
        "--add_reads_counts_csv",
        action="append",
        default=[],
        help="Reads counts csv(s) from newly generated barcodes; can be repeated",
    )
    p.add_argument("--outdir", default=".", help="Output directory")
    args = p.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    old_df = _load_reads_counts_csv(Path(args.old_counts_csv).resolve())
    drop = set(_read_barcode_list(Path(args.drop_barcodes).resolve()))
    add = set(_read_barcode_list(Path(args.add_barcodes).resolve()))

    old_df = old_df[~old_df["barcode"].isin(drop)].copy()
    old_df = old_df[~old_df["barcode"].isin(add)].copy()

    add_parts = []
    for pth in args.add_reads_counts_csv:
        add_parts.append(_load_reads_counts_csv(Path(pth).resolve()))
    if add_parts:
        add_df = pd.concat(add_parts, axis=0, ignore_index=True)
        add_df = add_df[add_df["barcode"].isin(add)]
        add_df = add_df.groupby("barcode", as_index=False)["reads_counts"].sum()
    else:
        add_df = pd.DataFrame({"barcode": [], "reads_counts": []})

    merged = pd.concat([old_df, add_df[["reads_counts", "barcode"]]], axis=0, ignore_index=True)
    merged["reads_counts"] = pd.to_numeric(merged["reads_counts"], errors="coerce").fillna(0).astype(int)
    merged["barcode"] = merged["barcode"].astype(str)
    merged = merged.drop_duplicates(subset=["barcode"], keep="last")
    merged = merged.sort_values(by="reads_counts", ascending=False)

    merged_path = outdir / "filtered_barcode_reads_counts.csv"
    merged[["reads_counts", "barcode"]].to_csv(merged_path, index=False)

    filtered_barcode_path = outdir / "filtered_barcode"
    with open(filtered_barcode_path, "w") as f:
        for bc in merged["barcode"].tolist():
            f.write(f"{bc}\n")


if __name__ == "__main__":
    main()
