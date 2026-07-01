#!/usr/bin/env python
import argparse
import glob
import pathlib
import shutil
import sys

import yaml
import xarray as xr

from ALLCools.mcds import MCDS
from ALLCools.mcds.utilities import update_dataset_config


def _load_allcools_config(dataset_dir):
    cfg_path = pathlib.Path(dataset_dir) / ".ALLCools"
    if not cfg_path.exists():
        return None
    with open(cfg_path) as f:
        return yaml.safe_load(f)


def _load_barcodes(path):
    lines = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s:
                lines.append(s)
    seen = set()
    ordered = []
    for s in lines:
        if s not in seen:
            seen.add(s)
            ordered.append(s)
    return ordered


def _expand_inputs(inputs):
    raw = []
    for item in inputs:
        if any(ch in item for ch in ["*", "?", "["]):
            raw.extend(glob.glob(item))
        else:
            raw.append(item)
    raw = [str(pathlib.Path(p).absolute()) for p in raw]
    raw = sorted(list(dict.fromkeys(raw)))
    if len(raw) == 0:
        raise ValueError("No input paths matched.")

    dataset_dirs = []
    for p in raw:
        cfg = pathlib.Path(p) / ".ALLCools"
        if cfg.exists():
            dataset_dirs.append(str(pathlib.Path(p).absolute()))
            continue
        found = glob.glob(str(pathlib.Path(p) / "**/.ALLCools"), recursive=True)
        if len(found) == 1:
            dataset_dirs.append(str(pathlib.Path(found[0]).parent.absolute()))
        elif len(found) == 0:
            raise ValueError(f"No .ALLCools found under: {p}")
        else:
            raise ValueError(f"Multiple .ALLCools found under: {p}. Please point inputs to specific dataset dirs.")
    return dataset_dirs


def _format_examples(values, limit=5):
    return ", ".join([str(v) for v in values[:limit]])


def _validate_no_duplicate_obs_within_shards(input_paths, first_cfg, var_dim_keys, obs_dim):
    for region_key in var_dim_keys:
        region_dim = first_cfg["ds_region_dim"][region_key]
        obs_dim_this = first_cfg.get("ds_sample_dim", {}).get(region_key, obs_dim)
        for shard_dir in input_paths:
            src = pathlib.Path(shard_dir) / region_dim
            if not src.exists():
                raise FileNotFoundError(f"Missing dataset subdir {region_dim} in shard: {shard_dir}")
            ds = xr.open_zarr(str(src))
            obs_index = list(ds.get_index(obs_dim_this))
            duplicate_in_shard = []
            seen_in_shard = set()
            for value in obs_index:
                if value in seen_in_shard and value not in duplicate_in_shard:
                    duplicate_in_shard.append(value)
                seen_in_shard.add(value)
            if len(duplicate_in_shard) > 0:
                raise ValueError(
                    f"Duplicate {obs_dim_this} values found within MCDS shard {shard_dir} "
                    f"for region {region_key}: {_format_examples(duplicate_in_shard)}. "
                    "Refusing to merge because duplicate-cell aggregation within a single shard is not implemented."
                )


def _validate_no_obs_overlap_across_shards(input_paths, region_dim, obs_dim):
    seen = set()
    for shard_dir in input_paths:
        src = pathlib.Path(shard_dir) / region_dim
        if not src.exists():
            raise FileNotFoundError(f"Missing dataset subdir {region_dim} in shard: {shard_dir}")
        ds = xr.open_zarr(str(src))
        obs_index = list(ds.get_index(obs_dim))
        duplicate_across_shards = []
        for value in obs_index:
            if value in seen and value not in duplicate_across_shards:
                duplicate_across_shards.append(value)
        if len(duplicate_across_shards) > 0:
            raise ValueError(
                f"Overlapping {obs_dim} values found across MCDS shards for region {region_dim}: "
                f"{_format_examples(duplicate_across_shards)}. "
                "Refusing to merge because MCDS-level aggregation cannot globally deduplicate "
                "UR-tagged molecules. Merge per-cell BAM/ALLC before generating MCDS."
            )
        seen.update(obs_index)


def _write_merged_region(input_paths, target, region_dim, obs_dim):
    for i, shard_dir in enumerate(input_paths):
        src = pathlib.Path(shard_dir) / region_dim
        if not src.exists():
            raise FileNotFoundError(f"Missing dataset subdir {region_dim} in shard: {shard_dir}")
        ds = xr.open_zarr(str(src))
        if i == 0:
            ds.to_zarr(str(target), mode="w")
        else:
            ds.to_zarr(str(target), append_dim=obs_dim)


def merge_mcds_shards(
    inputs,
    output_path,
    obs_dim="cell",
    chunks="auto",
    split_large_chunks=False,
    overwrite=False,
    merge_mode="stream",
):
    input_paths = _expand_inputs(inputs)
    first_cfg = _load_allcools_config(input_paths[0])
    if first_cfg is None:
        raise ValueError(
            f"Input does not look like an ALLCools MCDS directory (missing .ALLCools): {input_paths[0]}"
        )

    for p in input_paths[1:]:
        cfg = _load_allcools_config(p)
        if cfg is None:
            raise ValueError(f"Mixed input types; missing .ALLCools in: {p}")
        if cfg.get("ds_region_dim") != first_cfg.get("ds_region_dim") or cfg.get("ds_sample_dim") != first_cfg.get(
            "ds_sample_dim"
        ):
            raise ValueError(f"Input shard configs are not consistent: {input_paths[0]} vs {p}")

    var_dim_keys = list(first_cfg.get("ds_region_dim", {}).keys())
    if len(var_dim_keys) == 0:
        raise ValueError(f"Invalid .ALLCools config (empty ds_region_dim): {input_paths[0]}")
    var_dims = list(dict.fromkeys([first_cfg["ds_region_dim"][k] for k in var_dim_keys]))

    _validate_no_duplicate_obs_within_shards(input_paths, first_cfg, var_dim_keys, obs_dim)
    for region_key in var_dim_keys:
        region_dim = first_cfg["ds_region_dim"][region_key]
        obs_dim_this = first_cfg.get("ds_sample_dim", {}).get(region_key, obs_dim)
        _validate_no_obs_overlap_across_shards(input_paths, region_dim, obs_dim_this)

    output_path = pathlib.Path(output_path).absolute()
    if output_path.exists():
        if overwrite:
            if output_path.is_dir():
                shutil.rmtree(output_path)
            else:
                output_path.unlink()
        else:
            raise FileExistsError(f"output_path already exists: {output_path}")

    if merge_mode == "openwrite":
        merged = MCDS.open(
            input_paths, obs_dim=obs_dim, var_dim=var_dim_keys, chunks=chunks, split_large_chunks=split_large_chunks
        )
        merged.write_dataset(str(output_path), mode="w-", obs_dim=obs_dim, var_dims=var_dims, chunks=chunks)
    elif merge_mode == "stream":
        output_path.mkdir(parents=True, exist_ok=True)
        for region_key in var_dim_keys:
            region_dim = first_cfg["ds_region_dim"][region_key]
            obs_dim_this = first_cfg.get("ds_sample_dim", {}).get(region_key, obs_dim)
            target = output_path / region_dim
            _write_merged_region(input_paths, target, region_dim, obs_dim_this)
    else:
        raise ValueError(f"Unsupported merge_mode: {merge_mode}")

    chrom_sizes_src = pathlib.Path(input_paths[0]) / "chrom_sizes.txt"
    if chrom_sizes_src.exists():
        shutil.copyfile(chrom_sizes_src, output_path / "chrom_sizes.txt")

    update_dataset_config(str(output_path), config=first_cfg)
    return str(output_path)


def subset_mcds_by_barcodes(
    input_path,
    barcode_path,
    output_path,
    obs_dim="cell",
    chunks="auto",
    strict=False,
    mode="stream",
    batch_size=5000,
):
    input_path = pathlib.Path(input_path).absolute()
    output_path = pathlib.Path(output_path).absolute()
    cfg = _load_allcools_config(input_path)
    if cfg is None:
        raise ValueError(f"Input does not look like an ALLCools MCDS directory: {input_path}")
        
    var_keys = list(cfg.get("ds_region_dim", {}).keys())
    var_dims = list(cfg.get("ds_region_dim", {}).values())

    barcodes = _load_barcodes(barcode_path)
    if output_path.exists():
        raise FileExistsError(str(output_path))

    if mode == "openwrite":
        ds = MCDS.open(str(input_path), obs_dim=obs_dim, var_dim=var_keys, use_obs=barcodes, chunks=chunks)
        idx = set(ds.get_index(obs_dim))
        missing = [b for b in barcodes if b not in idx]
        if strict and len(missing) > 0:
            raise KeyError(f"Missing {len(missing)} barcodes")
        ds.write_dataset(str(output_path), mode="w-", obs_dim=obs_dim, var_dims=var_dims, use_obs=barcodes, chunks=chunks)
    elif mode == "stream":
        output_path.mkdir(parents=True, exist_ok=True)
        for region_key, region_dim in zip(var_keys, var_dims):
            src = input_path / region_dim
            if not src.exists():
                raise FileNotFoundError(f"Missing dataset subdir {region_dim} in input: {input_path}")
            ds = xr.open_zarr(str(src))
            obs_index = ds.get_index(obs_dim)
            exist = [b for b in barcodes if b in set(obs_index)]
            if strict:
                missing = [b for b in barcodes if b not in set(obs_index)]
                if len(missing) > 0:
                    raise KeyError(f"{region_dim}: Missing {len(missing)} barcodes")
            target = output_path / region_dim
            for i in range(0, len(exist), int(batch_size)):
                batch = exist[i : i + int(batch_size)]
                sub = ds.sel({obs_dim: batch}).load()
                if i == 0:
                    sub.to_zarr(str(target), mode="w")
                else:
                    sub.to_zarr(str(target), append_dim=obs_dim)
    else:
        raise ValueError(f"Unsupported mode: {mode}. Use one of ['stream', 'openwrite'].")

    chrom_sizes = input_path / "chrom_sizes.txt"
    if chrom_sizes.exists():
        shutil.copyfile(chrom_sizes, output_path / "chrom_sizes.txt")
    update_dataset_config(
        str(output_path),
        config={"region_dim": None, "ds_region_dim": cfg.get("ds_region_dim", {}), "ds_sample_dim": cfg.get("ds_sample_dim", {})},
    )
    return str(output_path)


def main():
    parser = argparse.ArgumentParser(description="MCDS Utility: Merge or Subset MCDS datasets")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Command to run")

    # Merge subcommand
    parser_merge = subparsers.add_parser("merge", help="Merge multiple MCDS shards")
    parser_merge.add_argument(
        "inputs",
        nargs="+",
        help="Input shard MCDS directories, supports glob (each contains a .ALLCools file).",
    )
    parser_merge.add_argument("--output_path", required=True, help="Output merged MCDS directory path.")
    parser_merge.add_argument("--obs_dim", default="cell")
    parser_merge.add_argument("--chunks", default="auto")
    parser_merge.add_argument("--split_large_chunks", action="store_true")
    parser_merge.add_argument("--overwrite", action="store_true")
    parser_merge.add_argument(
        "--merge_mode",
        choices=["stream", "openwrite"],
        default="stream",
        help="stream: append shards sequentially with low memory; openwrite: open all shards then write once.",
    )

    # Subset subcommand
    parser_subset = subparsers.add_parser("subset", help="Subset MCDS dataset by barcodes")
    parser_subset.add_argument("--input_path", required=True, help="Input MCDS directory path")
    parser_subset.add_argument("--barcode_path", required=True, help="Path to barcode list file")
    parser_subset.add_argument("--output_path", required=True, help="Output MCDS directory path")
    parser_subset.add_argument("--obs_dim", default="cell")
    parser_subset.add_argument("--chunks", default="auto")
    parser_subset.add_argument("--strict", action="store_true")
    parser_subset.add_argument("--mode", choices=["stream", "openwrite"], default="stream")
    parser_subset.add_argument("--batch_size", type=int, default=5000, help="Batch size for stream mode")

    args = parser.parse_args()
    
    chunks = args.chunks
    if chunks is not None and str(chunks).lower() in ["none", "null"]:
        chunks = None

    if args.command == "merge":
        merge_mcds_shards(
            inputs=args.inputs,
            output_path=args.output_path,
            obs_dim=args.obs_dim,
            chunks=chunks,
            split_large_chunks=args.split_large_chunks,
            overwrite=args.overwrite,
            merge_mode=args.merge_mode,
        )
    elif args.command == "subset":
        subset_mcds_by_barcodes(
            input_path=args.input_path,
            barcode_path=args.barcode_path,
            output_path=args.output_path,
            obs_dim=args.obs_dim,
            chunks=chunks,
            strict=args.strict,
            mode=args.mode,
            batch_size=args.batch_size,
        )

if __name__ == "__main__":
    main()
