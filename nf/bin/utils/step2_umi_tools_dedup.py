import os
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_cmd(cmd):
    subprocess.run(cmd, check=True)

def dedup_bam(bam_path, out_dir, samtools_threads, index_threads, skip_existing):
    bam_path = bam_path.strip()
    if not bam_path or bam_path.startswith('#'):
        return bam_path, False, 'skipped-empty'
    if not os.path.exists(bam_path):
        return bam_path, False, 'not-found'
    base = os.path.splitext(os.path.basename(bam_path))[0]
    out_dir_use = out_dir or os.path.dirname(bam_path) or '.'
    os.makedirs(out_dir_use, exist_ok=True)

    sorted_bam = os.path.join(out_dir_use, f"{base}_sorted.bam")
    sorted_bai = sorted_bam + ".bai"
    stats = os.path.join(out_dir_use, f"{base}.stats")
    log = os.path.join(out_dir_use, f"{base}.log")
    err_log = os.path.join(out_dir_use, f"{base}.err.log")
    dedup_bam_path = os.path.join(out_dir_use, f"{base}.dedup.bam")

    try:
        if skip_existing and os.path.exists(dedup_bam_path):
            try:
                run_cmd(["samtools", "index", "-@", str(index_threads), dedup_bam_path])
            except subprocess.CalledProcessError:
                pass
            return bam_path, True, 'done-exist'

        run_cmd(["samtools", "sort", "-@", str(samtools_threads), "-o", sorted_bam, bam_path])
        run_cmd(["samtools", "index", "-@", str(index_threads), sorted_bam])
        run_cmd([
            "umi_tools", "dedup",
            "--output-stats", stats,
            "--extract-umi-method", "tag",
            "--umi-tag", "UR",
            "--paired",
            "--ignore-tlen",
            "-I", sorted_bam,
            "-L", log,
            "-E", err_log,
            "-S", dedup_bam_path,
        ])

        for p in [sorted_bam, sorted_bai]:
            if os.path.exists(p):
                try:
                    os.remove(p)
                except Exception:
                    pass
        for name in os.listdir(out_dir_use):
            if name.startswith(base) and name.endswith('.tsv'):
                fp = os.path.join(out_dir_use, name)
                try:
                    os.remove(fp)
                except Exception:
                    pass
        for name in os.listdir(out_dir_use):
            if name.startswith(base) and name.endswith('.log'):
                fp = os.path.join(out_dir_use, name)
                try:
                    os.remove(fp)
                except Exception:
                    pass
        return bam_path, True, 'done'
    except subprocess.CalledProcessError as e:
        return bam_path, False, f'failed:{e}'

def main():
    parser = argparse.ArgumentParser(description='Deduplicate BAM files using umi_tools')
    parser.add_argument('--bam_list', type=str, required=True)
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--max_workers', type=int, default=max(4, (os.cpu_count() or 4)))
    parser.add_argument('--samtools_threads', type=int, default=1)
    parser.add_argument('--index_threads', type=int, default=1)
    parser.add_argument('--skip_existing', action='store_true')
    args = parser.parse_args()

    with open(args.bam_list, 'r') as fh:
        bam_paths = [line.strip() for line in fh if line.strip()]

    results = []
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        futs = [
            ex.submit(
                dedup_bam,
                bp,
                args.out_dir,
                args.samtools_threads,
                args.index_threads,
                args.skip_existing,
            )
            for bp in bam_paths
        ]
        for fut in as_completed(futs):
            results.append(fut.result())

    failed = [r for r in results if not r[1]]
    print(f"Total: {len(results)}, Succeeded: {len(results)-len(failed)}, Failed: {len(failed)}")
    for bp, ok, msg in failed:
        print(f"FAIL: {bp} -> {msg}")

if __name__ == '__main__':
    main()