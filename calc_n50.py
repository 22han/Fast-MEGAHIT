#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Assembly Pipeline: Contig Refinement + Assembly Statistics
Usage:
  # Full pipeline (refine contigs → calc stats on output):
  python asm_pipeline.py -i contigs.fa -r reads_R1.fq -r reads_R2.fq -o filtered.fa [--min-len 200]

  # Stats-only mode (skip refinement, just calculate N50 etc.):
  python asm_pipeline.py -s contigs.fa
"""

import sys
import os
import subprocess
import pandas as pd
import argparse
import gzip
from Bio import SeqIO


# ═══════════════════════════════════════════════════════════════
#  Utility: safely execute shell command
# ═══════════════════════════════════════════════════════════════
def run_cmd(cmd, desc=""):
    """安全执行 shell 命令，带错误捕获"""
    try:
        subprocess.run(cmd, shell=True, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] {desc} failed:\n{e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)


# ═══════════════════════════════════════════════════════════════
#  Part 1: Assembly Statistics (N50, N90, L50, L90 …)
# ═══════════════════════════════════════════════════════════════
def calc_assembly_stats(fasta_file):
    """从 FASTA 文件计算组装统计指标，并打印 & 返回字典"""
    lengths = []
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    with open_func(fasta_file, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            lengths.append(len(record.seq))

    if not lengths:
        print("⚠️  No contigs found in the input file!")
        return

    lengths.sort(reverse=True)
    total_length = sum(lengths)
    cumsum = 0

    n50 = n90 = l50 = l90 = None
    for i, length in enumerate(lengths, 1):
        cumsum += length
        if n50 is None and cumsum >= total_length * 0.5:
            n50 = length
            l50 = i
        if n90 is None and cumsum >= total_length * 0.9:
            n90 = length
            l90 = i
            break

    print(f"\n📊 Assembly Statistics for {fasta_file}")
    print(f"{'='*55}")
    print(f"  Total contigs      : {len(lengths):,}")
    print(f"  Total length (bp)  : {total_length:,}")
    print(f"  Longest contig     : {lengths[0]:,} bp")
    print(f"  Shortest contig    : {lengths[-1]:,} bp")
    print(f"  Mean length        : {total_length/len(lengths):.1f} bp")
    print(f"  Median length      : {lengths[len(lengths)//2]:,} bp")
    print(f"{'-'*55}")
    print(f"  N50                : {n50:,} bp  (L50 = {l50:,} contigs)")
    print(f"  N90                : {n90:,} bp  (L90 = {l90:,} contigs)")
    print(f"{'='*55}\n")

    return {
        'n_contigs': len(lengths),
        'total_length': total_length,
        'longest': lengths[0],
        'shortest': lengths[-1],
        'n50': n50, 'l50': l50,
        'n90': n90, 'l90': l90,
    }


# ═══════════════════════════════════════════════════════════════
#  Part 2: Mapping + Coverage calculation
# ═══════════════════════════════════════════════════════════════
def map_and_calc_cov(contigs_fa, reads_list, out_bam, threads=8):
    """比对 reads 并计算平均覆盖度（使用 samtools coverage）"""
    reads_str = " ".join(reads_list)

    cmd_map = (
        f"minimap2 -ax sr -t {threads} {contigs_fa} {reads_str} "
        f"| samtools view -bS - "
        f"| samtools sort -@ {threads} -o {out_bam}"
    )
    run_cmd(cmd_map, "Minimap2 mapping + sorting")
    run_cmd(f"samtools index -@ {threads} {out_bam}", "BAM indexing")

    cov_tsv = out_bam.replace(".bam", ".cov.tsv")
    run_cmd(f"samtools coverage -o {cov_tsv} -H {out_bam}", "Coverage calculation")
    return cov_tsv


# ═══════════════════════════════════════════════════════════════
#  Part 3: Adaptive contig refinement
# ═══════════════════════════════════════════════════════════════
def refine_contigs(cov_tsv, contigs_fa, out_fa,
                   min_len=200, rescue_factor=0.5):
    """基于长度的分层覆盖度过滤"""
    if not os.path.exists(cov_tsv) or os.path.getsize(cov_tsv) == 0:
        print("[ERROR] Coverage file is empty or not found.", file=sys.stderr)
        sys.exit(1)

    print("📏 Extracting contig lengths from FASTA...")
    contig_lengths = {}
    open_func = gzip.open if contigs_fa.endswith('.gz') else open
    with open_func(contigs_fa, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            contig_lengths[record.id.split()[0]] = len(record.seq)

    df = pd.read_csv(
        cov_tsv, sep="\t", comment="#", header=None,
        names=["rname", "startpos", "endpos", "numreads", "covbases",
               "coverage_pct", "meandepth", "meanbaseq", "meanmapq"]
    )
    df = df.rename(columns={"meandepth": "coverage", "rname": "contig"})
    df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce").fillna(0.0)

    df["length"] = df["contig"].map(contig_lengths)
    df = df.dropna(subset=["length"])
    df["length"] = df["length"].astype(int)
    total_before = len(df)

    print(f"📊 Coverage stats:  min={df['coverage'].min():.2f}x, "
          f"Q1={df['coverage'].quantile(0.25):.2f}x, "
          f"median={df['coverage'].median():.2f}x, "
          f"max={df['coverage'].max():.2f}x")

    df = df[df["length"] >= min_len]
    valid_cov = df[df["coverage"] > 0]["coverage"]
    dynamic_cutoff = (max(0.3, valid_cov.quantile(0.25) * 0.5)
                      if len(valid_cov) > 0 else 0.3)

    mask_main = df["coverage"] >= dynamic_cutoff
    mask_rescue = ((df["length"] >= 500) &
                   (df["coverage"] >= dynamic_cutoff * rescue_factor))
    final_df = df[mask_main | mask_rescue]
    keep_set = set(final_df["contig"])

    removed_len = total_before - len(df)
    removed_cov = len(df) - len(final_df)
    rescued    = int((mask_rescue & ~mask_main).sum())

    print(f"\n{'='*45}")
    print(f"📊 Contig Refinement Summary")
    print(f"{'='*45}")
    print(f"  Total contigs (≥ {min_len} bp)       : {total_before}")
    print(f"  Removed (length < {min_len} bp)      : {removed_len}")
    print(f"  Dynamic coverage threshold           : {dynamic_cutoff:.2f}x")
    print(f"  Removed (low cov & short)            : {removed_cov}")
    print(f"  Rescued (medium len + low cov)       : {rescued}")
    print(f"  ✅ Final retained contigs            : {len(final_df)}")
    print(f"{'='*45}\n")

    if len(final_df) == 0:
        print("⚠️  WARNING: All contigs were filtered out. "
              "Check mapping rate or lower the threshold.\n")

    with open_func(contigs_fa, "rt") as f_in, open(out_fa, "w") as f_out:
        for record in SeqIO.parse(f_in, "fasta"):
            seq_id = record.id.split()[0]
            if seq_id in keep_set:
                SeqIO.write(record, f_out, "fasta")


# ═══════════════════════════════════════════════════════════════
#  Utility: clean up temp files
# ═══════════════════════════════════════════════════════════════
def cleanup(files):
    for f in files:
        try:
            if os.path.exists(f):
                os.remove(f)
        except Exception:
            pass


# ═══════════════════════════════════════════════════════════════
#  Main entry
# ═══════════════════════════════════════════════════════════════
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="One-stop assembly pipeline: "
                    "coverage-aware contig refinement + assembly statistics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full pipeline  (refine → stats on output)
  python asm_pipeline.py -i contigs.fa \
      -r reads_R1.fq -r reads_R2.fq \
      -o filtered.fa

  # Stats only (skip refinement)
  python asm_pipeline.py -s contigs.fa
        """
    )
    parser.add_argument("-i", "--input", help="Input contigs FASTA (.fa / .fa.gz)")
    parser.add_argument("-r", "--reads", nargs="+", help="Read files (R1, or R1 R2)")
    parser.add_argument("-o", "--output", help="Output filtered FASTA")
    parser.add_argument("--min-len", type=int, default=200,
                        help="Minimum contig length (default: 200)")
    parser.add_argument("--threads", type=int, default=8,
                        help="CPU threads for minimap2/samtools (default: 8)")
    parser.add_argument("-s", "--stats", help="FASTA file to compute stats only "
                        "(mutually exclusive with -i/-r/-o)")

    args = parser.parse_args()

    # ── Mode 1: stats only ──────────────────────────────────
    if args.stats:
        if args.input or args.reads or args.output:
            parser.error("-s/--stats is mutually exclusive with -i/-r/-o")
        calc_assembly_stats(args.stats)
        sys.exit(0)

    # ── Mode 2: full pipeline ───────────────────────────────
    if not args.input or not args.reads or not args.output:
        parser.error("Full pipeline requires -i, -r, and -o "
                     "(or use -s for stats-only mode)")

    temp_bam = "temp_align.bam"
    temp_cov = "temp_align.cov.tsv"
    temp_files = [temp_bam, f"{temp_bam}.bai", temp_cov]

    try:
        print(f"🚀 [1/3] Mapping reads to {args.input} ...")
        cov_file = map_and_calc_cov(args.input, args.reads,
                                    temp_bam, args.threads)

        print(f"🔍 [2/3] Applying adaptive coverage filter ...")
        refine_contigs(cov_file, args.input, args.output,
                       args.min_len)

        print(f"📊 [3/3] Calculating assembly statistics on output ...")
        calc_assembly_stats(args.output)

        print(f"✨ Pipeline complete! Refined assembly saved to {args.output}")
    finally:
        cleanup(temp_files)
