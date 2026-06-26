#!/usr/bin/env bash
# Deduplicate a position-sorted, gene-annotated BAM by splitting per reference contig, running
# umi_tools dedup on each contig in parallel, and merging the results.
#
# umi_tools dedup --per-cell groups reads by cell barcode + start position, which is local to a
# single contig (reads on different contigs never interact). Splitting by contig therefore yields
# the identical set of deduplicated molecules as whole-BAM dedup: verified on real data to give an
# identical read-count-out and an identical (barcode, gene) multiset (i.e. identical count matrix).
#
# NB: umi_tools' choice of WHICH read represents a UMI group is stateful, so the merged dedup.bam
# keeps an equivalent-but-different representative read (same molecule: same barcode/UMI/position,
# hence same gene) for some groups. All counts and metrics are identical; only the representative
# reads in the published BAM differ.
#
# Usage: dedup_by_contig.sh <bam> <sample_id> <cpus>
set -euo pipefail

bam=$1
sample_id=$2
cpus=$3

# The dedup input is expected to be entirely mapped, gene-annotated reads. Per-contig splitting
# would silently drop unmapped reads, so fail loud if any are present (trust upstream; catch drift).
unmapped=$(samtools view -c -f 4 "$bam")
if [[ "$unmapped" -ne 0 ]]; then
  echo "ERROR dedup_by_contig: $bam has $unmapped unmapped reads; per-contig split would drop them" >&2
  exit 1
fi

# Contigs that carry reads (idxstats column 3 = mapped read count).
mapfile -t contigs < <(samtools idxstats "$bam" | awk '$3 > 0 {print $1}')
if [[ ${#contigs[@]} -eq 0 ]]; then
  echo "ERROR dedup_by_contig: no mapped contigs in $bam" >&2
  exit 1
fi

dedup_one() {
  local c="$1"
  samtools view -b "$bam" "$c" > "part.${c}.bam"
  samtools index "part.${c}.bam"
  umi_tools dedup --per-cell --in-sam -I "part.${c}.bam" --log="part.${c}.log" > "part.${c}.dedup.bam"
}
export -f dedup_one
export bam

# Dedup each contig in parallel, bounded by the allocated cpus.
printf '%s\n' "${contigs[@]}" | xargs -P "$cpus" -I {} bash -c 'set -e; dedup_one "$1"' _ {}

# Merge per-contig dedup bams in idxstats (header) order. Identical headers => samtools cat is exact
# and the result stays coordinate-sorted (each part is sorted, contigs in header order).
parts=()
for c in "${contigs[@]}"; do parts+=("part.${c}.dedup.bam"); done
samtools cat -o "${sample_id}.dedup.bam" "${parts[@]}"

# Rebuild the two log fields summary_statistics.py parses (it reads the last token of each line):
# reads_before/after_deduplication are the per-contig sums.
reads_in=$(cat part.*.log | awk '/Input Reads:/ {s += $NF} END {print s + 0}')
reads_out=$(cat part.*.log | awk '/Number of reads out:/ {s += $NF} END {print s + 0}')
printf 'INFO Reads: Input Reads: %s\nINFO Number of reads out: %s\n' "$reads_in" "$reads_out" > "${sample_id}.dedup.log"
