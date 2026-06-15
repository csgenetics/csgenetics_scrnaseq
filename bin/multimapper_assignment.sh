#!/usr/bin/env bash
# Name-hash-parallel multimapper gene assignment (WALL-CLOCK optimisation of the #1 critical-path step).
#
# The multimapper assignment is dominated by two `samtools sort -n` passes over the (large, ~1000x
# expanded) multimapper BAM. The gawk groups alignments by read name and picks a canonical, order-
# independent representative per gene, so splitting reads by a hash of the read NAME (all alignments
# of a read land in the same chunk), sorting+assigning each chunk in parallel, and merging yields the
# identical assigned/ambiguous read set as a single whole-BAM pass (verified). This parallelises the
# sort, bounded by the largest chunk.
#
# Usage: multimapper_assignment.sh <feature_count_bam> <sample_id> <gtf> <gawk_script> <cpus>
set -euo pipefail

feature_count_bam=$1
sample_id=$2
gtf=$(readlink -f "$3")
gawk_script=$(readlink -f "$4")
cpus=$5
N=$cpus
[ "$N" -lt 1 ] && N=1
export gtf gawk_script

# 1. filter_for_multimappers_mismatch: NH>1, <=3 mismatches
samtools view -@ "$cpus" -h -b -e '[NH]>1 && ([nM]==0 || [nM]==1 || [nM]==2 || [nM]==3)' -b "$feature_count_bam" > mm.featureCounts.bam
samtools view -H mm.featureCounts.bam > header.sam

# 2. split alignments by a deterministic hash of the read name (field 1). Same name -> same chunk, so
#    every read's alignments stay together (required for the per-name grouping the gawk relies on).
samtools view mm.featureCounts.bam | awk -v N="$N" '
  { h=0; s=$1; for(i=1;i<=length(s);i++){ h=(h*131 + index("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ:_-.", substr(s,i,1))) % 2000003 } ; print > ("chunk." (h%N) ".sam") }'

# 3. process each chunk in its OWN directory (the gawk writes fixed sam_body filenames, so chunks must
#    not share a CWD), in parallel: transcript assignment -> exon tie-break -> merge -> chunk bam.
assign_chunk() {
  c="$1"
  d="iso.$c"; mkdir -p "$d"
  if [ ! -s "chunk.$c.sam" ]; then samtools view -H -b mm.featureCounts.bam > "$d/annotated.bam"; return; fi
  ( cd "$d"
    # transcript pass: name-group (sort by field 1, where the gawk groups on $1) and gawk-split
    sort -k1,1 "../chunk.$c.sam" | gawk -f "$gawk_script"
    if [ -f assigned_reads.sam_body ]; then cat ../header.sam assigned_reads.sam_body | samtools view -b -h > t.assigned.bam; else samtools view -H -b ../mm.featureCounts.bam > t.assigned.bam; fi
    if [ -f ambiguous_reads.sam_body ]; then cat ../header.sam ambiguous_reads.sam_body | samtools view -h -b > t.unassigned.bam; else samtools view -H -b ../mm.featureCounts.bam > t.unassigned.bam; fi
    rm -f assigned_reads.sam_body ambiguous_reads.sam_body
    # exon tie-break on the still-ambiguous reads
    featureCounts -a "$gtf" -o exon.txt -R BAM t.unassigned.bam -T 1 -t exon -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1 -M
    samtools sort -n t.unassigned.bam.featureCounts.bam | samtools view | gawk -f "$gawk_script"
    if [ -f assigned_reads.sam_body ]; then cat ../header.sam assigned_reads.sam_body | samtools view -b -h > e.assigned.bam; else samtools view -H -b t.unassigned.bam.featureCounts.bam > e.assigned.bam; fi
    samtools merge -f annotated.bam t.assigned.bam e.assigned.bam
  )
}
export -f assign_chunk

seq 0 $((N-1)) | xargs -P "$cpus" -I {} bash -c 'set -e; assign_chunk "$1"' _ {}

# 4. merge the per-chunk annotated bams into the final output
parts=""
for c in $(seq 0 $((N-1))); do parts="$parts iso.$c/annotated.bam"; done
samtools merge -@ "$cpus" -f "${sample_id}.multimapped.annotated.bam" $parts
