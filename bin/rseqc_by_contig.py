#!/usr/bin/env python3
"""
Parallel RSeQC read_distribution by reference contig.

read_distribution.py is single-threaded; its per-read feature assignment is position-local, so
splitting the BAM by reference contig, running read_distribution.py on each contig in parallel, and
summing the per-feature Tag_counts is count-identical to running it on the whole BAM (verified:
Total Reads/Tags/Assigned + all 10 groups sum exactly). Total_bases is the genome-wide feature size
from the BED (constant across contigs); Tags/Kb is recomputed from the summed Tag_count.

This runs entirely inside the RSeQC container (pysam + read_distribution.py); the split is done with
pysam (the container has no samtools) and the input may be unsorted, so reads are routed by reference
name without needing an index.

Usage: rseqc_by_contig.py -i <bam> -r <bed> -p <cpus> > out.txt
"""
import sys, os, re, argparse, subprocess, tempfile, shutil
from concurrent.futures import ThreadPoolExecutor
import pysam

GROUPS = ["CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns",
          "TSS_up_1kb", "TSS_up_5kb", "TSS_up_10kb",
          "TES_down_1kb", "TES_down_5kb", "TES_down_10kb"]
SEP = "=" * 69


def parse(text):
    totals = {"Total Reads": 0, "Total Tags": 0, "Total Assigned Tags": 0}
    bases, counts = {}, {}
    for line in text.splitlines():
        m = re.match(r'(Total Reads|Total Tags|Total Assigned Tags)\s+(\d+)', line)
        if m:
            totals[m.group(1)] = int(m.group(2)); continue
        m = re.match(r'(\S+)\s+(\d+)\s+(\d+)\s', line)
        if m and m.group(1) in GROUPS:
            bases[m.group(1)] = int(m.group(2))
            counts[m.group(1)] = int(m.group(3))
    return totals, bases, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', required=True)
    ap.add_argument('-r', required=True)
    ap.add_argument('-p', type=int, default=4)
    a = ap.parse_args()

    tmp = tempfile.mkdtemp()
    try:
        # Split by reference contig (route reads; input may be unsorted so no index is used).
        inf = pysam.AlignmentFile(a.i, "rb", check_sq=False)
        writers, paths = {}, []
        for rec in inf.fetch(until_eof=True):
            c = rec.reference_name
            if c is None:
                continue
            w = writers.get(c)
            if w is None:
                p = os.path.join(tmp, f"part.{len(writers)}.bam")
                w = pysam.AlignmentFile(p, "wb", template=inf)
                writers[c] = w
                paths.append(p)
            w.write(rec)
        inf.close()
        for w in writers.values():
            w.close()

        if not paths:
            # No mapped reads: defer to read_distribution.py on the whole input (handles it correctly).
            subprocess.run(["read_distribution.py", "-i", a.i, "-r", a.r], check=True)
            return

        def run(p):
            # Fail LOUD: a non-zero exit on any contig would otherwise be summed as zeros and silently
            # under-count the totals.
            r = subprocess.run(["read_distribution.py", "-i", p, "-r", a.r],
                               capture_output=True, text=True)
            if r.returncode != 0:
                sys.stderr.write(r.stderr)
                raise RuntimeError(f"read_distribution.py failed on {p} (exit {r.returncode})")
            return r.stdout

        with ThreadPoolExecutor(max_workers=max(1, a.p)) as ex:
            outputs = list(ex.map(run, paths))

        sum_totals = {"Total Reads": 0, "Total Tags": 0, "Total Assigned Tags": 0}
        sum_counts = {g: 0 for g in GROUPS}
        bases = {}
        for out in outputs:
            t, b, c = parse(out)
            for k in sum_totals:
                sum_totals[k] += t.get(k, 0)
            for g in GROUPS:
                sum_counts[g] += c.get(g, 0)
            bases.update(b)

        o = sys.stdout
        o.write("%-30s%d\n" % ("Total Reads", sum_totals["Total Reads"]))
        o.write("%-30s%d\n" % ("Total Tags", sum_totals["Total Tags"]))
        o.write("%-30s%d\n" % ("Total Assigned Tags", sum_totals["Total Assigned Tags"]))
        o.write(SEP + "\n")
        o.write("%-20s%-20s%-20s%-20s\n" % ("Group", "Total_bases", "Tag_count", "Tags/Kb"))
        for g in GROUPS:
            b = bases.get(g, 0)
            cnt = sum_counts[g]
            # Exactly read_distribution.py's formula (note the +1 on bases) for byte-identical Tags/Kb.
            kb = cnt * 1000.0 / (b + 1)
            o.write("%-20s%-20d%-20d%-18.2f\n" % (g, b, cnt, kb))
        o.write(SEP + "\n")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
