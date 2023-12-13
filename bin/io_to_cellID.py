#!/usr/bin/env python
import sys
import csv
import pysam, argparse

# Parse the command line arguments
def parse_arguments(args):
    parser = argparse.ArgumentParser(description = "Arguments for io_to_cellID.py script to transform io_sequence to cellID_ioID")
    parser.add_argument("--sample", help="Sample ID")
    parser.add_argument("--white_list", help="path to white_list.csv file")
    parser.add_argument("--dedup", help="deduplication parameter")

    return parser.parse_args()


def transform_io(sam_file, bam_out_path, wl_dict, to_fastq=True,fastq_out_path=''):
    # Write bam file
    with pysam.AlignmentFile(sam_file,'r') as infile:
        with pysam.AlignmentFile(bam_out_path,"wb", template=infile) as outfile:
            for read in infile:
                readID = read.query_name
                # Convert IO sequence to corresponding cellID_ioID
                io = readID.split('_')[1]
                readID = readID.split('_')[0]
                cellID = wl_dict[io][0]
                ioID = wl_dict[io][1]
                read.query_name = '{readID}_{cellID}_io{ioID}'.format(readID=readID,cellID=cellID,ioID=ioID)
                outfile.write(read)
    pysam.index(bam_out_path)

    # Write fastq file
    if to_fastq:
        if not fastq_out_path:
            fastq_out_path = bam_out_path.replace('.bam','.fastq')

        with pysam.AlignmentFile(bam_out_path,'r') as infile:
            with open(fastq_out_path, "w") as fastq:
                for read in infile:
                    fastq.write("@%s\n" % read.qname)
                    fastq.write("%s\n" % read.seq)
                    fastq.write("+\n")
                    fastq.write("%s\n" % read.qual)

    return()


def launch(args):
    # load CSGX barcode_list of IOs
    wl_dict={}
    with open(args.white_list,'r') as inf:
        reader = csv.reader(inf)
        wl_dict = {rows[1]:[rows[0],rows[2]] for rows in reader}

    sample = args.sample
    dedup = args.dedup
    # UMR_sam (group_filtered.sam)
    umr_sam_path = '{sample}_group_filtered.sam'.format(sample=sample)
    out_umr_bam_path = '{sample}_UMR.bam'.format(sample=sample)
    transform_io(sam_file=umr_sam_path, bam_out_path=out_umr_bam_path,wl_dict=wl_dict)

    if dedup == 'true':
        # dedup_UMR_sam (dedup.sam)
        dumr_sam_path = '{sample}_dedup.sam'.format(sample=sample)
        out_dumr_bam_path = '{sample}_deduplicated_UMR.bam'.format(sample=sample)
        transform_io(sam_file=dumr_sam_path, bam_out_path=out_dumr_bam_path,wl_dict=wl_dict)


if __name__ == "__main__":
    args = parse_arguments(sys.argv[1:])
    launch(args)
