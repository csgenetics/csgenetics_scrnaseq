
# Parse a query-sorted bam file that contains multimapping reads that have been run through
# featureCounts to assign to geneic regions (-t gene).

# For those reads that are assigned to only one gene, revert the reads to be non-multipmapping
# (i.e NH:i:1  HI:i:1 ) and output the alignment to a file that can then be used to create a new bam.
# We will take account of the fact that a read may have multiple alignments (assigned) to the
# same gene by using the XT tag as the key to the assigned_array_corrected. However, for the assigned_array_uncorrected
# we do not use the XT tag as key because we want to keep all alignments to put into the ambiguous_reads.sam body.
# Instead we use the read name and the alignment number tag e.g. VH00671:469:AACFC22HV:1:1102:44836:17149_TGGTATGCAAGGC_HI:i:8.

# Because ambigously assigned alignments do not have an XT tag, we also use the read name + alignment number as key.

# The script will output 2 files. The first is a file containing the alignments that
# were assigned unambiguously. E.g. if (length(assigned_array_corrected)==1 && length(unassigned_ambiguity_array)==0).
# The second will contain alignments for reads in which unambigous assignment was not
# possible. E.g. the else statement of the above if statement. (if both arrays are empty, no reads will be output).

# Logic for handling the multipmappers on a read by read (multiple alignments):
#   If Assigned count == 1 -> output to the assigned_reads.sam_body sam_body
#   If Assigned count && Unassigned_Ambiguity count == 0 -> There is nothing to do. Do not output the alignments.
#   If Assigned count + Unassigned_Ambiguity count > 0 -> Then output the alignments that are responsible for the
#   counts (either of the counts) to ambiguous_reads.sam_body to be run through the exon tie breaking.

# We then run the ambigous_reads through this script again after having been through exon annotation using featureCounts.
BEGIN {
    current_query="";
    split("", assigned_array_corrected);
    split("", assigned_array_uncorrected);
    split("", unassigned_ambiguity_array);
    }

{
    if ($1 != current_query){
        # Then we are starting a count for a new query
        # Check to see if there is exactly 1 assigned genes
        # and 0 ambiguous_count
        if (length(assigned_array_corrected)==1 && length(unassigned_ambiguity_array)==0){
            # Then we have found a recoverable assigned read
            # Output the corrected alignment to the assigned_reads
            # NB only contains 1 item
            for (key in assigned_array_corrected){
                printf "%s\n",assigned_array_corrected[key] > "assigned_reads.sam_body";
            }
        }else{
            # Then output all alignments found in assigned_array_uncorrected and unassigned_ambiguity_array
            # N.B. this may be 0 alignments
            for (key in assigned_array_uncorrected){
                printf "%s\n",assigned_array_uncorrected[key] > "ambiguous_reads.sam_body";
            }
            for (key in unassigned_ambiguity_array){
                printf "%s\n",unassigned_ambiguity_array[key] > "ambiguous_reads.sam_body";
            }
        }

        # Update the current query
        current_query = $1;

        # Empty/initiate the arrays
        split("", assigned_array_corrected);
        split("", assigned_array_uncorrected);
        split("", unassigned_ambiguity_array);
    }

    if ($16 == "XS:Z:Assigned"){
        # This is an assigned alignment
        # We add the alignment to the assigned_array_corrected and assigned_array_uncorrected
        # The alignment that has been added assigned_array_corrected has had NH:i:X  HI:i:X set to NH:i:1  HI:i:1
        # Key is the  XT target and value is the full alignment
        # E.g. key = XT:Z:ENSG00000150093
        # value = VH00671:444:AACCJYTHV:1:1101:5696:28489_AAGCACCTATCCG_	256	5	67803287	0	34M	*	0	0	TTATTCACTATCACAAGAATAACACGGGAAAGAC	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	NH:i:X  HI:i:X	AS:i:33	nM:i:0	XS:Z:Assigned	XN:i:1	XT:Z:ENSG00000249364
        assigned_array_corrected[$18] = $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\tNH:i:1\tHI:i:1\t"$14"\t"$15"\t"$16"\t"$17"\t"$18;

        # For the uncorrected array we output without the featureCount derived tags because these alignments
        # will be run through feature counts again.
        # We also use the query_name + alignment number as a key so that all alignments are output
        # (including multiple assigned alignments that were assigned to the same gene (same XT tag)).
        assigned_array_uncorrected[$1 $13] = $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15;
    }

    if ($16 == "XS:Z:Unassigned_Ambiguity"){
        # This is an Unassigned_Ambiguity alignment
        # Key is the query_name + alignment number and value is the full alignment line uncorrected
        # E.g. key = VH00671:469:AACFC22HV:1:1102:44836:17149_TGGTATGCAAGGC_HI:i:8
        # value = VH00671:444:AACCJYTHV:1:1101:5696:28489_AAGCACCTATCCG_	256	5	67803287	0	34M	*	0	0	TTATTCACTATCACAAGAATAACACGGGAAAGAC	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	NH:i:X  HI:i:X	AS:i:33	nM:i:0	XS:Z:Unassigned_Ambiguity
        # We output the alignment without the featureCount derived tags because these alignments
        # will be run through feature counts again.
        unassigned_ambiguity_array[$1 $13] = $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15;
    }
    
}

END {
    # Then we have finished reading through the bam and we need to output for the last read
    if (length(assigned_array_corrected)==1 && length(unassigned_ambiguity_array)==0){
            # Then we have found a recoverable assigned read
            # Output the corrected alignment to the assigned_reads
            # NB only contains 1 item
            for (key in assigned_array_corrected){
                printf "%s\n",assigned_array_corrected[key] > "assigned_reads.sam_body";
            }
        }else{
            # Then output all alignments found in assigned_array_uncorrected and unassigned_ambiguity_array
            # N.B. this may be 0 alignments
            for (key in assigned_array_uncorrected){
                printf "%s\n",assigned_array_uncorrected[key] > "ambiguous_reads.sam_body";
            }
            for (key in unassigned_ambiguity_array){
                printf "%s\n",unassigned_ambiguity_array[key] > "ambiguous_reads.sam_body";
            }
        }
}