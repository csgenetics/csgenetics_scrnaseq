
# gawk -v r2="R2.fastq.gz" -v sample_id=${sample_id} -v bc_length=13 -f io_extract_all_outputs_interleaved_no_modulo.awk barcode_list.txt <(zcat R1.fastq.gz)

{
	if(ARGIND==1){ # Parse through the barcode list and store in the barcode_A array
		barcode_A[$0]=1;
		next;
	}
	if(ARGIND==2){
		# Start to parse through the fastqs.
		# HEADER
		# We want to get the header from R1 and store it
		# We read in and discard the R2 header
		r1_header_part_1=$1;
		r1_header_part_2=$2;
		"zcat " r2 | getline;
		r2_header=$0;
		
		# READ
		# NB we make a call to getline, because the sequence of the lines
		# is predictable and this saves us having to do modulo checks
		# to see what sort of line we are on.
		getline;

		# This is the sequence line and we need to check the barcode
		# against the barcode list array barcode_A.
		r1_read = $0
		"zcat " r2 | getline;
		r2_read = $0;
		bc=substr($0,1,bc_length);

		if(bc in barcode_A){
			# If the barcode is found in the barcode list, we set good to 1 showing
			# that it is a good read which will mean we evenutally write it out to
			# the good.fastq.gz file.
			good=1;
		}else{
			# If the barcode is not found in the barcode list, we set good to 0 showing
			# that it is a bad read which will be written to bad.fastq.gz.
			good=0;

		}
		
		# THIRD line
		getline;

		# This is the + line. We simply store it at this stage
		r1_third = $0;
		"zcat " r2 | getline;
		r2_third = $0;
		
		# QUAL
		getline;
		
		# This is the qual string and therefore the last of the four lines for the read
		# At this stage we write out to either the good or bad fastq depending
		# on the value of good
		r1_qual = $0;
		"zcat " r2 | getline;
		r2_qual = $0;
		if (good==1){
			good_count++;
			# Then this is a good read and we write out to the good.fastq.gz
			# With the modified header of read 1
			r1_new_header =r1_header_part_1 "_" bc "_ " r1_header_part_2;
			
			printf "%s\n%s\n%s\n%s\n",r1_new_header,r1_read,r1_third,r1_qual | "gzip -1 > " sample_id ".io_extract.R1.fastq.gz";
		}
	}
}
END{
	print "Reads output: "good_count > sample_id ".io_extract.log"
}
