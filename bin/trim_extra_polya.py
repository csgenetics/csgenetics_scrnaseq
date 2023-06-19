#!/usr/bin/env python

import sys, os, re, gzip
from itertools import islice
from numpy import mean, median

fastpOut = sys.argv[1]
sample = sys.argv[2]

counter = 0

short_read = 0
good_reads = 0
mid_polyA = 0
polyA_andASP = 0
r1_asp = 0
readLen = list()

# TODO remove read length output is not used.

outf2 = open('read_lengths_post_trimming.csv','a+')
outf = gzip.open(sample + '_R1.polyA_trimmed.fastq.gz','a+')

with gzip.open(fastpOut, 'rb') as f:
    
    while True:
        # read.data contains the 4 lines of data from fastq file for that read (@title + seq + '+' + quality_score)
        r1_data = list(islice(f, 4))        
        r1_data = [x.decode("utf-8") for x in r1_data]
        
        if not r1_data:
            break

        counter += 1
        read_id = r1_data[0].split(" ")[0]                

        polyA = "A"*15
        antislip_p = 'GCACG'
        antislip_p = 'A'*10 + antislip_p
        
        # Deal with polyA
        polya_flag = 0
        search_string_polya = re.compile(str(polyA))
        m = search_string_polya.search(r1_data[1])
        if m:
            polya_flag += 1            
            polya_pos = m.span()[0]
            mid_polyA += 1

            # Trim the read where polyA starts
            r1_data[1] = r1_data[1].rstrip()[:polya_pos] + '\n'
            r1_data[3] = r1_data[3].rstrip()[:polya_pos] + '\n'

        # Deal with PolyA + AS
        as_flag = 0                        
        search_string_asp = re.compile(str(antislip_p))
        m2 = search_string_asp.search(r1_data[1])
        if m2:
            as_flag += 1                                
            asp_start_pos = m2.span()[0]
            asp_end_pos = m2.span()[1]

            # Trim the read where polyA starts
            r1_data[1] = r1_data[1].rstrip()[:asp_start_pos] + '\n'
            r1_data[3] = r1_data[3].rstrip()[:asp_start_pos] + '\n'
            
            r1_asp += 1

            if polya_flag == 1:
                polyA_andASP +=1                        
        
        outf2.write(str(len(r1_data[1]))+'\n')
        
        readLen.append(len(r1_data[1]))

        if len(r1_data[1].strip()) <= 20:            
            short_read +=1
        else:
            good_reads += 1
            outf.write(''.join(r1_data).encode('utf-8'))

outf2.close()
outf.close()

with open(f'{sample}_trim_polyA_metrics.csv','w') as f:                
    f.write(str(good_reads)+'\n'+ str(mean(readLen)) + '\n' + str(median(readLen)))
