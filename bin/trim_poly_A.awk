{
	if(NR % 4 == 1){header=$0; next;}
	if(NR % 4 == 2){if(match($0, /[A]{15,}|[A]{13,}CG/)){ type=1; read=substr($0,1,RSTART-1); read_length=length(read); match_count++; next;}else{type=0; read=$0; read_length=length(read); next;}}
	if(NR % 4 == 3){third=$0;next;}
	if(NR % 4 == 0){
	if(type==0){qual=$0;}else{qual=substr($0,1,RSTART-1);}
	if(read_length <= 5){printf "%s\n%s\n%s\n%s\n",header,read,third,qual | "gzip -1 > bad.fastq.gz";}
	else{printf "%s\n%s\n%s\n%s\n",header,read,third,qual | "gzip -1 > good.fastq.gz";}}
}
