samtools depth $1 -b $2 | awk '{sum=sum+$3}END{print sum"\t"NR}' 
