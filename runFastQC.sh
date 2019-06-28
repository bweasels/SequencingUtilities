module load fastqc

for f in ./*/*.fastq.gz;
	do fastqc $f;
done
