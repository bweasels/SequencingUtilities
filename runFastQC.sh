module load fastqc
module load multiqc

for f in ./*/*.fastq.gz;
	do fastqc $f;
done

multiqc ./
