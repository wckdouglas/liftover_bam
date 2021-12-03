all: alignments/ref1.bam alignments/ref2.bam

ref/ref1.fa.amb:
	bwa index ref/ref1.fa

ref/ref2.fa.amb:
	bwa index ref/ref2.fa

alignments/ref1.bam: ref/ref1.fa.amb
	bwa mem -A 1 -T 5 -k4  ref/ref1.fa query.fa | samtools view -b > alignments/ref1.bam

alignments/ref2.bam: ref/ref2.fa.amb
	bwa mem -A 1 -T 5 -k4  ref/ref2.fa query.fa | samtools view -b > alignments/ref2.bam
