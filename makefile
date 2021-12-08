all: alignments/ref1.bam alignments/ref2.bam

ref/ref1.fa:
	python test/make.py ref/ref2.fa > ref/ref1.fa
	

ref/ref1.fa.amb: ref/ref1.fa
	bwa index ref/ref1.fa

ref/ref2.fa.amb:
	bwa index ref/ref2.fa

alignments/ref1.bam: ref/ref1.fa.amb
	bwa mem -A 1 -T 5 -k4  ref/ref1.fa query.fa | samtools view -b > alignments/ref1.bam

alignments/ref2.bam: ref/ref2.fa.amb
	bwa mem -A 1 -T 5 -k4  ref/ref2.fa query.fa | samtools view -b > alignments/ref2.bam

out.bam: alignments/ref1.bam alignments/ref2.bam
	poetry run python lifting_bam.py	

.PHONY: test
test: out.bam
	samtools view -h out.bam | grep -v "bwa\|samtools"  | diff - test/expected.sam

.PHONY: clean
clean: 
	rm ref/ref1* ref/ref2.fa.* alignments/*bam out.bam
