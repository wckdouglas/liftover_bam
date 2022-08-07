# Lifting over bam #

[![poetry CI](https://github.com/wckdouglas/liftover_bam/actions/workflows/ci.yml/badge.svg)](https://github.com/wckdouglas/liftover_bam/actions/workflows/ci.yml)

Sometimes for amplicon sequencings, we would like to map reads to the amplicon sequence only but bringing them back to genomic coordinates for easy variant calling and viewing.

Let's say we have a gene in `chr1:100-1000`, and we would first extract this locus from the genome fasta file to make a new fasta record with name `>chr1:100-1000`, this can be done with:
```
echo "chr1:100-1000" | samtools faidx -r - genome.fa > gene.fa 
```
and map the reads to this single gene fasta file with `bwa` or `bowtie2` to make a bam alignment file:
```
bwa mem gene.fa query.fq | samtools view -b > gene.bam
```

So what if you want to put these alignments back to the genomic coordinates after that?

The `liftover_bam.liftover` function is trying to solve this problem in pure python!

```
gene_bam="gene.bam"
genome_bam="any.bam_file_that_maps_to_the_genome"
out_bam="where_you_want_your_output_bam_file_to_be"
liftover(gene_bam, genome_bam, out_bam)
```
