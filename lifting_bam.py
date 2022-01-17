import logging
import re
from pathlib import Path
from typing import Annotated, Tuple  # type: ignore

import pysam  # type: ignore
from pydantic import FilePath, constr, validate_arguments

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Liftover")

NAME_REGEX = re.compile("^chr[0-9A-Za-z]+:[0-9]+-[0-9]+$|^[0-9]+:[0-9]+-[0-9]+$")
BamFileName: Annotated[str, constr] = constr(regex=r".*.bam$")
LocusString: Annotated[str, constr] = constr(regex=NAME_REGEX.pattern)


@validate_arguments
def make_ref_fasta(
    reference_fasta: FilePath, chrom: str, start: int, stop: int, padding: int = 5
) -> str:
    """
    extract sequence from reference fasta file

    :param str reference_fasta: file path to the reference fasta file
    :param str chrom: chromosome of the extracted locus
    :param int start: start position of the extracted locus
    :param int stop: stop position of the extracted locus
    :param int padding: padding position to add in front of start and after stop
    """
    if (start - padding) < 0:
        raise ValueError(
            f"start ({start}) must be larger than or equal to padding ({padding})"
        )

    if stop <= start:
        raise ValueError(f"stop ({stop}) must be larger than start {start}")

    padded_start = start - padding
    padded_stop = stop + padding
    fa = pysam.Fastafile(reference_fasta)
    if padded_stop > fa.get_reference_length(chrom):
        raise ValueError(
            f"stop ({stop}) + padding ({padding}) must be smaller than contig size ({chrom})"
        )

    sequence = fa.fetch(chrom, padded_start, padded_stop)
    fa.close()
    return f">{chrom}:{padded_start}-{padded_stop}\n{sequence}"


@validate_arguments
def parse_locus(locus_string: LocusString) -> Tuple[str, int, int]:
    """
    parsing a string to extract chromosome, start and end

    :param str ref_name: reference name must be in the pattern of "chromosom:start-end" (e.g. chr1:100-1000)
    :return: (chrom, start, end)
    :rtype: Tuple(str, int, int)
    """
    locus_string = str(locus_string)
    chrom, coordinates = locus_string.split(":")
    start, end = coordinates.split("-")
    return chrom, int(start), int(end)


def liftover_alignment(header, in_alignment):  # type: ignore
    """
    liftover a pysam alignedsegment
    we only support lifting the whole segment by changing the chromosome name and shifting the start position to the right
    lets say if a alignment is mapping at a contig "chr1:100-1000" at position 2, then we will lift that alignment to chr1:102
    if an alignment is not aligned, it will write out as it is
    :param pysam.libcalignmentfile.AlignmentHeader header: the header of the target genome
    :param pysam.libcalignedsegment.AlignedSegment in_alignment: the input alignment that needed to be lifted
    :return: lifted alignment
    :rtype: pysam.libcalignedsegment.AlignedSegment
    """
    if in_alignment.is_unmapped:
        return in_alignment

    chrom, start, end = parse_locus(in_alignment.reference_name)
    lifted_aln = pysam.AlignedSegment(header)
    lifted_aln.reference_name = chrom
    lifted_aln.query_name = in_alignment.query_name
    lifted_aln.query_sequence = in_alignment.seq
    lifted_aln.qual = in_alignment.qual
    lifted_aln.reference_start = in_alignment.reference_start + start
    lifted_aln.cigar = in_alignment.cigar
    lifted_aln.flag = in_alignment.flag
    lifted_aln.mapping_quality = in_alignment.mapping_quality
    lifted_aln.template_length = in_alignment.template_length
    for tag_id, tag_value in in_alignment.tags:
        lifted_aln.set_tag(tag_id, tag_value)

    # for paired end
    if in_alignment.next_reference_name:
        mate_chrom, mate_subseq_start, mate_subseq_end = parse_locus(
            in_alignment.next_reference_name
        )
        lifted_aln.next_reference_start = (
            in_alignment.next_reference_start + mate_subseq_start
        )
        lifted_aln.next_reference_name = mate_chrom
    return lifted_aln


@validate_arguments
def liftover(gene_bam: FilePath, genome_bam: FilePath, out_bam: BamFileName) -> None:
    """
    Lifting alignments mapping to gene segments to whole genome

    we did alignments to only the gene segment (e.g. chr1:100-1000), and the alignments in this gene_bam file
    would have position 1 that actually correspond to position 100 in the whole genome
    as such, we should shift these alignments to the corresponding positions, so they show up correctly
    in IGV or variant callings.

    The input alignments should map to a reference with the names in the pattern of "chromosom:start-end"

    :param str gene_bam: bam file path storing alignments mapping to the gene segment
    :param str genome_bam: bam file path storing genome alignments
    :param str out_bam: output bam file path to write the lifted over alignments to
    """
    with pysam.AlignmentFile(genome_bam) as genome_alignments, pysam.AlignmentFile(
        gene_bam
    ) as gene_alignments:
        with pysam.AlignmentFile(
            out_bam, "wb", template=genome_alignments
        ) as out_bam_fh:
            aln_count = 0
            for aln in gene_alignments:
                if not aln.is_unmapped or (
                    aln.is_paired and not aln.mate_is_unmapped
                ):  # keeping unmapped paired end alignments
                    lifted_aln = liftover_alignment(genome_alignments.header, aln)
                    out_bam_fh.write(lifted_aln)
                    aln_count += 1
    logger.info(f"Lifted {aln_count} alignments to {out_bam}")


if __name__ == "__main__":
    liftover(Path("alignments/ref1.bam"), Path("alignments/ref2.bam"), "out.bam")
