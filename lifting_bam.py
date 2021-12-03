import re
import pysam
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Liftover")

name_regex = re.compile("^chr[0-9A-Za-z]+:[0-9]+-[0-9]+$|^[0-9A-Za-z]+:[0-9]+-[0-9]+$")


def parse_ref_name(ref_name):
    """
    :param str ref_name: reference name must be in the pattern of "chromosom:start-end" (e.g. chr1:100-1000)
    :return: (chrom, start, end)
    :rtype: Tuple(str, int, int)
    """
    if not name_regex.search(ref_name):
        raise ValueError(
            f"reference name is not in the pattern of {name_regex.pattern}"
        )
    else:
        chrom, coordinates = ref_name.split(":")
        start, end = coordinates.split("-")
        return chrom, int(start), int(end)


def liftover_alignment(header, in_alignment):
    """
    :param pysam.libcalignmentfile.AlignmentHeader header: the header of the target genome
    :param pysam.libcalignedsegment.AlignedSegment in_alignment: the input alignment that needed to be lifted
    """
    chrom, start, end = parse_ref_name(in_alignment.reference_name)
    lifted_aln = pysam.AlignedSegment(header)
    lifted_aln.reference_name = chrom
    lifted_aln.query_name = in_alignment.query_name
    lifted_aln.query_sequence = in_alignment.seq
    lifted_aln.reference_start = in_alignment.reference_start + start
    lifted_aln.cigar = in_alignment.cigar
    lifted_aln.flag = in_alignment.flag
    lifted_aln.mapping_quality = in_alignment.mapping_quality
    for tag_id, tag_value in in_alignment.tags:
        lifted_aln.set_tag(tag_id, tag_value)
    return lifted_aln


def liftover(gene_bam, genome_bam, out_bam):
    """
    Lifting alignments mapping to gene segments to whole genome

    we did alignments to only the gene segment (e.g. chr1:100-1000), and the alignments in this gene_bam file
    would have position 1 that actually correspond to position 100 in the whole genome
    as such, we should shift these alignments to the corresponding positions, so they show up correctly
    in IGV or variant callings.

    The input alignments should map to a reference with the names in the pattern of "chromosom:start-end"

    :param str gene_bam: bam file storing alignments mapping to the gene segment
    :param str genome_bam: bam file storing genome alignments
    :param str out_bam: output bam file with lifted over alignments
    """
    with pysam.AlignmentFile(genome_bam) as genome_alignments, pysam.AlignmentFile(
        gene_bam
    ) as gene_alignments:
        with pysam.AlignmentFile(
            "out.bam", "wb", template=genome_alignments
        ) as out_bam_fh:
            for aln_count, aln in enumerate(gene_alignments):
                lifted_aln = liftover_alignment(genome_alignments.header, aln)
                out_bam_fh.write(lifted_aln)
    logger.info(f"Lifted {aln_count+1} to {out_bam}")


if __name__ == "__main__":
    liftover("alignments/ref1.bam", "alignments/ref2.bam", "out.bam")
