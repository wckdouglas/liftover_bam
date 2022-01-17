from collections import OrderedDict
from pathlib import Path
from typing import Dict

import pysam
import pytest

TEST_DIRECTORY = Path(__file__).parent


class PysamFakeBam:
    def __init__(self, header, reads):
        """
        a mock object that mimics the pysam.AlignmentFile object
        :param pysam.AlignmentHeader header: header of the mock sam file
        :param List[pysam.AlignedSegment] reads: reads of the mock sam file
        """
        self.header = header
        self.reads = reads

    def __iter__(self):
        return iter(self.reads)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self

    def close(self):
        return self


class PysamFakeFasta:
    def __init__(self, seqs: Dict):
        self.seqs = seqs

    def fetch(self, chrom, start, end):
        return self.seqs[chrom][start:end]

    def get_reference_length(self, chrom):
        return len(self.seqs[chrom])

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self

    def close(self):
        return self


def mock_bam_header(contig_list):
    """
    making a mock pysam.AlignmentHeader object
    Example::
        contigs = [("chr1", 10), ("chr2", 20)]
        mock_header = mock_bam_header(contigs)
    :param List[Tuple[str, int]] contig_list: a list of tuples of (contig name, contig length)
    :return: a pysam.AlignmentHeader object
    :rtype: pysam.AlignmentHeader
    """
    header_dict = OrderedDict(
        [
            ("SQ", [dict(SN=contig[0], LN=contig[1]) for contig in contig_list]),
        ]
    )
    return pysam.AlignmentHeader.from_dict(header_dict)


def mock_alignment(
    header,
    reference_name,
    query_name,
    query_sequence,
    reference_start,
    cigar,
    flag,
    mapping_quality,
    next_reference_name=None,
    next_reference_start=None,
):
    """
    making a mock pysam.AlignedSegment object
    :param pysam.AlignmentHeader header_dict: a pysam alignment header object (can be created by mock_bam_header)
    :param str reference_name: reference name
    :param str query_name: query name
    :param str query_sequence: query sequence
    :param int reference_start: reference start
    :param list cigar: cigar
    :param int flag: flag
    :param int mapping_quality: mapping quality
    :param str next_reference_name: reference name for the paired end alignment mapped
    :param int next_reference_start: reference start of the paired end alignment
    """
    alignment = pysam.AlignedSegment(header)
    alignment.reference_name = reference_name
    alignment.query_name = query_name
    alignment.query_sequence = query_sequence
    alignment.reference_start = reference_start
    alignment.cigar = cigar
    alignment.flag = flag
    alignment.mapping_quality = mapping_quality
    if (
        next_reference_name is not None
        and next_reference_start is not None
        and next_reference_start > 0
    ):
        alignment.next_reference_name = next_reference_name
        alignment.next_reference_start = next_reference_start
    return alignment


def mock_subseq_alignment(
    subseq_chrom,
    subseq_start,
    start,
    mate_reference_start=None,
):
    subseq_coordinate = f"{subseq_chrom}:{subseq_start}-{subseq_start+100}"
    header_dict = OrderedDict(
        [
            (
                "SQ",
                [
                    {"SN": subseq_coordinate, "LN": subseq_start + 100},
                ],
            ),
        ]
    )
    header = pysam.AlignmentHeader.from_dict(header_dict)
    return mock_alignment(
        header=header,
        reference_name=subseq_coordinate,
        query_name="aln1",
        query_sequence="A" * 10,
        reference_start=start,
        cigar=[(0, 10)],
        flag=0,
        mapping_quality=60,
        next_reference_name=subseq_coordinate
        if mate_reference_start is not None
        else mate_reference_start,
        next_reference_start=mate_reference_start,
    )


@pytest.fixture(scope="module")
def mock_header():
    return mock_bam_header([("chr1", 34), ("chr2", 33), ("chr3", 76)])
