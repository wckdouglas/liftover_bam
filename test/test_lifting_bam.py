import os
import sys
from collections import OrderedDict
from typing import Dict
from unittest.mock import patch

import pysam
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from lifting_bam import NAME_REGEX, liftover_alignment, make_ref_fasta, parse_locus


class FakeFastaFile:
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


@pytest.mark.parametrize(
    "seqs,chrom,start,stop,pad,expected_name,expected_seq",
    [
        (
            {"chr1": "ACTGACTGACTG", "chr2": "TTTGGGTTTGGG"},
            "chr1",
            2,
            4,
            2,
            "chr1:0-6",
            "ACTGAC",
        ),
        (
            {"chr1": "ACTGACTGACTG", "chr2": "TTTGGGTTTGGGTGC"},
            "chr2",
            5,
            10,
            4,
            "chr2:1-14",
            "TTGGGTTTGGGTG",
        ),
    ],
)
def test_make_ref_fasta(seqs, chrom, start, stop, pad, expected_name, expected_seq):
    expected_out = f">{expected_name}\n{expected_seq}"
    with patch("lifting_bam.pysam.Fastafile", return_value=FakeFastaFile(seqs)):
        assert (
            make_ref_fasta(__file__, chrom=chrom, start=start, stop=stop, padding=pad)
            == expected_out
        )


@pytest.mark.parametrize(
    "test_case, seqs,chrom,start,stop,pad,expected_error_message",
    [
        (
            "start - padding < 0",
            {"chr1": "ACTGACTGACTG", "chr2": "TTTGGGTTTGGG"},
            "chr1",
            2,
            4,
            3,
            "must be larger than or equal to padding",
        ),
        (
            "stop <= start",
            {"chr1": "ACTGACTGACTG", "chr2": "TTTGGGTTTGGGTGC"},
            "chr2",
            10,
            5,
            6,
            "must be larger than start",
        ),
        (
            "padding + stop > reference size",
            {"chr1": "ACTGACTGACTG", "chr2": "TTTGGGTTTGGGTGC"},
            "chr2",
            9,
            10,
            7,
            "must be smaller than contig size",
        ),
    ],
)
def test_make_ref_fasta_error(
    test_case, seqs, chrom, start, stop, pad, expected_error_message
):
    with patch(
        "lifting_bam.pysam.Fastafile", return_value=FakeFastaFile(seqs)
    ), pytest.raises(ValueError) as e:
        make_ref_fasta(__file__, chrom=chrom, start=start, stop=stop, padding=pad)

    assert expected_error_message in str(e), f"Fail to capture {test_case}"


@pytest.fixture(scope="module")
def mock_header():
    header_dict = OrderedDict(
        [
            (
                "SQ",
                [
                    {"SN": "chr1", "LN": 34},
                    {"SN": "chr2", "LN": 33},
                    {"SN": "chr3", "LN": 76},
                ],
            ),
            (
                "PG",
                [
                    {
                        "ID": "bwa",
                        "PN": "bwa",
                        "VN": "0.7.17-r1198-dirty",
                        "CL": "bwa mem -A 1 -T 5 -k4 ref/ref2.fa query.fa",
                    },
                    {
                        "ID": "samtools",
                        "PN": "samtools",
                        "PP": "bwa",
                        "VN": "1.13",
                        "CL": "samtools view -b",
                    },
                ],
            ),
        ]
    )
    return pysam.AlignmentHeader.from_dict(header_dict)


@pytest.mark.parametrize(
    "chrom,start,end",
    [("chr1", 10, 100), ("chr2", 1, 200), ("chr2", 0, 12)],
)
def test_parse_locus(chrom, start, end):
    out_chrom, out_start, out_end = parse_locus(f"{chrom}:{start}-{end}")
    assert out_chrom == chrom
    assert out_start == start
    assert out_end == end


@pytest.mark.parametrize(
    "test_case, locus_string",
    [
        ("Non standard chromosome name", "a:1-10"),
        ("Non-integer start", "chr1:1.1-100"),
        ("Non-integer end", "chr1:1-1.1"),
        ("Gene name", "TP53"),
    ],
)
def test_parse_locus_bad_input(test_case, locus_string):
    with pytest.raises(ValueError) as e:
        parse_locus(locus_string)
    assert (
        f"reference name is not in the pattern of {NAME_REGEX.pattern}: {locus_string}"
        in str(e)
    ), f"Didn't catch {test_case}"


def make_alignment(
    subseq_chrom,
    subseq_start,
    start,
    header,
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
    mock_alignment = pysam.AlignedSegment(header)
    mock_alignment.reference_name = subseq_coordinate
    mock_alignment.query_name = "aln1"
    mock_alignment.query_sequence = "A" * 10
    mock_alignment.reference_start = start
    mock_alignment.cigar = [(0, 10)]
    mock_alignment.flag = 0
    mock_alignment.mapping_quality = 60
    if mate_reference_start is not None:
        mock_alignment.next_reference_name = subseq_coordinate
        mock_alignment.next_reference_start = mate_reference_start
    return mock_alignment


@pytest.mark.parametrize(
    "test_case,subseq_chrom, subseq_start, start, mate_reference_start",
    [
        ("Singleton alignment", "chr1", 10, 10, None),
        ("Paired end alignment", "chr1", 10, 10, 12),
    ],
)
def test_liftover_alignment(
    mock_header,
    test_case,
    subseq_chrom,
    subseq_start,
    start,
    mate_reference_start,
):
    alignment = make_alignment(
        subseq_chrom,
        subseq_start,
        start,
        mock_header,
        mate_reference_start=mate_reference_start,
    )
    lifted_alignment = liftover_alignment(mock_header, alignment)
    assert lifted_alignment.reference_name == subseq_chrom, f"Failed for {test_case}"
    assert (
        lifted_alignment.reference_start == start + subseq_start
    ), f"Failed for {test_case}"
    if mate_reference_start is not None:
        assert (
            lifted_alignment.next_reference_start == mate_reference_start + subseq_start
        ), f"Failed for {test_case}"
        assert (
            lifted_alignment.next_reference_name == subseq_chrom
        ), f"Failed for {test_case}"
