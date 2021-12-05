import os
import sys
from collections import OrderedDict

import pysam
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from lifting_bam import liftover_alignment, parse_locus, NAME_REGEX


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


def make_alignment(subseq_chrom, subseq_start, start, header):
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
    return mock_alignment


@pytest.mark.parametrize(
    "subseq_chrom, subseq_start, start",
    [("chr1", 10, 10), ("chr2", 20, 5), ("chr3", 30, 12)],
)
def test_liftover_alignment(mock_header, subseq_chrom, subseq_start, start):
    alignment = make_alignment(subseq_chrom, subseq_start, start, mock_header)
    lifted_alignment = liftover_alignment(mock_header, alignment)
    assert lifted_alignment.reference_name == subseq_chrom
    assert lifted_alignment.reference_start == start + subseq_start - 1
