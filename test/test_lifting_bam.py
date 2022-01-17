import os
import sys
from unittest.mock import MagicMock, patch

import pytest
from conftest import (
    PysamFakeBam,
    PysamFakeFasta,
    mock_alignment,
    mock_bam_header,
    mock_subseq_alignment,
)
from pydantic import ValidationError

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from lifting_bam import liftover, liftover_alignment, make_ref_fasta, parse_locus


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
    with patch("lifting_bam.pysam.Fastafile", return_value=PysamFakeFasta(seqs)):
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
        "lifting_bam.pysam.Fastafile", return_value=PysamFakeFasta(seqs)
    ), pytest.raises(ValueError) as e:
        make_ref_fasta(__file__, chrom=chrom, start=start, stop=stop, padding=pad)

    assert expected_error_message in str(e), f"Fail to capture {test_case}"


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
    with pytest.raises(ValidationError) as e:
        parse_locus(locus_string)
    assert f"string does not match regex" in str(e), f"Didn't catch {test_case}"


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
    alignment = mock_subseq_alignment(
        subseq_chrom,
        subseq_start,
        start,
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


def test_liftover(mock_header):
    gene_bam_header = mock_bam_header([("chr1:1-100", 100)])
    in_gene_alignment = mock_alignment(
        header=gene_bam_header,
        reference_name="chr1:1-100",
        query_name="aln1",
        query_sequence="A" * 10,
        reference_start=10,
        cigar=[(0, 10)],
        flag=0,
        mapping_quality=30,
    )

    out_genome_alignment = mock_alignment(
        header=mock_header,
        reference_name="chr1",
        query_name="aln1",
        query_sequence="A" * 10,
        reference_start=11,
        cigar=[(0, 10)],
        flag=0,
        mapping_quality=30,
    )

    with patch("pysam.AlignmentFile") as pysam_bam:
        mock_in_gene_bam = PysamFakeBam(gene_bam_header, [in_gene_alignment])
        mock_in_genome_bam = PysamFakeBam(mock_header, [])
        mock_out_bam = MagicMock()
        pysam_bam.return_value.__enter__.side_effect = [
            mock_in_genome_bam,
            mock_in_gene_bam,
            mock_out_bam,
        ]
        liftover(__file__, __file__, "/path/to/outbam")
        mock_out_bam.write.assert_called_once_with(out_genome_alignment)


def test_liftover_bad_name():
    with pytest.raises(ValidationError) as e:
        liftover(__file__, __file__, "out.sam")
    assert f"string does not match regex" in str(e)
