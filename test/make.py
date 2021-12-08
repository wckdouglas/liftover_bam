import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from lifting_bam import make_ref_fasta

make_ref_fasta(sys.argv[1], chrom="chr3", sstart=30, stop=46, padding=2)
