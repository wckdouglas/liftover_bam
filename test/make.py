import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from lifting_bam import make_ref_fasta

print(make_ref_fasta(sys.argv[1], chrom="chr3", start=29, stop=46, padding=2))
