import sys
from pathlib import Path

from liftover_bam import make_ref_fasta

print(make_ref_fasta(Path(sys.argv[1]), chrom="chr3", start=29, stop=46, padding=2))
