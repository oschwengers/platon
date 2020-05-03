
import re

MIN_CONTIG_LENGTH = 1000
MAX_CONTIG_LENGTH = 500000
MIN_PROTEIN_IDENTITY = 90.0
RDS_SENSITIVITY_THRESHOLD = -7.7  # sensitivity => 95 %
RDS_CONSERVATIVE_THRESHOLD = 0.1  # (-0.9) specificity >= 95 %
RDS_SPECIFICITY_THRESHOLD = 0.4  # (-0.5) specificity >= 99.9 %
MIN_CIRC_BASEPAIR_OVERLAP = 100

SPADES_CONTIG_PATTERN = re.compile(r'NODE_\d+_length_\d+_cov_(\d+\.\d+)')
HEADER = 'ID\tLength\tCoverage\t# ORFs\tRDS\tCircular\tInc Type(s)\t# Replication\t# Mobilization\t# OriT\t# Conjugation\t# AMRs\t# rRNAs\t# Plasmid Hits'
CITATION = '''Schwengers O., Barth P., Falgenhauer L., Hain T., Chakraborty T., Goesmann A. (2020)
Platon: identification and characterization of bacterial plasmid contigs in short-read draft assemblies exploiting protein-sequence-based replicon distribution scores.
bioRxiv 2020.04.21.053082; doi: https://doi.org/10.1101/2020.04.21.053082'''
