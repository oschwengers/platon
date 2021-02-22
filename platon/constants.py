
import re

MIN_CONTIG_LENGTH = 1000
MAX_CONTIG_LENGTH = 500000
MIN_PROTEIN_IDENTITY = 90.0
RDS_SENSITIVITY_THRESHOLD = -7.9  # sensitivity => 95 %
RDS_CONSERVATIVE_THRESHOLD = 0.1  # highest accuracy
RDS_SPECIFICITY_THRESHOLD = 0.7  # specificity >= 99.9 %
MIN_CIRC_BASEPAIR_OVERLAP = 100

SPADES_CONTIG_PATTERN = re.compile(r'NODE_\d+_length_\d+_cov_(\d+\.\d+)')
UNICYCLER_CONTIG_PATTERN = re.compile(r'\d+ length=\d+ depth=(\d+\.\d{2})x( circular=true)?')
HEADER = 'ID\tLength\tCoverage\t# ORFs\tRDS\tCircular\tInc Type(s)\t# Replication\t# Mobilization\t# OriT\t# Conjugation\t# AMRs\t# rRNAs\t# Plasmid Hits'
CITATION = '''Schwengers O., Barth P., Falgenhauer L., Hain T., Chakraborty T., & Goesmann A. (2020).
Platon: identification and characterization of bacterial plasmid contigs in short-read draft assemblies exploiting protein sequence-based replicon distribution scores.
Microbial Genomics, 95, 295. https://doi.org/10.1099/mgen.0.000398'''
