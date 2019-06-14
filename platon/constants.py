
import re

MIN_PROTEIN_IDENTITY = 90.0
MIN_PROTEIN_SCORE_THRESHOLD = -6.8  # sensitivity => 99 %
PROTEIN_SCORE_PENALTY = -1.0
PROTEIN_SCORE_CONSERVATIVE_THRESHOLD = -0.9  # specificity >= 95 %
PROTEIN_SCORE_TRUSTED_THRESHOLD = -0.5  # specificity >= 99.99 %

SPADES_CONTIG_PATTERN = re.compile( 'NODE_\d+_length_\d+_cov_(\d+\.\d+)' )
HEADER = 'ID\tLength\tCoverage\t# ORFs\tProtein Score\tCircular\tInc Type(s)\t# Replication\t# Mobilization\t# Conjugation\t# AMRs\t# rRNAs\t# Plasmid Hits'