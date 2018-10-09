
# Conjugation proteins
grep Tra[^ICG] $1 | cut -f2,3,4,8 >> conjugation.tmp.tsv
grep Trb[A-Z] $1 | cut -f2,3,4,8 >> conjugation.tmp.tsv
grep Trw[^ABC] $1 | cut -f2,3,4,8 >> conjugation.tmp.tsv
grep VirB[0-9] $1 | cut -f2,3,4,8 >> conjugation.tmp.tsv

sort conjugation.tmp.tsv | uniq > conjugation.tsv
rm conjugation.tmp.tsv
