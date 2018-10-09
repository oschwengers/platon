
# Relaxases
grep MobC $1 | cut -f2,3,4,8 >> mobilization.tmp.tsv
grep TrwC $1 | cut -f2,3,4,8 >> mobilization.tmp.tsv
grep TraI $1 | cut -f2,3,4,8 >> mobilization.tmp.tsv

# Relaxosome proteins
grep -i plasmid $1 | grep -i mobilization | cut -f2,3,4,8 >> mobilization.tmp.tsv
grep -i relaxosome $1 | cut -f2,3,4,8 >> mobilization.tmp.tsv

# T4CP proteins
grep TraD $1 | cut -f2,3,4,8 >> mobilization.tmp.tsv
grep VirD4 $1 | cut -f2,3,4,8 >> mobilization.tmp.tsv
grep TrwB $1 | cut -f2,3,4,8 >> mobilization.tmp.tsv
grep -i conjugative $1 | grep -i coupling | cut -f2,3,4,8 >> mobilization.tmp.tsv

sort mobilization.tmp.tsv | uniq > mobilization.tsv
rm mobilization.tmp.tsv
