
# Replication proteins
grep RepH $1 | cut -f2,3,4,8 >> replication.tmp.tsv
grep -i plasmid $1 | grep -i rep | cut -f2,3,4,8 >> replication.tmp.tsv

# Partitioning proteins

grep SopA $1 | cut -f2,3,4,8 >> replication.tmp.tsv
grep KorB $1 | cut -f2,3,4,8 >> replication.tmp.tsv
grep ParM $1 | cut -f2,3,4,8 >> replication.tmp.tsv
grep ParR $1 | cut -f2,3,4,8 >> replication.tmp.tsv
grep -i plasmid $1 | grep -i partition | cut -f2,3,4,8 >> replication.tmp.tsv

sort replication.tmp.tsv | uniq > replication.tsv
rm replication.tmp.tsv
