#!/bin/bash

mkdir db
cd db

# download incompatibility group signatures from CGE PlasmidFinder DB
printf "1/10: download incompatibility group signatures from CGE PlasmidFinder DB...\n"
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git
grep '\S' plasmidfinder_db/gram_positive.fsa > inc-types-raw.fasta
grep '\S' plasmidfinder_db/enterobacteriaceae.fsa >> inc-types-raw.fasta
cat inc-types-raw.fasta | sed -r "s/^>([a-zA-Z0-9()]+)_([0-9]+)_(.*)_(.*)$/>\1 \4/" > inc-types-mixed.fasta
awk '{if(!/>/){print toupper($0)}else{print $1}}' inc-types-mixed.fasta > inc-types.fasta
rm -rf plasmidfinder_db inc-types-raw.fasta inc-types-mixed.fasta


# download rRNA covariance models from Rfam
printf "\n2/10: download rRNA covariance models from Rfam...\n"
mkdir Rfam
cd Rfam
wget -q -nH ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.tar.gz
tar -xzf Rfam.tar.gz
cd ..
cat Rfam/RF00001.cm >  rRNA
cat Rfam/RF00177.cm >> rRNA
cat Rfam/RF02541.cm >> rRNA
cmpress rRNA
rm -r Rfam rRNA


# download AMR HMMs from NCBIfams
printf "\n3/10: download AMR HMMs from NCBIfams...\n"
wget -q -nH ftp://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMR/latest/NCBIfam-AMR.LIB
mv NCBIfam-AMR.LIB ncbifam-amr
hmmpress ncbifam-amr
rm ncbifam-amr
wget -q -nH ftp://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMR/latest/NCBIfam-AMR.tsv
mv NCBIfam-AMR.tsv ncbifam-amr.tsv


# download RefSeq reference plasmids
printf "\n4/10: download RefSeq reference plasmids...\n"
wget -q -O refseq-plasmids-raw.tsv ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/plasmids.txt
grep 'Bacteria' refseq-plasmids-raw.tsv | cut -f1,5,6,8,9,10,11,12,15 > refseq-plasmids.tsv
grep 'Bacteria' refseq-plasmids-raw.tsv | cut -f3 refseq-plasmids.tsv > refseq-plasmids-ids.txt
mkdir refseq-plasmids-dir
cd refseq-plasmids-dir
wget -q -nH ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.*.1.genomic.fna.gz
cd ..
gzip -dc refseq-plasmids-dir/plasmid.*.1.genomic.fna.gz | seqtk subseq - refseq-plasmids-ids.txt | seqtk seq -CU > refseq-plasmids
makeblastdb -dbtype nucl -in refseq-plasmids -title 'RefSeq Plasmids'
mv refseq-plasmids refseq-plasmids.fna
rm -r refseq-plasmids-dir refseq-plasmids-raw.tsv refseq-plasmids-ids.txt


# download RefSeq nonredundant proteins and clusters
printf "\n5/10: download RefSeq nonredundant proteins and clusters...\n"
wget -q -nH ftp://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_proteins.txt
mkdir refseq-bacteria
cd refseq-bacteria
wget -q -nH ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.nonredundant_protein.*.protein.faa.gz
cd ..
gzip -dc refseq-bacteria/bacteria.nonredundant_protein.* | seqtk seq -CU > refseq-bacteria-nrp.trimmed.faa
rm -r refseq-bacteria


# extract NRP cluster representatives and build db
printf "\n6/10: extract NRP cluster representatives and build database...\n"
python3 $PLATON_HOME/db-scripts/select-cluster-rep.py --cluster PCLA_proteins.txt --proteins refseq-bacteria-nrp.trimmed.faa > refseq-bacteria-nrpc-reps.txt
seqtk subseq refseq-bacteria-nrp.trimmed.faa refseq-bacteria-nrpc-reps.txt > refseq-bacteria-nrpc-reps.faa
ghostz db -i refseq-bacteria-nrpc-reps.faa -o refseq-bacteria-nrpc-reps -L 6
rm refseq-bacteria-nrpc-reps.txt refseq-bacteria-nrpc-reps.faa


# build HMMs
printf "\n7/10: build HMMs...\n"
wget -q -nH ftp://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_clusters.txt
touch clusters.lst
for type in replication mobilization conjugation; do
    sh $PLATON_HOME/db-scripts/extract-${type}.sh PCLA_clusters.txt
    cut -f1 ${type}.tsv >> clusters.lst
done
grep -f clusters.lst PCLA_proteins.txt | cut -f2 > clusters-proteins.lst
seqtk subseq refseq-bacteria-nrp.trimmed.faa clusters-proteins.lst > clusters-proteins.faa
for type in replication mobilization conjugation; do
    printf "\nbuild ${type}...\n"
    nextflow run $PLATON_HOME/db-scripts/build-hmms.nf \
        --nrpc ${type}.tsv \
        --nrpcMapping PCLA_proteins.txt \
        --nrp clusters-proteins.faa \
        --geneClass ${type}
    mv ${type}.hmm ${type}
    hmmpress ${type}
    rm -rf work .nextflow* ${type}.tsv ${type}
done
rm refseq-bacteria-nrp.trimmed.faa PCLA_clusters.txt clusters-proteins.faa clusters.lst clusters-proteins.lst


# download RefSeq chromosomes
printf "\n8/10: download RefSeq chromosomes...\n"
mkdir nf-tmp
nextflow run $PLATON_HOME/db-scripts/download-chromosomes.nf
rm -rf work .nextflow* nf-tmp


# calculate protein scores
printf "\n9/10: calculate protein scores...\n"
nextflow run $PLATON_HOME/db-scripts/calculate-scores.nf \
    --plasmids ./refseq-plasmids.fna \
    --chromosomes ./refseq-chromosomes.fna \
    --protClusterMapping ./PCLA_proteins.txt \
    --nrpcDB ./refseq-bacteria-nrpc-reps.inf
rm -rf work .nextflow* PCLA_proteins.txt
cut -f 1,2,3,10 rds-full.tsv > rds.tsv


# calculate protein score cutoffs
printf "\n10/10: calculate protein score cutoffs...\n"
export NXF_OPTS="-Xms32G -Xmx64G"
nextflow run $PLATON_HOME/db-scripts/test-scores.nf \
    --plasmids ./refseq-plasmids.fna \
    --chromosomes ./refseq-chromosomes.fna \
    --protein_scores ./rds.tsv \
    --nrpcDB ./refseq-bacteria-nrpc-reps.inf
rm -rf work .nextflow* refseq-plasmids.fna refseq-chromosomes.fna

cd ..
mv db $PLATON_HOME/
