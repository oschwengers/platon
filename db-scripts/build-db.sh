#!/bin/bash

mkdir db
cd db

# download incompatibility group signatures from CGE PlasmidFinder DB
printf "1/13: download incompatibility group signatures from CGE PlasmidFinder DB...\n"
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git
for f in plasmidfinder_db/*.fsa; do
    cat $f >> inc-types-raw.fasta
done
#grep '\S' plasmidfinder_db/gram_positive.fsa > inc-types-raw.fasta
#grep '\S' plasmidfinder_db/enterobacteriaceae.fsa >> inc-types-raw.fasta
cat inc-types-raw.fasta | sed -r "s/^>([a-zA-Z0-9()]+)_([0-9]+)_(.*)_(.*)$/>\1 \4/" > inc-types-mixed.fasta
awk '{if(!/>/){print toupper($0)}else{print $1}}' inc-types-mixed.fasta > inc-types.fasta
rm -rf plasmidfinder_db inc-types-raw.fasta inc-types-mixed.fasta


# download rRNA covariance models from Rfam
printf "\n2/13: download rRNA covariance models from Rfam...\n"
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
printf "\n3/13: download AMR HMMs from NCBIfams...\n"
wget -q -nH ftp://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.LIB
mv NCBIfam-AMRFinder.LIB ncbifam-amr
hmmpress ncbifam-amr
rm ncbifam-amr
wget -q -nH ftp://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.tsv
mv NCBIfam-AMRFinder.tsv ncbifam-amr.tsv


# build HMMs
printf "\n4/13: build HMMs...\n"
wget -q -nH ftp://ftp.ncbi.nlm.nih.gov/genomes/CLUSTERS/PCLA_clusters.txt
touch clusters.lst
for type in replication conjugation; do
    sh $PLATON_HOME/db-scripts/extract-${type}.sh PCLA_clusters.txt
    cut -f1 ${type}.tsv >> clusters.lst
done
grep -f clusters.lst PCLA_proteins.txt | cut -f2 > clusters-proteins.lst
seqtk subseq refseq-bacteria-nrp.trimmed.faa clusters-proteins.lst > clusters-proteins.faa
for type in replication conjugation; do
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


# build conjugation HMMs
printf "\n5/15: build HMMs...\n"
wget -q -nH https://castillo.dicom.unican.es/mobscan_about/MOBfamDB.gz
gunzip MOBfamDB.gz
mv MOBfamDB mobilization
hmmpress mobilization
rm mobilization


# build oriT blastn database
printf "\n6/15: build oriT blastn database...\n"
mkdir mobsuite
cd mobsuite
wget -q -nH https://castillo.dicom.unican.es/mobscan_about/MOBfamDB.gz
unzip mobsuitedb.zip
gunzip orit.fas.gz
mv orit.fas ../orit
cd ..
makeblastdb -dbtype nucl -in orit -title 'OriT'
rm -r orit mobsuite


# download RefSeq reference plasmids
printf "\n7/15: download RefSeq reference plasmids...\n"
wget -q -O refseq-plasmids-raw.tsv ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/plasmids.txt
grep 'Bacteria' refseq-plasmids-raw.tsv | cut -f1,5,6,8,9,10,11,12,15 > refseq-plasmids.tsv
grep 'Bacteria' refseq-plasmids-raw.tsv | cut -f3 refseq-plasmids.tsv > refseq-plasmids-ids.txt
mkdir refseq-plasmids-dir
cd refseq-plasmids-dir
for i in {1..8}; do
    wget -q -nH ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.$i.1.genomic.fna.gz
done
cd ..
gzip -dc refseq-plasmids-dir/plasmid.*.1.genomic.fna.gz | seqtk subseq - refseq-plasmids-ids.txt | seqtk seq -CU > refseq-plasmids
makeblastdb -dbtype nucl -in refseq-plasmids -title 'RefSeq Plasmids'
mv refseq-plasmids refseq-plasmids.fna
rm -r refseq-plasmids-dir refseq-plasmids-raw.tsv refseq-plasmids-ids.txt


# download RefSeq chromosomes
printf "\n8/15: download RefSeq chromosomes...\n"
mkdir nf-tmp
nextflow run $PLATON_HOME/db-scripts/download-chromosomes.nf
rm -rf work .nextflow* nf-tmp


# download NCBI taxonomy
printf "\n9/15: download NCBI taxonomy...\n"
mkdir taxonomy
cd taxonomy
wget -q -nv ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
cd ..


# download UniProt UniRef90 clusters
printf "\n10/15: download UniProt UniRef90 clusters...\n"
wget -q -nv ftp://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref90/uniref90.xml.gz
wget -q -nv ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/uniparc_active.fasta.gz
wget -q -nv ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
python3 $PLATON_HOME/db-scripts/uniref-extract-bacteria.py \
    --taxonomy taxonomy/nodes.dmp \
    --xml uniref90.xml.gz \
    --uniprotkb uniprot_trembl.fasta.gz \
    --uniparc uniparc_active.fasta.gz \
    --fasta uniref90.faa \
    --tsv uniref90.tsv
rm -r uniref90.xml.gz taxonomy/ uniparc_active.fasta.gz uniprot_trembl.fasta.gz


# build complete protein cluster db
printf "\n11/15: build complete protein cluster database...\n"
diamond makedb --in uniref90.faa --db uniref90


# count protein hits
printf "\n12/15: count protein hits...\n"
nextflow run $PLATON_HOME/db-scripts/count-protein-hits.nf \
    --plasmids ./refseq-plasmids.fna \
    --chromosomes ./refseq-chromosomes.fna \
    --pcDb ./uniref90.dmnd
rm -rf work .nextflow*


# calculate RDS, extract found MPS and rebuild database
printf "\n13/15: calculate RDS, extract found MPS and rebuild database...\n"
python3 $PLATON_HOME/db-scripts/setup-mps-db.py \
    --plasmids refseq-plasmids.fna \
    --chromosomes refseq-chromosomes.fna \
    --counts protein-hit-counts.tsv \
    --cluster-seqs uniref90.faa \
    --cluster-info uniref90.tsv
diamond makedb --in mps.faa --db mps
rm uniref90.*


# create artificial contigs
printf "\n14/15: create artificial contigs...\n"
export NXF_OPTS="-Xms256G -Xmx512G"
nextflow run $PLATON_HOME/db-scripts/generate-artificial-contigs.nf \
    --plasmids refseq-plasmids.fna \
    --chromosomes refseq-chromosomes.fna
rm -rf work .nextflow* refseq-plasmids.fna refseq-chromosomes.fna


# compute RDS thresholds
printf "\n15/15: compute RDS thresholds...\n"
export NXF_OPTS="-Xms32G -Xmx256G"
nextflow run $PLATON_HOME/db-scripts/compute-rds-thresholds.nf \
    --contigs artificial-contigs.fna \
    --mps mps.tsv \
    --pcDb mps.dmnd
rm -rf work .nextflow* artificial-contigs.fna

cd ..
