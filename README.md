[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/platon/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/cb-platon.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/platon.svg)
![PyPI](https://img.shields.io/pypi/v/cb-platon.svg)
![PyPI - Status](https://img.shields.io/pypi/status/cb-platon.svg)
![Conda](https://img.shields.io/conda/v/bioconda/platon.svg)
![Conda](https://img.shields.io/conda/pn/bioconda/platon.svg)

# Platon: identification and characterization of bacterial plasmid contigs from short-read draft assemblies.

## Contents
- [Description](#description)
- [Input/Output](#inputoutput)
- [Installation](#installation)
  - [BioConda](#bioconda)
  - [GitHub](#github)
  - [Pip](#pip)
- [Usage](#usage)
- [Examples](#examples)
- [Database](#database)
- [Dependencies](#dependencies)
- [Citation](#citation)
- [Issues](#issues)

## Description
**TL;DR**
Platon detects plasmid contigs within bacterial draft genomes from WGS short-read assemblies.
Therefore, Platon analyses the natural distribution biases of certain protein coding genes between
chromosomes and plasmids. This analysis is complemented by comprehensive contig characterizations
upon which the software applies several heuristics.

Platon conducts three analysis steps:
1. It predicts and searches coding sequences against a custom and pre-computed
database comprising marker protein sequences (**MPS**) and related replicon
distribution scores (**RDS**). These scores express the empirically measured
frequency biases of protein sequence distributions between plasmids and chromosomes
pre-computed on complete NCBI RefSeq replicons. Platon calculates the mean RDS for
each contig and either classifies them as chromosome if the RDS is below a
sensitivity cutoff determined to 95% sensitivity or as plasmid if the RDS is
above a specificity cutoff determined to 99.9% specificity.
Exact values for these thresholds have been computed based on Monte Carlo simulations of
artifical replicon fragments created from complete RefSeq chromosome and plasmid sequences.
2. Contigs passing the sensitivity filter get comprehensivley characterized.
Hereby, Platon tries to circularize the contig sequences, searches for rRNA,
replication, mobilization and conjugation genes, oriT sequences, incompatibility group
DNA probes and finally performs a BLAST+ search against the NCBI plasmid database.
3. Finally, to increase the overall sensitivity, Platon classifies all remaining contigs based on the gathered information
by several heuristics.

| ![Replicon distribution and alignment hit frequencies of marker protein sequences](rds-ratio-counts.web.png?raw=true) |
| -- |
| *Fig: Replicon distribution and alignment hit frequencies of marker protein sequences. Shown are summed plasmid and chromosome alignment hit frequencies per marker protein sequences plotted against plasmid/chromosome hit count ratios scaled to [-1, 1]; Hue: normalized replicon distribution score values (min=-100, max=100), hit count outliers below 10-4 and above 1 are discarded for the sake of readability.* |

## Input/Output

### Input
Platon accepts draft assemblies in fasta format. If contigs have been assembled with
SPAdes, Platon is able to extract the coverage information from the contig names.

### Output
For each contig classified as plasmid sequence the following columns are printed
to `STDOUT` as tab separated values:
-   Contig ID
-   Length
-   Coverage
-   \# ORFs
-   Protein Score
-   Circularity
-   Incompatibility Type(s)
-   \# Replication Genes
-   \# Mobilization Genes
-   \# OriT Sequences
-   \# Conjugation Genes
-   \# rRNA Genes
-   \# Plasmid Database Hits

In addition, Platon writes the following files into the output directory:
-   `<prefix>`.plasmid.fasta: contigs classified as plasmids or plasmodal origin
-   `<prefix>`.chromosome.fasta: contigs classified as chromosomal origin
-   `<prefix>`.tsv: dense information as printed to STDOUT (see above)
-   `<prefix>`.json: comprehensive results and information on each single plasmid contig.
All files are prefixed (`<prefix>`) as the input genome fasta file.

## Installation
Platon can be installed in 3 different ways, though we advise to use Conda/BioConda.

In all cases, the custom database must be downloaded which we provide for download:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3349651.svg)](https://doi.org/10.5281/zenodo.3349651)

### BioConda
1.  install Platon via [Conda](https://conda.io/docs/install/quick.html) and the [Bioconda](https://bioconda.github.io/) channel
2.  download & extract the database

Example:
```
$ conda install -c conda-forge -c bioconda -c defaults platon
$ wget https://zenodo.org/record/3358926/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ platon --db ./db genome.fasta
```

### GitHub
1.  clone the the repository
2.  download & extract the database

Example:
```
$ git clone git@github.com:oschwengers/platon.git
$ wget https://zenodo.org/record/3358926/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ platon/bin/platon --db ./db genome.fasta
```

Info: Just move the extracted database directory into the platon directory.
Platon will automatically recognise it and thus, the database path doesn't need
to be specified:
```
$ git clone git@github.com:oschwengers/platon.git
$ wget https://zenodo.org/record/3358926/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ mv db/ platon
$ platon/bin/platon genome.fasta
```

### Pip
1.  install Platon per pip
2.  download and extract the database
3.  install 3rd party binaries

Platon/database (1./2.):
```
$ pip3 install cb-platon
$ wget https://zenodo.org/record/3358926/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ platon --db ./db genome.fasta
```

3rd party dependencies on Ubuntu (3.):
```
$ sudo apt install ncbi-blast+ prodigal diamond-aligner infernal hmmer mummer
```

## Usage
Usage:
```
usage: platon [-h] [--db DB] [--threads THREADS] [--verbose] [--characterize]
              [--output OUTPUT] [--version]
              <genome>

Identification and characterization of bacterial plasmid contigs from short-read draft assemblies.

positional arguments:
  <genome>              draft genome in fasta format

optional arguments:
  -h, --help            show this help message and exit
  --db DB, -d DB        database path (default = <platon_path>/db)
  --threads THREADS, -t THREADS
                        number of threads to use (default = number of
                        available CPUs)
  --verbose, -v         print verbose information
  --characterize, -c    deactivate filters; characterize all contigs
  --output OUTPUT, -o OUTPUT
                        output directory (default = current working directory)
  --version, -V         show program's version number and exit
```

## Examples
Simple:
```
$ platon genome.fasta
```

Expert: writing results to `results` directory with verbose output using 8 threads:
```
$ platon -db ~/db --output results/ --verbose --threads 8 genome.fasta
```

## Database
Platon depends on a custom database based on MPS, RDS, RefSeq Plasmid database,
PlasmidFinder db as well as manually curated MOB HMM models from MOBscan,
custom conjugation and replication HMM models and oriT sequences from MOB-suite.
This database based on UniProt UniRef90 release 2020_01 can be downloaded here:
(zipped 1.4 Gb, unzipped 2.4 Gb)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3349651.svg)](https://doi.org/10.5281/zenodo.3349651)
-   [https://zenodo.org/record/3358926/files/db.tar.gz](https://zenodo.org/record/3358926/files/db.tar.gz)

## Dependencies
Platon was developed and tested in Python 3.5 and depends on BioPython (>=1.71).

Additionally, it depends on the following 3rd party executables:
-   Prodigal (2.6.3) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648> <https://github.com/hyattpd/Prodigal>
-   Diamond (0.9.31) <https://pubmed.ncbi.nlm.nih.gov/25402007> <http://www.diamondsearch.org>
-   Blast+ (2.7.1) <https://www.ncbi.nlm.nih.gov/pubmed/2231712> <https://blast.ncbi.nlm.nih.gov>
-   MUMmer (4.0.0-beta2) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC395750/> <https://github.com/gmarcais/mummer>
-   HMMER (3.2.1) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3695513/> <http://hmmer.org/>
-   INFERNAL (1.1.2) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3810854> <http://eddylab.org/infernal>

## Citation
A manuscript is submitted.
To temporarily cite our work, please transitionally refer to:
> Schwengers O., Barth P., Falgenhauer L., Hain T., Chakraborty T., Goesmann A. (2020) Platon: identification and characterization of bacterial plasmid contigs in short-read draft assemblies exploiting protein-sequence-based replicon distribution scores. GitHub https://github.com/oschwengers/platon

As Platon takes advantage of PlasmidFinder's incompatibility database, please also cite:
> Carattoli A., Zankari E., Garcia-Fernandez A., Voldby Larsen M., Lund O., Villa L., Aarestrup F.M., Hasman H. (2014) PlasmidFinder and pMLST: in silico detection and typing of plasmids. Antimicrobial Agents and Chemotherapy, https://doi.org/10.1128/AAC.02412-14

As Platon takes advantage of MOBscan's MOB HMM profiles, please also cite:
> Garcillán-Barcia M. P., Redondo-Salvo S., Vielva L., de la Cruz F. (2020) MOBscan: Automated Annotation of MOB Relaxases. Methods in Molecular Biology, https://doi.org/10.1007/978-1-4939-9877-7_21

As Platon takes advantage of MOB-suite's oriT sequences, please also cite:
> Robertson J., Nash J. H. E. (2018) MOB-suite: Software Tools for Clustering, Reconstruction and Typing of Plasmids From Draft Assemblies. Microbial Genomics, https://doi.org/10.1099/mgen.0.000206


## Issues
If you run into any issues with Platon, we'd be happy to hear about it!
Please, start the pipeline with `-v` (verbose) and do not hesitate
to file an issue including as much of the following as possible:
- a detailed description of the issue
- the platon cmd line output
- the `<prefix>.json` file if possible
- A reproducible example of the issue with a small dataset that you can share
(helps us identify whether the issue is specific to a particular computer, operating system, and/or dataset).
