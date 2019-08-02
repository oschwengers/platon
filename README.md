[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/platon/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/cb-platon.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/platon.svg)
![PyPI](https://img.shields.io/pypi/v/cb-platon.svg)
![PyPI - Status](https://img.shields.io/pypi/status/cb-platon.svg)
![Conda](https://img.shields.io/conda/v/bioconda/platon.svg)
![Conda](https://img.shields.io/conda/pn/bioconda/platon.svg)

# Platon: Plasmid contig detection and characterization for short read draft assemblies.

## Contents
-   [Description](#description)
-   [Input/Output](#inputoutput)
-   [Installation](#installation)
-   [Usage](#usage)
-   [Examples](#examples)
-   [Database](#database)
-   [Dependencies](#dependencies)
-   [Citation](#citation)

## Description
Platon detects plasmid contigs from bacterial WGS short read assemblies.
Therefore, Platon computes replicon distribution scores (**RDS**) of marker proteins
sequences (**MPS**) per contig based on pre-computed protein distribution statistics
and tests them against specific thresholds. Contigs whose mean RDS does not reach
the defined thresholds are comprehensively characterized and finally classified
by heuristic filters.

Platon conducts three analysis steps. First, it predicts and searches coding
sequences against a custom and pre-computed database comprising MPS and RDS.
These scores express the measured bias in plasmid/chromosome distributions
based on complete NCBI RefSeq genomes and plasmids.
Platon then calculates the mean RDS for each contig and either classifies them
as chromosome if the RDS is below a sensitivity cutoff (95% sensitivity) or as
plasmid if the RDS is above a specificity cutoff (99.99% specificity).
These thresholds have been set based on Monte Carlo simulations of artifical
subsequences created from complete RefSeq chromosome and plasmid sequences. In a second
step contigs passing the sensitivity filter get comprehensivley characterized.
Hereby, Platon tries to circularize the contig sequences, searches for rRNA,
replication, mobilization and conjugation genes as well as incompatibility group
DNA probes and finally performs a BLAST+ search against the NCBI plasmid database.
In a third step, Platon finally classifies all remaining contigs based on an heuristic
approach, i.e. following a set of heuristic filters.

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
Platon can be installed/used in 2 different ways.

In all cases, the custom database must be downloaded which we provide for download:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3349652.svg)](https://doi.org/10.5281/zenodo.3349652)


### GitHub
1.  clone the the repository
2.  download & extract the database

Example:
```
$ git clone git@github.com:oschwengers/platon.git
$ wget https://zenodo.org/record/3349652/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ platon/bin/platon --db ./db genome.fasta
```

Info: Just move the extracted database directory into the platon directory.
Platon will automatically recognise it and thus, the database path doesn't need
to be specified:
```
$ git clone git@github.com:oschwengers/platon.git
$ wget https://zenodo.org/record/3349652/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ mv db/ platon
$ platon/bin/platon genome.fasta
```

### Conda
1.  install Platon via Conda
2.  download & extract the database

Example:
```
$ conda install -c conda-forge -c bioconda -c defaults platon
$ wget https://zenodo.org/record/3349652/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ platon --db ./db genome.fasta
```

### Pip
1.  install Platon per pip
2.  download and extract the database
3.  install 3rd party binaries

Platon/database (1./2.):
```
$ pip3 install cb-platon
$ wget https://zenodo.org/record/3349652/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
$ platon --db ./db genome.fasta
```

3rd party dependencies on Ubuntu (3.):
```
$ sudo apt install ncbi-blast+ prodigal infernal hmmer mummer
$ wget http://www.bi.cs.titech.ac.jp/ghostz/releases/ghostz-1.0.2.tar.gz
$ tar -xzf ghostz-1.0.2.tar.gz
$ cd ghostz-1.0.2/
$ make
$ sudo cp ghostz /usr/bin/
```
If there are any issues compiling ghostz, please make sure you have everything
correctly setup, e.g. `$ sudo apt install build-essential`.

## Usage
Usage:
```
usage: platon [-h] [--threads THREADS] [--verbose] [--output OUTPUT]
              [--version]
              <genome>

Plasmid contig classification and characterization

positional arguments:
  <genome>              draft genome in fasta format

optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS, -t THREADS
                        number of threads to use (default = number of
                        available CPUs)
  --verbose, -v         print verbose information
  --output OUTPUT, -o OUTPUT
                        output directory (default = current working directory)
  --version             show program's version number and exit
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
PlasmidFinder db as well as custom HMM models. This database based on
RefSeq release 95 can be downloaded here:
(zipped 1.8 Gb, unzipped 2.6 Gb)
-   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3349652.svg)](https://doi.org/10.5281/zenodo.3349652)
-   [https://zenodo.org/record/3349652/files/db.tar.gz](https://zenodo.org/record/3349652/files/db.tar.gz)

## Dependencies
Platon was developed and tested in Python 3.5 and depends on BioPython (>=1.71).

Additionally, it depends on the following 3rd party executables:
-   Prodigal (2.6.3) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648> <https://github.com/hyattpd/Prodigal>
-   Ghostz (1.0.2) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4393512> <http://www.bi.cs.titech.ac.jp/ghostz>
-   Blast+ (2.7.1) <https://www.ncbi.nlm.nih.gov/pubmed/2231712> <https://blast.ncbi.nlm.nih.gov>
-   MUMmer (4.0.0-beta2) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC395750/> <https://github.com/gmarcais/mummer>
-   HMMER (3.2.1) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3695513/> <http://hmmer.org/>
-   INFERNAL (1.1.2) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3810854> <http://eddylab.org/infernal>

## Citation
A manuscript is in preparation... stay tuned!
To temporarily cite our work, please transitionally refer to:
> Schwengers O., Barth P., Falgenhauer L., Hain T., Chakraborty T., Goesmann A. (2019) Platon: Plasmid contig classification and characterization for short read draft assemblies. GitHub https://github.com/oschwengers/platon

As Platon takes advantage of PlasmidFinder's incompatibility database, please also cite:
> Carattoli A., Zankari E., Garcia-Fernandez A., Voldby Larsen M., Lund O., Villa L., Aarestrup F.M., Hasman H. (2014) PlasmidFinder and pMLST: in silico detection and typing of plasmids. Antimicrobial Agents and Chemotherapy, https://doi.org/10.1128/AAC.02412-14
