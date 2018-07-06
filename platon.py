#!/usr/bin/env python3

import argparse
import functools as ft
import json
import multiprocessing as mp
import os
import re
import sys
import shutil
import subprocess as sp
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from Bio import SeqIO

PLATON_HOME = os.getenv( 'PLATON_HOME', None )
if( PLATON_HOME is None ): sys.exit( 'ERROR: PLATON_HOME not set!' )

from platon_functions import *


# parse arguments
parser = argparse.ArgumentParser( prog='platon',
    description='PLAsmid conTig extractiOn' )
parser.add_argument( 'genome', metavar='<genome>', help='Draft genome in fasta format' )
parser.add_argument( '--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)' )
parser.add_argument( '--verbose', '-v', action='store_true', help='Print verbose information' )
parser.add_argument( '--output', '-o', help='Output directory' )
parser.add_argument( '--version', action='version', version='%(prog)s 0.1' )
args = parser.parse_args()


# global constants
MIN_PROTEIN_IDENTITY = 80.0
MIN_PROTEIN_SCORE_THRESHOLD = -18.4  # counts for sensitivity => 99 %
PROTEIN_SCORE_PENALTY = -1.0
SPADES_CONTIG_PATTERN = re.compile( 'NODE_\d+_length_\d+_cov_(\d+\.\d+)' )


# check parameters & environment variables
if( args.verbose ): print( 'Options, parameters and arguments:' )


PLATON_HOME = os.path.abspath( PLATON_HOME )
if( not os.access( PLATON_HOME, os.R_OK & os.X_OK ) ): sys.exit( 'ERROR: PLATON_HOME ('+PLATON_HOME+') not readable/accessible!' )
if( args.verbose ): print( '\tPLATON_HOME: ' + PLATON_HOME )

genomePath = os.path.abspath( args.genome )
if( not os.access( genomePath, os.R_OK ) ): sys.exit( 'ERROR: genome file not readable!' )
if( os.stat( genomePath ).st_size == 0 ): sys.exit( 'ERROR: genome file is empty!' )
if( args.verbose ): print( '\tgenome path: ' + genomePath )

outputPath = os.path.abspath(args.output) if args.output else os.getcwd()
if( args.verbose ): print( '\toutput path: ' + outputPath )

tmpPath = tempfile.mkdtemp()
if( args.verbose ): print( '\ttmp path: ' + tmpPath )

if( args.verbose ): print( '\t# threads: ' + str(args.threads) )

fhFNULL = open( os.devnull, 'w' )


# parse draft genome
if( args.verbose ): print( 'parse draft genome...' )
contigs = {}
rawContigs = []
filteredDraftGenomePath = tmpPath + '/genome.fasta'
with open( filteredDraftGenomePath, 'w' ) as fh:
    try:
        for record in SeqIO.parse( genomePath, 'fasta' ):
            contig = {
                'id': record.id,
                'length': len(record.seq),
                'sequence': str(record.seq),
                'orfs': {},
                'inc_types': [],
                'mobilization_hits': [],
                'replication_hits': [],
                'conjugation_hits': [],
                'rrnas': [],
                'plasmid_hits': []
            }
            rawContigs.append(contig)

            # read coverage from contig names if they were assembled with SPAdes
            match = re.fullmatch( SPADES_CONTIG_PATTERN, record.id )
            contig['coverage'] = 0 if (match is None) else float( match.group(1) )

            # only include contigs with reasonable lengths
            if( contig['length'] >= 1000  and  contig['length'] < 500000 ):
                contigs[ record.id ] = contig
                fh.write( '>' + contig['id'] + '\n' )
                fh.write( contig['sequence'] + '\n' )
            else:
                if( args.verbose ):
                    if( contig['length'] < 1000 ):
                        print( '\texclude contig \'%s\', too short (%d)' % (contig['id'], contig['length']) )
                    elif( contig['length'] >= 500000 ):
                        print( '\texclude contig \'%s\', too long (%d)' % (contig['id'], contig['length']) )
    except:
        sys.exit( 'ERROR: wrong genome file format!' )


if( len(rawContigs) == 0 ):
    sys.exit(0)


# predict ORFs
if( args.verbose ): print( 'predict ORFs...' )
proteinsPath = predict_orfs( tmpPath, contigs, filteredDraftGenomePath )
if( args.verbose ):
    noOrfs = ft.reduce( lambda x,y: x+y , map(lambda k: len(contigs[k]['orfs']), contigs))
    print( '\tfound %d ORFs' % noOrfs )


# exclude contigs without ORFs
tmpContigs = {}
for contig in filter( lambda k: len(k['orfs']) > 0, contigs.values() ):
    tmpContigs[ contig['id'] ] = contig
if( args.verbose ):
    noRemovedContigs = len(contigs) - len(tmpContigs)
    print( '\tremoved %d contig(s), no ORFs found' % (noRemovedContigs) )
contigs = tmpContigs


# find marker genes
if( args.verbose ): print( 'search marker proteins...' )
outPath = tmpPath + '/ghostz.tsv'
sp.check_call( [ PLATON_HOME + '/share/ghostz',
        'aln',
        '-d', PLATON_HOME + '/db/refseq-bacteria-nrpc-reps',
        '-b', '1', # max 1 result per query
        '-a', str(args.threads), # threads
        '-i', proteinsPath,
        '-o', outPath
    ],
    cwd = tmpPath,
    stdout = fhFNULL,
    stderr = sp.STDOUT
)


# parse ghostz output
noProteinsIdentified = 0
with open( outPath, 'r' ) as fh:
        for line in fh:
            cols = line.split( '\t' )
            locus = cols[0].split( ' ' )[0].rpartition( '_' )
            if( float(cols[2]) >= MIN_PROTEIN_IDENTITY ):
                contig = contigs[ locus[0] ]
                orf = contig['orfs'][ locus[2] ]
                orf['protein_id'] = cols[1]
                noProteinsIdentified += 1
if( args.verbose ): print( '\tfound %d marker proteins' % (noProteinsIdentified) )


# parse protein score file
if( args.verbose ): print( 'score and prefilter contigs...' )
markerProteins = {}
with open( PLATON_HOME + '/db/protein-scores.tsv', 'r' ) as fh:
    for line in fh:
        cols = line.split( '\t' )
        markerProteins[ cols[0] ] = {
            'id': cols[0],
            'product': cols[1],
            'length': cols[2],
            'score': float(cols[3]),
        }


# calculate protein score per contig
for contig in contigs.values():
    scoreSum = 0.0
    for orf in contig['orfs'].values():
        if( ('protein_id' in orf)  and  (orf['protein_id'] in markerProteins) ):
            markerProtein = markerProteins[ orf['protein_id'] ]
            score = markerProtein['score']
            orf['score'] = score
            orf['product'] = markerProtein['product']
            scoreSum += score
        else:
            scoreSum += PROTEIN_SCORE_PENALTY
    contig['protein_score'] = scoreSum / len(contig['orfs']) if len(contig['orfs']) > 0 else 0


# filter contigs based on conservative protein score threshold
# MIN_PROTEIN_SCORE_THRESHOLD and execute per contig analyses in parallel
scoredContigs = { k:v for (k,v) in contigs.items() if v['protein_score'] >= MIN_PROTEIN_SCORE_THRESHOLD }
if( args.verbose ):
    noExcludedContigs = len(contigs) - len(scoredContigs)
    print( '\texcluded %d contigs by min protein score threshold (%2.1f)' % (noExcludedContigs, MIN_PROTEIN_SCORE_THRESHOLD) )
    print( 'analyze contigs...' )


# extract proteins from potential plasmid contigs for subsequent analyses
filteredProteinsPath = tmpPath + '/proteins-filtered.faa'
with open( filteredProteinsPath, 'w' ) as fh:
    for record in SeqIO.parse( proteinsPath, 'fasta' ):
        orfName = str(record.id).split()[0]
        contigId = orfName.rsplit('_', 1 )[0]
        if( contigId in scoredContigs ):
            fh.write( '>' + orfName + '\n' )
            fh.write( str(record.seq) + '\n' )


# write contig sequences to fasta files for subsequent parallel analyses
for cId, contig in scoredContigs.items():
    contigPath = tmpPath + '/' + contig['id'] + '.fasta'
    with open( contigPath, 'w' ) as fh:
        fh.write( '>' + contig['id'] + '\n' )
        fh.write( contig['sequence'] + '\n' )


# init thread pool to parallize io bound analyses
with ThreadPoolExecutor( max_workers=args.threads ) as tpe:
    # start search for replication, mobilization and conjugation genes
    for fn in (search_replication_genes, search_mobilization_genes, search_conjugation_genes):
        tpe.submit( fn, tmpPath, scoredContigs, filteredProteinsPath )
    # analyse contigs
    for cId, contig in scoredContigs.items():
        tpe.submit( test_circularity, tmpPath, contig )
        tpe.submit( search_inc_type, tmpPath, contig )
        tpe.submit( search_rrnas, tmpPath, contig )
        tpe.submit( search_reference_plasmids, tmpPath, contig )


# remove tmp dir
shutil.rmtree( tmpPath )
fhFNULL.close()


# filter contigs
filteredContigs = {k:v for (k,v) in scoredContigs.items() if filter_contig(v) }


# get file prefix
prefix = Path(genomePath).name
if( '.' in prefix ): prefix = prefix.split('.')[0]


# print results to tsv file and STDOUT
header = 'ID\tLength\tCoverage\t# ORFs\tProtein Score\tCircular\tInc Type(s)\t# Replication\t# Mobilization\t# Conjugation\t# rRNAs\t# Plasmid Hits'
print( header )
outPath = outputPath + '/' + prefix + '.tsv'
with open( outPath, 'w' ) as fh:
    fh.write( header + '\n' )
    for cId in sorted(filteredContigs, key=lambda k: -filteredContigs[k]['length']):
        c = filteredContigs[cId]
        line = '%s\t%d\t%4.1f\t%d\t%3.1f\t%s\t%s\t%d\t%d\t%d\t%d\t%d' % (c['id'], c['length'], c['coverage'], len(c['orfs']), c['protein_score'], 'yes' if c['is_circular'] else 'no', c['inc_types'][0]['type'] if len(c['inc_types'])>0 else '-', len(c['replication_hits']), len(c['mobilization_hits']), len(c['conjugation_hits']), len(c['rrnas']), len(c['plasmid_hits']))
        print( line )
        fh.write( line + '\n' )


# write comprehensive results to JSON file
outPath = outputPath + '/' + prefix + '.json'
with open( outPath, 'w' ) as fh:
    indent = '\t' if args.verbose else None
    separators = ( ', ', ': ' ) if args.verbose else ( ',', ':' )
    json.dump( filteredContigs, fh, indent=indent, separators=separators )


# write chromosome contigs to fasta file
outPath = outputPath + '/' + prefix + '.chromosome.fasta'
with open( outPath, 'w' ) as fh:
    for contig in rawContigs:
        if( not contig['id'] in filteredContigs ):
            fh.write( '>' + contig['id'] + '\n' )
            fh.write( contig['sequence'] + '\n' )


# write plasmid contigs to fasta file
outPath = outputPath + '/' + prefix + '.plasmid.fasta'
with open( outPath, 'w' ) as fh:
    for contig in rawContigs:
        if( contig['id'] in filteredContigs ):
            fh.write( '>' + contig['id'] + '\n' )
            fh.write( contig['sequence'] + '\n' )
