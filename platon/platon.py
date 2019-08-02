
import argparse
import functools as ft
import json
import multiprocessing as mp
import os
import re
import sys
import shutil
import subprocess as sp
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from Bio import SeqIO

import platon
import platon.functions as pf
import platon.constants as pc


def main():
    # parse arguments
    parser = argparse.ArgumentParser(
        prog='platon',
        description='Plasmid contig classification and characterization'
    )
    parser.add_argument('genome', metavar='<genome>', help='draft genome in fasta format')
    parser.add_argument('--db', '-d', action='store', help='database path (default = <platon_path>/db)')
    parser.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='number of threads to use (default = number of available CPUs)')
    parser.add_argument('--verbose', '-v', action='store_true', help='print verbose information')
    parser.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize & classify all contigs')
    parser.add_argument('--output', '-o', help='output directory (default = current working directory)')
    parser.add_argument('--version', '-V', action='version', version='%(prog)s '+platon.__version__)
    args = parser.parse_args()

    # check parameters and test/setup runtime configuration
    config = pf.setup_configuration()

    if(args.db):
        config['db'] = os.path.abspath(args.db)
    pf.test_database(config)

    if('bundled-binaries' not in config):
        pf.test_binaries()

    genome_path = os.path.abspath(args.genome)
    if(not os.access(genome_path, os.R_OK)):
        sys.exit('ERROR: genome file not readable!')
    if(os.stat(genome_path).st_size == 0):
        sys.exit('ERROR: genome file is empty!')

    output_path = os.path.abspath(args.output) if args.output else os.getcwd()
    if(args.verbose):
        print('Options, parameters and arguments:')
        print('\tdb path: ' + config['db'])
        print('\tuse bundled binaries: ' + str(config['bundled-binaries']))
        print('\tgenome path: ' + genome_path)
        print('\toutput path: ' + output_path)
        print('\tcharacterize: ' + str(args.characterize))
        print('\ttmp path: ' + config['tmp'])
        print('\t# threads: ' + str(args.threads))

    # parse draft genome
    if(args.verbose):
        print('parse draft genome...')
    contigs = {}
    raw_contigs = []
    try:
        for record in SeqIO.parse(genome_path, 'fasta'):
            contig = {
                'id': record.id,
                'length': len(record.seq),
                'sequence': str(record.seq),
                'orfs': {},
                'inc_types': [],
                'amr_hits': [],
                'mobilization_hits': [],
                'replication_hits': [],
                'conjugation_hits': [],
                'rrnas': [],
                'plasmid_hits': []
            }
            raw_contigs.append(contig)

            # read coverage from contig names if they were assembled with SPAdes
            match = re.fullmatch(pc.SPADES_CONTIG_PATTERN, record.id)
            contig['coverage'] = 0 if (match is None) else float(match.group(1))

            # only include contigs with reasonable lengths
            if(args.characterize
                    or contig['length'] >= pc.MIN_CONTIG_LENGTH
                    and contig['length'] < pc.MAX_CONTIG_LENGTH):
                contigs[record.id] = contig
            elif (args.verbose):
                if(contig['length'] < pc.MIN_CONTIG_LENGTH):
                    print('\texclude contig \'%s\', too short (%d)' % (contig['id'], contig['length']))
                elif(contig['length'] >= pc.MAX_CONTIG_LENGTH):
                    print('\texclude contig \'%s\', too long (%d)' % (contig['id'], contig['length']))
    except:
        sys.exit('ERROR: wrong genome file format!')

    if(len(raw_contigs) == 0):
        sys.exit('Error: input file contains no valid contigs.')

    if(len(contigs) == 0):
        print(pc.HEADER)
        if(args.verbose):
            print('No potential plasmid contigs found. Please, check contig lengths. Maybe you passed a finished or pseudo genome?')
        sys.exit(0)

    # predict ORFs
    if(args.verbose):
        print('predict ORFs...')
    proteins_path = pf.predict_orfs(config, contigs, genome_path)
    if(args.verbose):
        no_orfs = ft.reduce(lambda x, y: x+y, map(lambda k: len(contigs[k]['orfs']), contigs))
        print('\tfound %d ORFs' % no_orfs)

    # exclude contigs without ORFs
    if(not args.characterize):
        tmp_contigs = {}
        for contig in filter(lambda k: len(k['orfs']) > 0, contigs.values()):
            tmp_contigs[contig['id']] = contig
        if(args.verbose):
            no_removed_contigs = len(contigs) - len(tmp_contigs)
            print('\tremoved %d contig(s), no ORFs found' % (no_removed_contigs))
        contigs = tmp_contigs

    # find marker genes
    if(args.verbose):
        print('search marker proteins...')
    tmp_output_path = config['tmp'] + '/ghostz.tsv'
    sp.check_call(
        [
            'ghostz',
            'aln',
            '-d', config['db'] + '/refseq-bacteria-nrpc-reps',
            '-b', '1',  # max 1 result per query
            '-a', str(args.threads),  # threads
            '-i', proteins_path,
            '-o', tmp_output_path
        ],
        cwd=config['tmp'],
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT
    )

    # parse ghostz output
    no_proteins_identified = 0
    with open(tmp_output_path, 'r') as fh:
        for line in fh:
            cols = line.split('\t')
            locus = cols[0].split(' ')[0].rpartition('_')
            if((float(cols[2]) >= pc.MIN_PROTEIN_IDENTITY) and (locus[0] in contigs)):
                contig = contigs[locus[0]]
                orf = contig['orfs'][locus[2]]
                orf['protein_id'] = cols[1]
                no_proteins_identified += 1
    if(args.verbose):
        print('\tfound %d marker proteins' % (no_proteins_identified))

    # parse protein score file
    if(args.verbose):
        print('score contigs...')
    marker_proteins = {}
    with open(config['db'] + '/rds.tsv', 'r') as fh:
        for line in fh:
            cols = line.split('\t')
            marker_proteins[cols[0]] = {
                'id': cols[0],
                'product': cols[1],
                'length': cols[2],
                'score': float(cols[3]),
            }

    # calculate protein score per contig
    for contig in contigs.values():
        score_sum = 0.0
        for orf in contig['orfs'].values():
            if(('protein_id' in orf) and (orf['protein_id'] in marker_proteins)):
                marker_protein = marker_proteins[orf['protein_id']]
                score = marker_protein['score']
                orf['score'] = score
                orf['product'] = marker_protein['product']
                score_sum += score
            else:
                score_sum += pc.PROTEIN_SCORE_PENALTY
        contig['protein_score'] = score_sum / len(contig['orfs']) if len(contig['orfs']) > 0 else 0

    # filter contigs based on conservative protein score threshold
    # RDS_SENSITIVITY_THRESHOLD and execute per contig analyses in parallel
    scored_contigs = None
    if(args.characterize):
        scored_contigs = contigs
    else:
        if(args.verbose):
            print('prefilter contigs...')
        scored_contigs = {k: v for (k, v) in contigs.items() if v['protein_score'] >= pc.RDS_SENSITIVITY_THRESHOLD}
        if(args.verbose):
            no_excluded_contigs = len(contigs) - len(scored_contigs)
            print('\texcluded %d contigs by min protein score threshold (%2.1f)' % (no_excluded_contigs, pc.RDS_SENSITIVITY_THRESHOLD))
            print('analyze contigs...')

    # extract proteins from potential plasmid contigs for subsequent analyses
    filtered_proteins_path = config['tmp'] + '/proteins-filtered.faa'
    with open(filtered_proteins_path, 'w') as fh:
        for record in SeqIO.parse(proteins_path, 'fasta'):
            orf_name = str(record.id).split()[0]
            contig_id = orf_name.rsplit('_', 1)[0]
            if(contig_id in scored_contigs):
                fh.write('>' + orf_name + '\n')
                fh.write(str(record.seq) + '\n')

    # write contig sequences to fasta files for subsequent parallel analyses
    for id, contig in scored_contigs.items():
        contig_path = config['tmp'] + '/' + contig['id'] + '.fasta'
        with open(contig_path, 'w') as fh:
            fh.write('>' + contig['id'] + '\n')
            fh.write(contig['sequence'] + '\n')

    # init thread pool to parallize io bound analyses
    with ThreadPoolExecutor(max_workers=args.threads) as tpe:
        # start search for replication, mobilization and conjugation genes
        for fn in (pf.search_replication_genes, pf.search_mobilization_genes, pf.search_conjugation_genes, pf.search_amr_genes):
            tpe.submit(fn, config, scored_contigs, filtered_proteins_path)
        # analyse contigs
        for id, contig in scored_contigs.items():
            tpe.submit(pf.test_circularity, config, contig)
            tpe.submit(pf.search_inc_type, config, contig)
            tpe.submit(pf.search_rrnas, config, contig)
            tpe.submit(pf.search_reference_plasmids, config, contig)

    # lookup AMR genes
    amr_genes = {}
    with open(config['db'] + '/ncbifam-amr.tsv', 'r') as fh:
        for line in fh:
            cols = line.split('\t')
            amr_genes[cols[0]] = {
                'gene': cols[4],
                'product': cols[8]
            }
    for id, contig in scored_contigs.items():
        for hit in contig['amr_hits']:
            amr_gene = amr_genes[hit['hmm-id']]
            hit['gene'] = amr_gene['gene']
            hit['product'] = amr_gene['product']

    # remove tmp dir
    shutil.rmtree(config['tmp'])

    # filter contigs
    filtered_contigs = None
    if(args.characterize):  # skip protein score based filtering
        filtered_contigs = scored_contigs
    else:
        filtered_contigs = {k: v for (k, v) in scored_contigs.items() if pf.filter_contig(v)}

    # get file prefix
    prefix = Path(genome_path).name
    tmp = prefix.split('.')
    if(len(tmp) > 1):  # remove a potential suffix, remaining inner name periods ('.')
        prefix = '.'.join(tmp[:-1:])

    # print results to tsv file and STDOUT
    print(pc.HEADER)
    tmp_output_path = output_path + '/' + prefix + '.tsv'
    with open(tmp_output_path, 'w') as fh:
        fh.write(pc.HEADER + '\n')
        for id in sorted(filtered_contigs, key=lambda k: -filtered_contigs[k]['length']):
            c = filtered_contigs[id]
            if(len(c['inc_types']) == 0):
                inc_types = '-'
            else:
                inc_types = c['inc_types'][0]['type']
                i = 1
                while(i < len(c['inc_types'])):
                    inc_types += ',' + c['inc_types'][i]['type']
                    i += 1
            line = '%s\t%d\t%4.1f\t%d\t%3.1f\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d' % (c['id'], c['length'], c['coverage'], len(c['orfs']), c['protein_score'], 'yes' if c['is_circular'] else 'no', inc_types, len(c['replication_hits']), len(c['mobilization_hits']), len(c['conjugation_hits']), len(c['amr_hits']), len(c['rrnas']), len(c['plasmid_hits']))
            print(line)
            fh.write(line + '\n')

    # write comprehensive results to JSON file
    tmp_output_path = output_path + '/' + prefix + '.json'
    with open(tmp_output_path, 'w') as fh:
        indent = '\t' if args.verbose else None
        separators = (', ', ': ') if args.verbose else (',', ':')
        json.dump(filtered_contigs, fh, indent=indent, separators=separators)

    # write chromosome contigs to fasta file
    tmp_output_path = output_path + '/' + prefix + '.chromosome.fasta'
    with open(tmp_output_path, 'w') as fh:
        for contig in raw_contigs:
            if(contig['id'] not in filtered_contigs):
                fh.write('>' + contig['id'] + '\n')
                fh.write(contig['sequence'] + '\n')

    # write plasmid contigs to fasta file
    tmp_output_path = output_path + '/' + prefix + '.plasmid.fasta'
    with open(tmp_output_path, 'w') as fh:
        for contig in raw_contigs:
            if(contig['id'] in filtered_contigs):
                fh.write('>' + contig['id'] + '\n')
                fh.write(contig['sequence'] + '\n')


if __name__ == '__main__':
    main()
