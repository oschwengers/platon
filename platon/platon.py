
import argparse
import functools as ft
import json
import logging
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
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Identification and characterization of bacterial plasmid contigs from short-read draft assemblies.',
        epilog="Citation:\n%s\n\nGitHub:\nhttps://github.com/oschwengers/platon" % pc.CITATION
    )
    parser.add_argument('genome', metavar='<genome>', help='draft genome in fasta format')
    parser.add_argument('--db', '-d', action='store', help='database path (default = <platon_path>/db)')
    parser.add_argument('--mode', '-m', action='store', type=str, choices=['sensitivity', 'accuracy', 'specificity'], default='accuracy', help='applied filter mode: sensitivity: RDS only (>= 95%% sensitivity); specificity: RDS only (>=99.9%% specificity); accuracy: RDS & characterization heuristics (highest accuracy) (default = accuracy)')
    parser.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize all contigs')
    parser.add_argument('--output', '-o', help='output directory (default = current working directory)')
    parser.add_argument('--prefix', '-p', action='store', default='', help='file prefix (default = input file name)')
    parser.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='number of threads to use (default = number of available CPUs)')
    parser.add_argument('--verbose', '-v', action='store_true', help='print verbose information')
    parser.add_argument('--version', '-V', action='version', version='%(prog)s ' + platon.__version__)
    args = parser.parse_args()

    # check input file
    genome_path = Path(args.genome).resolve()
    if(not os.access(str(genome_path), os.R_OK)):
        sys.exit('ERROR: genome file (%s) not readable!' % genome_path)
    if(genome_path.stat().st_size == 0):
        sys.exit('ERROR: genome file (%s) is empty!' % genome_path)

    # check output file
    output_path = Path(args.output) if args.output else Path.cwd()
    if(not output_path.exists()):
        output_path.mkdir(parents=True, exist_ok=True)
    output_path = output_path.resolve()

    # get file prefix
    prefix = args.prefix if args.prefix != '' else genome_path.stem

    # setup logging
    logging.basicConfig(
        filename='%s/%s.log' % (str(output_path), prefix),
        filemode='w',
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO
    )
    log = logging.getLogger('main')
    log.info('version %s', platon.__version__)

    # check parameters and test/setup runtime configuration
    config = pf.setup_configuration(args)

    if(args.db):
        db_path = Path(args.db).resolve()
        config['db'] = db_path
    pf.test_database(config)

    if('bundled-binaries' not in config):
        pf.test_binaries()

    log.info('configuration: db-path=%s', config['db'])
    log.info('configuration: bundled binaries=%s', config['bundled-binaries'])
    log.info('configuration: tmp-path=%s', config['tmp'])
    log.info('parameters: genome=%s', genome_path)
    log.info('parameters: mode=%s', args.mode)
    log.info('parameters: output=%s', output_path)
    log.info('parameters: prefix=%s', prefix)
    log.info('options: characterize=%s', args.characterize)
    log.info('options: threads=%d', args.threads)
    if(args.verbose):
        print('Options, parameters and arguments:')
        print('\tdb path: ' + str(config['db']))
        print('\tuse bundled binaries: ' + str(config['bundled-binaries']))
        print('\tgenome path: ' + str(genome_path))
        print('\toutput path: ' + str(output_path))
        print('\tprefix: ' + prefix)
        print('\tmode: ' + args.mode)
        print('\tcharacterize: ' + str(args.characterize))
        print('\ttmp path: ' + str(config['tmp']))
        print('\t# threads: ' + str(args.threads))

    # parse draft genome
    if(args.verbose):
        print('parse draft genome...')
    contigs = {}
    raw_contigs = []
    try:
        for record in SeqIO.parse(str(genome_path), 'fasta'):
            id = record.id
            length = len(record.seq)
            contig = {
                'id': id,
                'length': length,
                'sequence': str(record.seq),
                'orfs': {},
                'is_circular': False,
                'inc_types': [],
                'amr_hits': [],
                'mobilization_hits': [],
                'orit_hits': [],
                'replication_hits': [],
                'conjugation_hits': [],
                'rrnas': [],
                'plasmid_hits': []
            }
            raw_contigs.append(contig)

            # read coverage from contig names if they were assembled with SPAdes
            match = re.fullmatch(pc.SPADES_CONTIG_PATTERN, id)
            contig['coverage'] = 0 if (match is None) else float(match.group(1))

            # only include contigs with reasonable lengths except of
            # platon runs in characterization mode
            if(args.characterize):
                contigs[id] = contig
            else:
                if(length < pc.MIN_CONTIG_LENGTH):
                    log.info('exclude contig: too short: id=%s, length=%d', id, length)
                    if (args.verbose):
                        print('\texclude contig \'%s\', too short (%d)' % (id, length))
                elif(length >= pc.MAX_CONTIG_LENGTH):
                    log.info('exclude contig: too long: id=%s, length=%d', id, length)
                    if (args.verbose):
                        print('\texclude contig \'%s\', too long (%d)' % (id, length))
                else:
                    contigs[id] = contig
    except:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')

    if(args.verbose):
        print('\tparsed %d raw contigs' % len(raw_contigs))
        print('\texcluded %d contigs by size filter' % (len(raw_contigs) - len(contigs)))
        print('\tanalyze %d contigs' % len(contigs))
    log.info(
        'length contig filter: # input=%d, # discarded=%d, # remaining=%d',
        len(raw_contigs), (len(raw_contigs) - len(contigs)), len(contigs)
    )

    if(len(raw_contigs) == 0):
        log.warning('no valid contigs!')
        sys.exit('Error: input file contains no valid contigs.')

    if(len(contigs) == 0):
        print(pc.HEADER)
        log.warning('no potential plasmid contigs!')
        if(args.verbose):
            print('No potential plasmid contigs found. Please, check contig lengths. Maybe you passed a finished or pseudo genome?')
        sys.exit(0)

    # predict ORFs
    if(args.verbose):
        print('predict ORFs...')
    proteins_path = pf.predict_orfs(config, contigs, genome_path)
    if(proteins_path is None):
        sys.exit('Error: ORF prediction failed!')
    no_orfs = ft.reduce(lambda x, y: x + y, map(lambda k: len(contigs[k]['orfs']), contigs))
    log.info('ORF detection: # ORFs=%d', no_orfs)
    if(args.verbose):
        print('\tfound %d ORFs' % no_orfs)

    # retain / exclude contigs without ORFs
    if(args.characterize or args.mode == 'sensitivity' or args.mode == 'accuracy'):
        log.info('ORF contig filter disabled! # passed contigs=%s', len(contigs))
    else:  # exclude contigs without ORFs in specificity mode
        tmp_contigs = {}
        for contig in filter(lambda k: len(k['orfs']) > 0, contigs.values()):
            tmp_contigs[contig['id']] = contig
        no_removed_contigs = len(contigs) - len(tmp_contigs)
        if(args.verbose):
            print('\tdiscarded %d contig(s) without ORFs' % no_removed_contigs)
        contigs = tmp_contigs
        log.info('ORF contig filter: # discarded=%s, # remaining=%s', no_removed_contigs, len(contigs))

    # find marker genes
    if(args.verbose):
        print('search marker protein sequences (MPS)...')
    tmp_output_path = config['tmp'].joinpath('ghostz.tsv')
    cmd = [
        'diamond',
        'blastp',
        '--db', str(config['db'].joinpath('mps.dmnd')),
        '--query', str(proteins_path),
        '--out', str(tmp_output_path),
        '--max-target-seqs', '1',  # max 1 result per query
        '--id', '90',  # min alignment identity 90%
        '--query-cover', '80',  # min query cov 80%
        '--subject-cover', '80',  # min subjetc cov 80%
        '--threads', str(args.threads),  # threads
        '--tmpdir', str(config['tmp'])
    ]
    proc = sp.run(
        cmd,
        cwd=str(config['tmp']),
        env=config['env'],
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.error(
            'diamond execution failed! diamond-error-code=%d',
            proc.returncode
        )
        log.debug(
            'diamond execution: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        sys.exit('Marker protein search failed!')

    # parse diamond output
    proteins_identified = 0
    with tmp_output_path.open() as fh:
        for line in fh:
            cols = line.split('\t')
            locus = cols[0].rpartition('_')
            contig_id = locus[0]
            orf_id = locus[2]
            if((float(cols[2]) >= pc.MIN_PROTEIN_IDENTITY) and (contig_id in contigs)):
                contig = contigs[contig_id]
                orf = contig['orfs'][orf_id]
                orf['protein_id'] = cols[1]
                proteins_identified += 1
    log.info('MPS detection: # MPS=%d', proteins_identified)
    if(args.verbose):
        print('\tfound %d MPS' % proteins_identified)

    # parse protein score file
    if(args.verbose):
        print('compute replicon distribution scores (RDS)...')
    marker_proteins = {}
    with config['db'].joinpath('mps.tsv').open() as fh:
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
        contig['protein_score'] = score_sum / len(contig['orfs']) if len(contig['orfs']) > 0 else 0
        log.info(
            'contig RDS: contig=%s, RDS=%f, score-sum=%f, #ORFs=%d',
            contig['id'], contig['protein_score'], score_sum, len(contig['orfs'])
        )

    # filter contigs based on conservative protein score threshold
    # RDS_SENSITIVITY_THRESHOLD and execute per contig analyses in parallel
    scored_contigs = None
    if(args.characterize):
        scored_contigs = contigs
    else:
        if(args.verbose):
            print('apply RDS sensitivity threshold (SNT=%2.1f) filter...' % pc.RDS_SENSITIVITY_THRESHOLD)
        scored_contigs = {k: v for (k, v) in contigs.items() if v['protein_score'] >= pc.RDS_SENSITIVITY_THRESHOLD}
        no_excluded_contigs = len(contigs) - len(scored_contigs)
        log.info('RDS SNT filter: # discarded contigs=%d, # remaining contigs=%d', no_excluded_contigs, len(scored_contigs))
        if(args.verbose):
            print('\texcluded %d contigs by SNT filter' % no_excluded_contigs)
            print('characterize contigs...')

    # extract proteins from potential plasmid contigs for subsequent analyses
    filtered_proteins_path = config['tmp'].joinpath('proteins-filtered.faa')
    with filtered_proteins_path.open(mode='w') as fh:
        for record in SeqIO.parse(str(proteins_path), 'fasta'):
            orf_name = str(record.id).split()[0]
            contig_id = orf_name.rsplit('_', 1)[0]
            if(contig_id in scored_contigs):
                fh.write('>' + orf_name + '\n')
                fh.write(str(record.seq) + '\n')

    # write contig sequences to fasta files for subsequent parallel analyses
    for id, contig in scored_contigs.items():
        contig_path = config['tmp'].joinpath(contig['id'] + '.fasta')
        with contig_path.open(mode='w') as fh:
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
            tpe.submit(pf.search_orit_sequences, config, contig)
            tpe.submit(pf.search_reference_plasmids, config, contig)

    # lookup AMR genes
    amr_genes = {}
    with config['db'].joinpath('ncbifam-amr.tsv').open() as fh:
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
    shutil.rmtree(str(config['tmp']))
    log.debug('removed tmp dir: %s', config['tmp'])

    # filter contigs
    filtered_contigs = None
    if(args.characterize):  # skip protein score based filtering
        filtered_contigs = scored_contigs
    elif(args.mode == 'sensitivity'):  # skip protein score based filtering but apply rRNA filter
        filtered_contigs = {k: v for (k, v) in scored_contigs.items() if len(v['rrnas']) == 0}
    elif(args.mode == 'specificity'):
        filtered_contigs = {k: v for (k, v) in scored_contigs.items() if v['protein_score'] >= pc.RDS_SPECIFICITY_THRESHOLD}
    else:
        filtered_contigs = {k: v for (k, v) in scored_contigs.items() if pf.filter_contig(v)}

    # print results to tsv file and STDOUT
    print(pc.HEADER)
    tmp_output_path = output_path.joinpath(prefix + '.tsv')
    log.debug('output: tsv=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        fh.write(pc.HEADER + '\n')
        for id in sorted(filtered_contigs, key=lambda k: -filtered_contigs[k]['length']):
            c = filtered_contigs[id]
            cov = 'NA' if c['coverage'] == 0 else "%4.1f" % c['coverage']
            line = '%s\t%d\t%s\t%d\t%3.1f\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % (c['id'], c['length'], cov, len(c['orfs']), c['protein_score'], 'yes' if c['is_circular'] else 'no', len(c['inc_types']), len(c['replication_hits']), len(c['mobilization_hits']), len(c['orit_hits']), len(c['conjugation_hits']), len(c['amr_hits']), len(c['rrnas']), len(c['plasmid_hits']))
            print(line)
            log.info(
                'plasmid: id=%s, len=%d, cov=%s, ORFs=%d, RDS=%f, circ=%s, incs=%s, # rep=%d, # mob=%d, #oriT=%d, # con=%d, # AMRs=%d, # rRNAs=%d, # refs=%d',
                c['id'], c['length'], cov, len(c['orfs']), c['protein_score'], c['is_circular'], len(c['inc_types']), len(c['replication_hits']), len(c['mobilization_hits']), len(c['orit_hits']), len(c['conjugation_hits']), len(c['amr_hits']), len(c['rrnas']), len(c['plasmid_hits'])
            )
            fh.write(line + '\n')

    # write comprehensive results to JSON file
    tmp_output_path = output_path.joinpath(prefix + '.json')
    log.debug('output: json=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        indent = '\t' if args.verbose else None
        separators = (', ', ': ') if args.verbose else (',', ':')
        json.dump(filtered_contigs, fh, indent=indent, separators=separators)

    # write chromosome contigs to fasta file
    tmp_output_path = output_path.joinpath(prefix + '.chromosome.fasta')
    log.debug('output: chromosomes=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        for contig in raw_contigs:
            if(contig['id'] not in filtered_contigs):
                fh.write('>' + contig['id'] + '\n')
                fh.write(contig['sequence'] + '\n')

    # write plasmid contigs to fasta file
    tmp_output_path = output_path.joinpath(prefix + '.plasmid.fasta')
    log.debug('output: plasmids=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        for contig in raw_contigs:
            if(contig['id'] in filtered_contigs):
                fh.write('>' + contig['id'] + '\n')
                fh.write(contig['sequence'] + '\n')


if __name__ == '__main__':
    main()
