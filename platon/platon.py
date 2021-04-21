import functools as ft
import json
import logging
import os
import re
import sys
import shutil
import subprocess as sp
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from Bio import SeqIO

import platon
import platon.db as db
import platon.config as cfg
import platon.constants as pc
import platon.functions as pf
import platon.utils as pu


def main():
    # parse arguments
    args = pu.parse_arguments()

    ############################################################################
    # Setup logging
    ############################################################################
    cfg.prefix = args.prefix if args.prefix else Path(args.genome).stem
    try:
        output_path = Path(args.output) if args.output else Path.cwd()
        if(not output_path.exists()):
            output_path.mkdir(parents=True, exist_ok=True)
        elif(not os.access(str(output_path), os.X_OK)):
            sys.exit(f'ERROR: output path ({output_path}) not accessible!')
        elif(not os.access(str(output_path), os.W_OK)):
            sys.exit(f'ERROR: output path ({output_path}) not writable!')
        output_path = output_path.resolve()
        cfg.output_path = output_path
    except:
        sys.exit(f'ERROR: could not resolve or create output directory ({args.output})!')
    logging.basicConfig(
        filename=str(output_path.joinpath(f'{cfg.prefix}.log')),
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO
    )
    log = logging.getLogger('MAIN')
    log.info('version %s', platon.__version__)
    log.info('command line: %s', ' '.join(sys.argv))

    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################
    cfg.setup(args)  # check parameters and prepare global configuration
    db.check()  # test database
    pu.test_dependencies()  # test dependencies
    
    if(cfg.verbose):
        print(f'Platon v{platon.__version__}')
        print('Options and arguments:')
        print(f'\tinput: {cfg.genome_path}')
        print(f"\tdb: {cfg.db_path}")
        print(f'\toutput: {cfg.output_path}')
        print(f'\tprefix: {cfg.prefix}')
        print(f'\tmode: {cfg.mode}')
        print(f'\tcharacterize: {cfg.characterize}')
        print(f'\ttmp path: {cfg.tmp_path}')
        print(f'\t# threads: {cfg.threads}')

    # parse draft genome
    if(args.verbose):
        print('parse draft genome...')
    contigs = {}
    raw_contigs = []
    try:
        for record in SeqIO.parse(str(cfg.genome_path), 'fasta'):
            length = len(record.seq)
            contig = {
                'id': record.id,
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
            match_spades = re.fullmatch(pc.SPADES_CONTIG_PATTERN, record.id)
            match_unicycler = re.fullmatch(pc.UNICYCLER_CONTIG_PATTERN, record.description)
            if(match_spades is not None):
                contig['coverage'] = float(match_spades.group(1))
            elif(match_unicycler is not None):
                contig['coverage'] = float(match_unicycler.group(1))
                if(match_unicycler.group(2) is not None):
                    contig['is_circular'] = True
            else:
                contig['coverage'] = 0

            # only include contigs with reasonable lengths except of
            # platon runs in characterization mode
            if(args.characterize):
                contigs[record.id] = contig
            else:
                if(length < pc.MIN_CONTIG_LENGTH):
                    log.info('exclude contig: too short: id=%s, length=%d', record.id, length)
                    if (args.verbose):
                        print(f'\texclude contig \'{record.id}\', too short ({length})')
                elif(length >= pc.MAX_CONTIG_LENGTH):
                    log.info('exclude contig: too long: id=%s, length=%d', record.id, length)
                    if (args.verbose):
                        print(f'\texclude contig \'{record.id}\', too long ({length})')
                else:
                    contigs[record.id] = contig
    except:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')

    if(args.verbose):
        print(f'\tparsed {len(raw_contigs)} raw contigs')
        print(f'\texcluded {len(raw_contigs) - len(contigs)} contigs by size filter')
        print(f'\tanalyze {len(contigs)} contigs')
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
    proteins_path = pf.predict_orfs(contigs, cfg.genome_path)
    if(proteins_path is None):
        sys.exit('Error: ORF prediction failed!')
    no_orfs = ft.reduce(lambda x, y: x + y, map(lambda k: len(contigs[k]['orfs']), contigs))
    log.info('ORF detection: # ORFs=%d', no_orfs)
    if(args.verbose):
        print(f'\tfound {no_orfs} ORFs')

    # retain / exclude contigs without ORFs
    if(args.characterize or args.mode == 'sensitivity' or args.mode == 'accuracy'):
        log.info('ORF contig filter disabled! # passed contigs=%s', len(contigs))
    else:  # exclude contigs without ORFs in specificity mode
        tmp_contigs = {}
        for contig in filter(lambda k: len(k['orfs']) > 0, contigs.values()):
            tmp_contigs[contig['id']] = contig
        no_removed_contigs = len(contigs) - len(tmp_contigs)
        if(args.verbose):
            print(f'\tdiscarded {no_removed_contigs} contig(s) without ORFs')
        contigs = tmp_contigs
        log.info('ORF contig filter: # discarded=%s, # remaining=%s', no_removed_contigs, len(contigs))

    # find marker genes
    if(args.verbose):
        print('search marker protein sequences (MPS)...')
    tmp_output_path = cfg.tmp_path.joinpath('diamond.tsv')
    cmd = [
        'diamond',
        'blastp',
        '--db', str(cfg.db_path.joinpath('mps.dmnd')),
        '--query', str(proteins_path),
        '--out', str(tmp_output_path),
        '--max-target-seqs', '1',  # max 1 result per query
        '--id', '90',  # min alignment identity 90%
        '--query-cover', '80',  # min query cov 80%
        '--subject-cover', '80',  # min subjetc cov 80%
        '--threads', str(args.threads),  # threads
        '--tmpdir', str(cfg.tmp_path)
    ]
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
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
        print(f'\tfound {proteins_identified} MPS')

    # parse protein score file
    if(args.verbose):
        print('compute replicon distribution scores (RDS)...')
    marker_proteins = {}
    with cfg.db_path.joinpath('mps.tsv').open() as fh:
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
            print(f'apply RDS sensitivity threshold (SNT={pc.RDS_SENSITIVITY_THRESHOLD:2.1f}) filter...')
        scored_contigs = {k: v for (k, v) in contigs.items() if v['protein_score'] >= pc.RDS_SENSITIVITY_THRESHOLD}
        no_excluded_contigs = len(contigs) - len(scored_contigs)
        log.info('RDS SNT filter: # discarded contigs=%d, # remaining contigs=%d', no_excluded_contigs, len(scored_contigs))
        if(args.verbose):
            print(f'\texcluded {no_excluded_contigs} contigs by SNT filter')
            print('characterize contigs...')

    # extract proteins from potential plasmid contigs for subsequent analyses
    filtered_proteins_path = cfg.tmp_path.joinpath('proteins-filtered.faa')
    with filtered_proteins_path.open(mode='w') as fh:
        for record in SeqIO.parse(str(proteins_path), 'fasta'):
            orf_name = str(record.id).split()[0]
            contig_id = orf_name.rsplit('_', 1)[0]
            if(contig_id in scored_contigs):
                fh.write(f'>{orf_name}\n')
                fh.write(f'{record.seq}\n')

    # write contig sequences to fasta files for subsequent parallel analyses
    for id, contig in scored_contigs.items():
        contig_path = cfg.tmp_path.joinpath(f"{contig['id']}.fasta")
        with contig_path.open(mode='w') as fh:
            fh.write(f">{contig['id']}\n")
            fh.write(f"{contig['sequence']}\n")

    # init thread pool to parallize io bound analyses
    with ThreadPoolExecutor(max_workers=args.threads) as tpe:
        # start search for replication, mobilization and conjugation genes
        for fn in (pf.search_replication_genes, pf.search_mobilization_genes, pf.search_conjugation_genes, pf.search_amr_genes):
            tpe.submit(fn, scored_contigs, filtered_proteins_path)
        # analyse contigs
        for id, contig in scored_contigs.items():
            tpe.submit(pf.test_circularity, contig)
            tpe.submit(pf.search_inc_type, contig)
            tpe.submit(pf.search_rrnas, contig)
            tpe.submit(pf.search_orit_sequences, contig)
            tpe.submit(pf.search_reference_plasmids, contig)

    # lookup AMR genes
    amr_genes = {}
    with cfg.db_path.joinpath('ncbifam-amr.tsv').open() as fh:
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
    shutil.rmtree(str(cfg.tmp_path))
    log.debug('removed tmp dir: %s', cfg.tmp_path)

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
    tmp_output_path = output_path.joinpath(f'{cfg.prefix}.tsv')
    log.debug('output: tsv=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        fh.write(pc.HEADER + '\n')
        for id in sorted(filtered_contigs, key=lambda k: -filtered_contigs[k]['length']):
            c = filtered_contigs[id]
            cov = 'NA' if c['coverage'] == 0 else f"{c['coverage']:4.1f}"
            circ = 'yes' if c['is_circular'] else 'no'
            line = f"{c['id']}\t{c['length']}\t{cov}\t{len(c['orfs'])}\t{c['protein_score']:3.1f}\t{circ}\t{len(c['inc_types'])}\t{len(c['replication_hits'])}\t{len(c['mobilization_hits'])}\t{len(c['orit_hits'])}\t{len(c['conjugation_hits'])}\t{len(c['amr_hits'])}\t{len(c['rrnas'])}\t{len(c['plasmid_hits'])}"
            print(line)
            log.info(
                'plasmid: id=%s, len=%d, cov=%s, ORFs=%d, RDS=%f, circ=%s, incs=%s, # rep=%d, # mob=%d, #oriT=%d, # con=%d, # AMRs=%d, # rRNAs=%d, # refs=%d',
                c['id'], c['length'], cov, len(c['orfs']), c['protein_score'], c['is_circular'], len(c['inc_types']), len(c['replication_hits']), len(c['mobilization_hits']), len(c['orit_hits']), len(c['conjugation_hits']), len(c['amr_hits']), len(c['rrnas']), len(c['plasmid_hits'])
            )
            fh.write(f'{line}\n')

    # write comprehensive results to JSON file
    tmp_output_path = output_path.joinpath(f'{cfg.prefix}.json')
    log.debug('output: json=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        indent = '\t' if args.verbose else None
        separators = (', ', ': ') if args.verbose else (',', ':')
        json.dump(filtered_contigs, fh, indent=indent, separators=separators)

    # write chromosome contigs to fasta file
    tmp_output_path = output_path.joinpath(f'{cfg.prefix}.chromosome.fasta')
    log.debug('output: chromosomes=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        for contig in raw_contigs:
            if(contig['id'] not in filtered_contigs):
                fh.write(f">{contig['id']}\n{contig['sequence']}\n")

    # write plasmid contigs to fasta file
    tmp_output_path = output_path.joinpath(f'{cfg.prefix}.plasmid.fasta')
    log.debug('output: plasmids=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        for contig in raw_contigs:
            if(contig['id'] in filtered_contigs):
                fh.write(f">{contig['id']}\n{contig['sequence']}\n")


if __name__ == '__main__':
    main()
