
import logging
import os
import sys
import tempfile
import subprocess as sp
from pathlib import Path

import platon.constants as pc

log = logging.getLogger('functions')


def test_circularity(config, contig):
    """Test if this contig can be circularized."""

    contig_split_position = int(contig['length'] / 2)

    contig_fragment_a_path = config['tmp'].joinpath(contig['id'] + '-a.fasta')
    contig_fragment_a_seq = contig['sequence'][:contig_split_position]
    with contig_fragment_a_path.open(mode='w') as fh:
        fh.write('>a\n')
        fh.write(contig_fragment_a_seq + '\n ')

    contig_fragment_b_path = config['tmp'].joinpath(contig['id'] + '-b.fasta')
    contig_fragment_b_seq = contig['sequence'][contig_split_position:]
    with contig_fragment_b_path.open(mode='w') as fh:
        fh.write('>b\n')
        fh.write(contig_fragment_b_seq + '\n ')
    log.debug(
        'circularity: contig=%s, len=%d, seq-a-len=%d, seq-b-len=%d',
        contig['id'], contig['length'], len(contig_fragment_a_seq), len(contig_fragment_b_seq)
    )

    cmd = [
        'nucmer',
        '-f',  # only forward strand
        '-l', '40',  # increase min match length to 40 bp
        '--threads=1',
        '-p', contig['id'],
        str(contig_fragment_b_path),
        str(contig_fragment_a_path)
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
        log.warning(
            'circularity failed! contig=%s, nucmer-error-code=%d',
            contig['id'], proc.returncode
        )
        log.debug(
            'circularity: contig=%s, cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            contig['id'], cmd, proc.stdout, proc.stderr
        )
        return

    has_match = False
    with config['tmp'].joinpath(contig['id'] + '.delta').open() as fh:
        for line in fh:
            line = line.rstrip()
            if(line[0] == '>'):
                has_match = True
            elif(has_match):
                cols = line.split(' ')
                if(len(cols) == 7):
                    start_b = int(cols[0])
                    end_b = int(cols[1])
                    start_a = int(cols[2])
                    end_a = int(cols[3])
                    mismatches = int(cols[4])
                    alignment_a = end_a - start_a + 1
                    alignment_b = end_b - start_b + 1
                    if(alignment_a == alignment_b
                            and alignment_a > pc.MIN_CIRC_BASEPAIR_OVERLAP
                            and (mismatches / alignment_a) < 0.05
                            and end_b == len(contig_fragment_b_seq)
                            and start_a == 1):
                        contig['is_circular'] = True
                        link = {
                            'length': alignment_a,
                            'mismatches': mismatches,
                            'prime5Start': 1,
                            'prime5End': alignment_a,
                            'prime3Start': contig['length'] - alignment_b + 1,
                            'prime3End': contig['length'],
                        }
                        contig['circular_link'] = link
                        log.info(
                            'circularity: link! id=%s, length=%d, # mismatches=%d, linking-region=[1-%d]...[%d-%d]',
                            contig['id'], link['length'], link['mismatches'], link['prime5End'], link['prime3Start'], link['prime3End']
                        )
                        break
    log.info('circularity: contig=%s, is-circ=%s', contig['id'], contig['is_circular'])
    return


def search_inc_type(config, contig):
    """Search for incompatibility motifs."""

    contig_path = config['tmp'].joinpath(contig['id'] + '.fasta')
    tmp_output_path = config['tmp'].joinpath(contig['id'] + '.inc.blast.out')

    cmd = [
        'blastn',
        '-query', str(config['db'].joinpath('inc-types.fasta')),
        '-subject', str(contig_path),
        '-num_threads', '1',
        '-perc_identity', '90',
        '-culling_limit', '1',
        '-outfmt', '6 qseqid sstart send sstrand pident qcovs bitscore',
        '-out', str(tmp_output_path)
    ]
    proc = sp.run(
        cwd=str(config['tmp']),
        env=config['env'],
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.warning(
            'inc-type failed! contig=%s, blastn-error-code=%d',
            contig['id'], proc.returncode
        )
        log.debug(
            'inc-type: contig=%s, cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            contig['id'], cmd, proc.stdout, proc.stderr
        )
        return

    hits_per_pos = {}
    with tmp_output_path.open() as fh:
        for line in fh:
            cols = line.rstrip().split('\t')
            hit = {
                'type': cols[0],
                'start': int(cols[1]),
                'end': int(cols[2]),
                'strand': '+' if cols[3] == 'plus' else '-',
                'identity': float(cols[4]) / 100,
                'coverage': float(cols[5]) / 100,
                'bitscore': int(cols[6])
            }
            if(hit['coverage'] >= 0.6):
                hit_pos = hit['end'] if hit['strand'] == '+' else hit['start']
                if(hit_pos in hits_per_pos):
                    former_hit = hits_per_pos[hit_pos]
                    if(hit['bitscore'] > former_hit['bitscore']):
                        hits_per_pos[hit_pos] = hit
                        log.info(
                            'inc-type: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                            contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                        )
                else:
                    hits_per_pos[hit_pos] = hit
                    log.info(
                        'inc-type: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )

    contig['inc_types'] = list(hits_per_pos.values())
    log.info('inc-type: contig=%s, # inc-types=%s', contig['id'], len(contig['inc_types']))
    return


def search_rrnas(config, contig):
    """Search for ribosomal RNA sequences."""

    contig_path = config['tmp'].joinpath(contig['id'] + '.fasta')
    tmp_output_path = config['tmp'].joinpath(contig['id'] + '.rrna.cmscan.tsv')

    cmd = [
        'cmscan',
        '--noali',
        '--cut_tc',
        '--cpu', '1',
        '--tblout', str(tmp_output_path),
        str(config['db'].joinpath('rRNA')),
        str(contig_path)
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
        log.warning(
            'rRNAs failed! contig=%s, cmscan-error-code=%d',
            contig['id'], proc.returncode
        )
        log.debug(
            'rRNAs: contig=%s, cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            contig['id'], cmd, proc.stdout, proc.stderr
        )
        return

    with tmp_output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.rstrip().split()
                hit = {
                    'type': cols[0],
                    'start': int(cols[7]),
                    'end': int(cols[8]),
                    'strand': cols[9],
                    'bitscore': float(cols[14]),
                    'evalue': float(cols[15])
                }
                contig['rrnas'].append(hit)
                log.info(
                    'rRNAs: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                    contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                )
    log.info('rRNAs: contig=%s, # rRNAs=%s', contig['id'], len(contig['rrnas']))
    return


def search_amr_genes(config, contigs, filteredProteinsPath):
    """Search for AMR genes."""

    tmp_output_path = config['tmp'].joinpath('amr.hmm.out')
    cmd = [
        'hmmsearch',
        '--noali',
        '--cpu', '1',
        '--cut_tc',
        '--tblout', str(tmp_output_path),
        str(config['db'].joinpath('ncbifam-amr')),
        str(filteredProteinsPath)
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
        log.warning('AMRs failed! hmmsearch-error-code=%d', proc.returncode)
        log.debug(
            'AMRs: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        return

    hits = set()
    with tmp_output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.rstrip().split()
                if(cols[0] not in hits):
                    tmp = cols[0].rsplit('_', 1)
                    contig_id = tmp[0]
                    contig = contigs[contig_id]
                    orf = contig['orfs'][tmp[1]]
                    hit = {
                        'type': cols[2],
                        'hmm-id': cols[3],
                        'start': int(orf['start']),
                        'end': int(orf['end']),
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add(cols[0])
                    contig['amr_hits'].append(hit)
                    log.info(
                        'AMRs: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )
    log.info('AMRs: contig=%s, # AMRs=%s', contig['id'], len(contig['amr_hits']))
    return


def search_reference_plasmids(config, contig):
    """Search for reference plasmid hits."""

    # Reduce blastn word size to overcome segmentation faults due to too many
    # HSPs. As filtered contigs are at least 1k bp long, word size cannot be
    # smaller than 10.
    blast_word_size = int(contig['length'] / 100)

    contig_path = config['tmp'].joinpath(contig['id'] + '.fasta')
    tmp_output_path = config['tmp'].joinpath(contig['id'] + '.refplas.blast.out')

    cmd = [
        'blastn',
        '-query', str(contig_path),
        '-db', str(config['db'].joinpath('refseq-plasmids')),
        '-num_threads', '1',
        '-culling_limit', '1',
        '-perc_identity', '80',
        '-word_size', str(blast_word_size),
        '-outfmt', '6 sseqid qstart qend sstart send slen length nident',
        '-out', str(tmp_output_path)
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
        log.warning(
            'ref plasmids failed! contig=%s, blastn-error-code=%d',
            contig['id'], proc.returncode
        )
        log.debug(
            'ref plasmids: id=%s, cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            contig['id'], cmd, proc.stdout, proc.stderr
        )
        return

    with tmp_output_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            cols = line.split('\t')
            hit = {
                'contig_start': int(cols[1]),
                'contig_end': int(cols[2]),
                'plasmid_start': int(cols[3]),
                'plasmid_end': int(cols[4]),
                'plasmid': {
                    'id': cols[0],
                    'length': int(cols[5])
                },
                'coverage': float(cols[6]) / contig['length'],
                'identity': float(cols[7]) / float(cols[6])
            }
            if(hit['coverage'] >= 0.8):
                contig['plasmid_hits'].append(hit)
                log.info(
                    'ref plasmids: hit! contig=%s, id=%s, c-start=%d, c-end=%d, coverage=%f, identity=%f',
                    contig['id'], hit['plasmid']['id'], hit['contig_start'], hit['contig_end'], hit['coverage'], hit['identity']
                )
    log.info('ref plasmids: contig=%s, # ref plasmids=%s', contig['id'], len(contig['plasmid_hits']))
    return


def search_orit_sequences(config, contig):
    """Search for oriT sequence hits."""

    contig_path = config['tmp'].joinpath(contig['id'] + '.fasta')
    tmp_output_path = config['tmp'].joinpath(contig['id'] + '.orit.blast.out')

    cmd = [
        'blastn',
        '-query', str(contig_path),
        '-db', str(config['db'].joinpath('orit')),
        '-num_threads', '1',
        '-culling_limit', '1',
        '-perc_identity', '90',
        '-evalue', '1E-5',
        '-outfmt', '6 sseqid qstart qend sstart send slen length nident',
        '-out', str(tmp_output_path)
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
        log.warning(
            'oriT failed! contig=%s, blastn-error-code=%d',
            contig['id'], proc.returncode
        )
        log.debug(
            'oriT: id=%s, cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            contig['id'], cmd, proc.stdout, proc.stderr
        )
        return

    with tmp_output_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            cols = line.split('\t')
            hit = {
                'contig_start': int(cols[1]),
                'contig_end': int(cols[2]),
                'orit_start': int(cols[3]),
                'orit_end': int(cols[4]),
                'orit': {
                    'id': cols[0],
                    'length': int(cols[5])
                },
                'coverage': float(cols[6]) / int(cols[5]),
                'identity': float(cols[7]) / float(cols[6])
            }
            if(hit['coverage'] >= 0.9 and hit['identity'] >= 0.9):
                contig['orit_hits'].append(hit)
                log.info(
                    'oriT: hit! contig=%s, id=%s, c-start=%d, c-end=%d, coverage=%f, identity=%f',
                    contig['id'], hit['orit']['id'], hit['contig_start'], hit['contig_end'], hit['coverage'], hit['identity']
                )
    log.info('oriT: contig=%s, # oriT=%s', contig['id'], len(contig['orit_hits']))
    return


def filter_contig(contig):
    """Apply heuristic filters based on contig information."""

    # include all circular contigs
    if(contig['is_circular']):
        log.debug('filter: is circ! contig=%s', contig['id'])
        return True

    # include all contigs with Inc type signatures
    if(len(contig['inc_types']) > 0):
        log.debug('filter: has inc types! contig=%s', contig['id'])
        return True

    # include all contigs containing replication genes
    if(len(contig['replication_hits']) > 0):
        log.debug('filter: has rep hits! contig=%s', contig['id'])
        return True

    # include all contigs containing mobilization genes
    if(len(contig['mobilization_hits']) > 0):
        log.debug('filter: has mob hits! contig=%s', contig['id'])
        return True

    # include all contigs containing oriT sequences
    if(len(contig['orit_hits']) > 0):
        log.debug('filter: has oriT hits! contig=%s', contig['id'])
        return True

    # include all contigs with high confidence protein scores
    if(contig['protein_score'] > pc.RDS_SPECIFICITY_THRESHOLD):
        log.debug('filter: RDS > SPT! contig=%s', contig['id'])
        return True

    # include all contigs with mediocre protein scores but additional blast hit evidence without rRNAs
    if(contig['protein_score'] > pc.RDS_CONSERVATIVE_THRESHOLD
            and len(contig['plasmid_hits']) > 0
            and len(contig['rrnas']) == 0):
        log.debug('filter: RDS > CT & plasmid hits & no rRNAs! contig=%s', contig['id'])
        return True

    return False


def predict_orfs(config, contigs, filteredDraftGenomePath):
    """Predict open reading frames with Prodigal."""

    proteins_path = config['tmp'].joinpath('proteins.faa')
    gff_path = config['tmp'].joinpath('prodigal.gff')
    cmd = [
        'prodigal',
        '-i', str(filteredDraftGenomePath),
        '-a', str(proteins_path),
        '-c',  # closed ends
        '-f', 'gff',  # GFF output
        '-o', str(gff_path)  # prodigal output
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
        log.warning(
            'ORFs failed! prodigal-error-code=%d', proc.returncode
        )
        log.debug(
            'ORFs: cmd=%s stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        return

    # parse orfs
    with gff_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.split('\t')
                orf_id = cols[8].split(';')[0].split('=')[1].split('_')[1]
                orf = {
                    'start': int(cols[3]),
                    'end': int(cols[4]),
                    'strand': cols[6],
                    'id': orf_id
                }
                contig = contigs.get(cols[0], None)
                if(contig is not None):
                    contig['orfs'][orf_id] = orf
    return proteins_path


def search_replication_genes(config, contigs, filteredProteinsPath):
    """Search for replication genes, e.g. repA by custom HMMs."""

    tmp_output_path = config['tmp'].joinpath('rep.hmm.out')
    cmd = [
        'hmmsearch',
        '--noali',
        '--cpu', '1',
        '-E', '1E-100',
        '--tblout', str(tmp_output_path),
        str(config['db'].joinpath('replication')),
        str(filteredProteinsPath)
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
        log.warning(
            'rep genes failed! hmmsearch-error-code=%d', proc.returncode
        )
        log.debug(
            'rep genes: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        return

    hits = set()
    with tmp_output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.rstrip().split()
                if(cols[0] not in hits):
                    tmp = cols[0].rsplit('_', 1)
                    contig_id = tmp[0]
                    contig = contigs[contig_id]
                    orf = contig['orfs'][tmp[1]]
                    hit = {
                        'type': cols[2],
                        'start': int(orf['start']),
                        'end': int(orf['end']),
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add(cols[0])
                    contig['replication_hits'].append(hit)
                    log.info(
                        'rep genes: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )
    log.info('rep genes: contig=%s, # rep-genes=%s', contig['id'], len(contig['replication_hits']))
    return


def search_mobilization_genes(config, contigs, filteredProteinsPath):
    """Search for mobilization genes, e.g. mob by custom HMMs."""

    tmp_output_path = config['tmp'].joinpath('mob.hmm.out')
    cmd = [
        'hmmsearch',
        '--noali',
        '--cpu', '1',
        '-E', '1E-100',
        '--tblout', str(tmp_output_path),
        str(config['db'].joinpath('mobilization')),
        str(filteredProteinsPath)
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
        log.warning(
            'mob genes failed! hmmsearch-error-code=%d', proc.returncode
        )
        log.debug(
            'mob genes: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        return

    hits = set()
    with tmp_output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.rstrip().split()
                if(cols[0] not in hits):
                    tmp = cols[0].rsplit('_', 1)
                    contig_id = tmp[0]
                    contig = contigs[contig_id]
                    orf = contig['orfs'][tmp[1]]
                    hit = {
                        'type': cols[2],
                        'start': int(orf['start']),
                        'end': int(orf['end']),
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add(cols[0])
                    contig['mobilization_hits'].append(hit)
                    log.info(
                        ' mob genes: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )
    log.info('mob genes: contig=%s, # mob-genes=%s', contig['id'], len(contig['mobilization_hits']))
    return


def search_conjugation_genes(config, contigs, filteredProteinsPath):
    """Search for conjugation genes by custom HMMs."""

    tmp_output_path = config['tmp'].joinpath('conj.hmm.out')
    cmd = [
        'hmmsearch',
        '--noali',
        '--cpu', '1',
        '-E', '1E-100',
        '--tblout', str(tmp_output_path),
        str(config['db'].joinpath('conjugation')),
        str(filteredProteinsPath)
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
        log.warning(
            'conj genes failed! hmmsearch-error-code=%d', proc.returncode
        )
        log.debug(
            'conj genes: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        return

    hits = set()
    with tmp_output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.rstrip().split()
                if(cols[0] not in hits):
                    tmp = cols[0].rsplit('_', 1)
                    contig_id = tmp[0]
                    contig = contigs[contig_id]
                    orf = contig['orfs'][tmp[1]]
                    hit = {
                        'type': cols[2],
                        'start': int(orf['start']),
                        'end': int(orf['end']),
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add(cols[0])
                    contig['conjugation_hits'].append(hit)
                    log.info(
                        'conj genes: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )
    log.info('conj genes: contig=%s, # conj-genes=%s', contig['id'], len(contig['conjugation_hits']))
    return


def setup_configuration():
    """Test environment and build a runtime configuration."""

    config = {
        'env': os.environ.copy(),
        'tmp': Path(tempfile.mkdtemp()),
        'bundled-binaries': False
    }
    base_dir = Path(__file__).parent.parent
    share_dir = base_dir.joinpath('share')
    log.debug('config: base-dir=%s', base_dir)
    log.debug('config: share-dir=%s', share_dir)
    if(os.access(str(share_dir), os.R_OK & os.X_OK)):
        config['env']["PATH"] = str(share_dir) + ':' + config['env']["PATH"]
        config['bundled-binaries'] = True
    log.debug('config: bundled binaries=%s', config['bundled-binaries'])

    db_path = base_dir.joinpath('db')
    if(os.access(str(db_path), os.R_OK & os.X_OK)):
        config['db'] = db_path
        log.debug('config: use bundled db, db-path=%s', db_path)
    return config


def test_database(config):
    """Test if database directory exists, is accessible and contains necessary files."""

    if('db' not in config):
        log.error('database directory not detected!')
        sys.exit('ERROR: database directory not detected! Please provide a valid path to the database directory.')

    if(not os.access(str(config['db']), os.R_OK & os.X_OK)):
        log.error('database directory (%s) not readable/accessible!', config['db'])
        sys.exit('ERROR: database directory (%s) not readable/accessible!' % config['db'])

    file_names = [
        'conjugation.h3f',
        'conjugation.h3i',
        'conjugation.h3m',
        'conjugation.h3p',
        'mobilization.h3f',
        'mobilization.h3i',
        'mobilization.h3m',
        'mobilization.h3p',
        'orit.nhr',
        'orit.nin',
        'orit.nsq',
        'replication.h3f',
        'replication.h3i',
        'replication.h3m',
        'replication.h3p',
        'ncbifam-amr.h3f',
        'ncbifam-amr.h3i',
        'ncbifam-amr.h3m',
        'ncbifam-amr.h3p',
        'ncbifam-amr.tsv',
        'rRNA.i1f',
        'rRNA.i1i',
        'rRNA.i1m',
        'rRNA.i1p',
        'mps.tsv',
        'mps.dmnd',
        'refseq-plasmids.tsv',
        'inc-types.fasta',
        'refseq-plasmids.nhr',
        'refseq-plasmids.nin',
        'refseq-plasmids.nsq'
    ]

    for file_name in file_names:
        path = config['db'].joinpath(file_name)
        if(not os.access(str(path), os.R_OK) or not path.is_file()):
            log.error('database file not readable! file=%s', file_name)
            sys.exit('ERROR: database file (%s) not readable!' % file_name)
    return


def test_binaries():
    """Test the proper installation of necessary 3rd party executables."""

    # test prodigal
    try:
        sp.check_call(
            ['prodigal', '-v'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('prodigal not found!')
        sys.exit('ERROR: \'prodigal\' not executable!')
    except:
        pass

    # test diamond
    try:
        sp.check_call(
            ['diamond', 'version'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('diamond not found!')
        sys.exit('ERROR: \'diamond\' not executable!')
    except:
        pass

    # test blastn
    try:
        sp.check_call(
            ['blastn', '-version'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('blastn not found!')
        sys.exit('ERROR: \'blastn\' not executable!')
    except:
        pass

    # test hmmsearch
    try:
        sp.check_call(
            ['hmmsearch', '-h'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('hmmsearch not found!')
        sys.exit('ERROR: \'hmmsearch\' not executable!')
    except:
        pass

    # test nucmer
    try:
        sp.check_call(
            ['nucmer', '-V'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('nucmer not found!')
        sys.exit('ERROR: \'nucmer\' not executable!')
    except:
        pass

    # test cmscan
    try:
        sp.check_call(
            ['cmscan', '-h'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        log.exception('cmscan not found!')
        sys.exit('ERROR: \'cmscan\' not executable!')
    except:
        pass

    return
