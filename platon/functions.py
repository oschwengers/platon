
import os
import sys
import tempfile
import subprocess as sp

import platon.constants as pc


def test_circularity(config, contig):
    """Test if this contig can be circularized."""

    contig_split_position = int(contig['length'] / 2)

    contig_fragment_a_path = config['tmp'] + '/' + contig['id'] + '-a.fasta'
    with open(contig_fragment_a_path, 'w') as fh:
        fh.write('>a\n')
        fh.write(contig['sequence'][:contig_split_position] + '\n ')

    contig_fragment_b_path = config['tmp'] + '/' + contig['id'] + '-b.fasta'
    contig_fragment_b_seq = contig['sequence'][contig_split_position:]
    with open(contig_fragment_b_path, 'w') as fh:
        fh.write('>b\n ')
        fh.write(contig_fragment_b_seq + '\n ')

    sp.check_call(
        [
            'nucmer',
            '-f',  # only forward strand
            '-l', '40',  # increase min match length to 40 bp
            '--threads=1',
            '-p', contig['id'],
            contig_fragment_b_path,
            contig_fragment_a_path
        ],
        cwd=config['tmp'],
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT
    )

    contig['is_circular'] = False
    has_match = False
    with open(config['tmp'] + '/' + contig['id'] + '.delta', 'r') as fh:
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
                            and alignment_a > 50
                            and (mismatches/alignment_a) < 0.05
                            and end_b == len(contig_fragment_b_seq)
                            and start_a == 1):
                        contig['is_circular'] = True
                        contig['circular_link'] = {
                            'length': alignment_a,
                            'mismatches': mismatches,
                            'prime5Start': 1,
                            'prime5End': alignment_a,
                            'prime3Start': contig['length'] - alignment_b + 1,
                            'prime3End': contig['length'],
                        }
                        break
    return


def search_inc_type(config, contig):
    """Search for incompatibility motifs."""

    contig_path = config['tmp'] + '/' + contig['id'] + '.fasta'
    tmp_output_path = config['tmp'] + '/' + contig['id'] + '.inc.blast.out'
    sp.check_call(
        [
            'blastn',
            '-query', config['db'] + '/inc-types.fasta',
            '-subject', contig_path,
            '-num_threads', '1',
            '-perc_identity', '90',
            '-culling_limit', '1',
            '-outfmt', '6 qseqid sstart send sstrand pident qcovs bitscore',
            '-out', tmp_output_path
        ],
        cwd=config['tmp'],
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT
    )

    hits_per_pos = {}
    with open(tmp_output_path, 'r') as fh:
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
                else:
                    hits_per_pos[hit_pos] = hit

    contig['inc_types'] = list(hits_per_pos.values())
    return


def search_rrnas(config, contig):
    """Search for ribosomal RNA sequences."""

    contig_path = config['tmp'] + '/' + contig['id'] + '.fasta'
    tmp_output_path = config['tmp'] + '/' + contig['id'] + '.rrna.cmscan.tsv'
    sp.check_call(
        [
            'cmscan',
            '--noali',
            '--cut_tc',
            '--cpu', '1',
            '--tblout', tmp_output_path,
            config['db'] + '/rRNA',
            contig_path,
        ],
        cwd=config['tmp'],
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT
    )

    with open(tmp_output_path, 'r') as fh:
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
    return


def search_amr_genes(config, contigs, filteredProteinsPath):
    """Search for AMR genes."""

    tmp_output_path = config['tmp'] + '/amr.hmm.out'
    sp.check_call(
        [
            'hmmsearch',
            '--noali',
            '--cpu', '1',
            '--cut_tc',
            '--tblout', tmp_output_path,
            config['db'] + '/ncbifam-amr',
            filteredProteinsPath
        ],
        cwd=config['tmp'],
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT
    )

    hits = set()
    with open(tmp_output_path, 'r') as fh:
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
                        'start': orf['start'],
                        'end': orf['end'],
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add(cols[0])
                    contig['amr_hits'].append(hit)

    return


def search_reference_plasmids(config, contig):
    """Search for reference plasmid hits."""

    # Reduce blastn word size to overcome segmentation faults due to too many
    # HSPs. As filtered contigs are at least 1k bp long, word size cannot be
    # smaller than 10.
    blast_word_size = int(contig['length'] / 100)

    contig_path = config['tmp'] + '/' + contig['id'] + '.fasta'
    tmp_output_path = config['tmp'] + '/' + contig['id'] + '.refplas.blast.out'
    sp.check_call(
        [
            'blastn',
            '-query', contig_path,
            '-db', config['db'] + '/refseq-plasmids',
            '-num_threads', '1',
            '-culling_limit', '1',
            '-perc_identity', '80',
            '-word_size', str(blast_word_size),
            '-outfmt', '6 sseqid qstart qend sstart send slen length nident',
            '-out', tmp_output_path
        ],
        cwd=config['tmp'],
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT
    )

    with open(tmp_output_path, 'r') as fh:
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
    return


def filter_contig(contig):
    """Apply heuristic filters based on contig information."""

    # include all circular contigs
    if(contig['is_circular']):
        return True

    # include all contigs with Inc type signatures
    if(len(contig['inc_types']) > 0):
        return True

    # include all contigs containing replication genes
    if(len(contig['replication_hits']) > 0):
        return True

    # include all contigs containing mobilization genes
    if(len(contig['mobilization_hits']) > 0):
        return True

    # include all contigs with high confidence protein scores
    if(contig['protein_score'] > pc.RDS_SPECIFICITY_THRESHOLD):
        return True

    # include all contigs with mediocre protein scores but additional blast hit evidence without rRNAs
    if(contig['protein_score'] > pc.RDS_CONSERVATIVE_THRESHOLD
            and len(contig['plasmid_hits']) > 0
            and len(contig['rrnas']) == 0):
        return True

    return False


def predict_orfs(config, contigs, filteredDraftGenomePath):
    """Predict open reading frames with Prodigal."""

    proteins_path = config['tmp'] + '/proteins.faa'
    gff_path = config['tmp'] + '/prodigal.gff'
    sp.check_call(
        [
            'prodigal',
            '-i', filteredDraftGenomePath,
            '-a', proteins_path,
            '-c',  # closed ends
            '-f', 'gff',  # GFF output
            '-o', gff_path  # prodigal output
        ],
        cwd=config['tmp'],
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT
    )

    # parse orfs
    with open(gff_path, 'rU') as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.split('\t')
                orf_id = cols[8].split(';')[0].split('=')[1].split('_')[1]
                orf = {
                    'start': cols[3],
                    'end': cols[4],
                    'strand': cols[6],
                    'id': orf_id
                }
                contig = contigs.get(cols[0], None)
                if(contig is not None):
                    contig['orfs'][orf_id] = orf
    return proteins_path


def search_replication_genes(config, contigs, filteredProteinsPath):
    """Search for replication genes, e.g. repA by custom HMMs."""

    tmp_output_path = config['tmp'] + '/rep.hmm.out'
    sp.check_call(
        [
            'hmmsearch',
            '--noali',
            '--cpu', '1',
            '-E', '1E-100',
            '--tblout', tmp_output_path,
            config['db'] + '/replication',
            filteredProteinsPath
        ],
        cwd=config['tmp'],
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT
    )

    hits = set()
    with open(tmp_output_path, 'r') as fh:
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
                        'start': orf['start'],
                        'end': orf['end'],
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add(cols[0])
                    contig['replication_hits'].append(hit)

    return


def search_mobilization_genes(config, contigs, filteredProteinsPath):
    """Search for mobilization genes, e.g. mob by custom HMMs."""

    tmp_output_path = config['tmp'] + '/mob.hmm.out'
    sp.check_call(
        [
            'hmmsearch',
            '--noali',
            '--cpu', '1',
            '-E', '1E-100',
            '--tblout', tmp_output_path,
            config['db'] + '/mobilization',
            filteredProteinsPath
        ],
        cwd=config['tmp'],
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT
    )

    hits = set()
    with open(tmp_output_path, 'r') as fh:
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
                        'start': orf['start'],
                        'end': orf['end'],
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add(cols[0])
                    contig['mobilization_hits'].append(hit)

    return


def search_conjugation_genes(config, contigs, filteredProteinsPath):
    """Search for conjugation genes by custom HMMs."""

    tmp_output_path = config['tmp'] + '/conj.hmm.out'
    sp.check_call(
        [
            'hmmsearch',
            '--noali',
            '--cpu', '1',
            '-E', '1E-100',
            '--tblout', tmp_output_path,
            config['db'] + '/conjugation',
            filteredProteinsPath
        ],
        cwd=config['tmp'],
        env=config['env'],
        stdout=sp.DEVNULL,
        stderr=sp.STDOUT
    )

    hits = set()
    with open(tmp_output_path, 'r') as fh:
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
                        'start': orf['start'],
                        'end': orf['end'],
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add(cols[0])
                    contig['conjugation_hits'].append(hit)

    return


def setup_configuration():
    """Test environment and build a runtime configuration."""

    config = {
        'env': os.environ.copy(),
        'tmp': tempfile.mkdtemp(),
        'bundled-binaries': False
    }
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
    share_dir = os.path.join(base_dir, 'share')
    if(os.access(share_dir, os.R_OK & os.X_OK)):
        config['env']["PATH"] = share_dir + ':' + config['env']["PATH"]
        config['bundled-binaries'] = True

    db_path = os.path.join(base_dir, 'db')
    if(os.access(db_path, os.R_OK & os.X_OK)):
        config['db'] = db_path

    return config


def test_database(config):
    """Test if database directory exists, is accessible and contains necessary files."""

    if('db' not in config):
        sys.exit('ERROR: database directory not detected! Please provide a valid path to the database directory.')

    if(not os.access(config['db'], os.R_OK & os.X_OK)):
        sys.exit('ERROR: database directory ('+config['db']+') not readable/accessible!')

    file_names = [
        'conjugation.h3f',
        'conjugation.h3i',
        'conjugation.h3m',
        'conjugation.h3p',
        'mobilization.h3f',
        'mobilization.h3i',
        'mobilization.h3m',
        'mobilization.h3p',
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
        'rds.tsv',
        'refseq-plasmids.tsv',
        'inc-types.fasta',
        'refseq-plasmids.nhr',
        'refseq-plasmids.nin',
        'refseq-plasmids.nsq',
        'refseq-bacteria-nrpc-reps.inf',
        'refseq-bacteria-nrpc-reps_0.inf',
        'refseq-bacteria-nrpc-reps_0.nam',
        'refseq-bacteria-nrpc-reps_0.off',
        'refseq-bacteria-nrpc-reps_0.seq',
        'refseq-bacteria-nrpc-reps_0.src'
    ]

    for file_name in file_names:
        path = config['db'] + '/' + file_name
        if(not os.access(path, os.R_OK) or not os.path.isfile(path)):
            sys.exit('ERROR: database file ('+file_name+') not readable!')
    return


def test_binaries():
    """Test the proper installation of necessary 3rd party executables."""

    # test prodigal
    try:
        sp.check_call(
            [
                'prodigal',
                '-v'
            ],
            stdout=sp.DEVNULL,
            stderr=sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'prodigal\' not executable!')
    except:
        pass

    # test ghostz
    try:
        sp.check_call(
            ['ghostz'],
            stdout=sp.DEVNULL,
            stderr=sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'ghostz\' not executable!')
    except:
        pass

    # test blastn
    try:
        sp.check_call(
            [
                'blastn',
                '-version'
            ],
            stdout=sp.DEVNULL,
            stderr=sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'blastn\' not executable!')
    except:
        pass

    # test hmmsearch
    try:
        sp.check_call(
            [
                'hmmsearch'
                '-h'
            ],
            stdout=sp.DEVNULL,
            stderr=sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'hmmsearch\' not executable!')
    except:
        pass

    # test nucmer
    try:
        sp.check_call(
            [
                'nucmer',
                '-V'
            ],
            stdout=sp.DEVNULL,
            stderr=sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'nucmer\' not executable!')
    except:
        pass

    # test cmscan
    try:
        sp.check_call(
            [
                'cmscan',
                '-h'
            ],
            stdout=sp.DEVNULL,
            stderr=sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'cmscan\' not executable!')
    except:
        pass

    return
