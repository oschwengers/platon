
import os
import sys
import tempfile
import subprocess as sp

from platon.constants import *

# global constants


# test for circularity
def test_circularity( config, contig ):

    contigSplitPosition = int(contig['length'] / 2)

    contigFragmentAPath = config['tmp'] + '/' + contig['id'] + '-a.fasta'
    with open( contigFragmentAPath, 'w' ) as fh:
        fh.write( '>a\n')
        fh.write( contig['sequence'][:contigSplitPosition] + '\n ')

    contigFragmentBPath = config['tmp'] + '/' + contig['id'] + '-b.fasta'
    contigFragmentBSequence = contig['sequence'][contigSplitPosition:]
    with open( contigFragmentBPath, 'w' ) as fh:
        fh.write( '>b\n ')
        fh.write( contigFragmentBSequence + '\n ')

    sp.check_call( ['nucmer',
            '-f', # only forward strand
            '-l', '40', # increase min match length to 40 bp
            '--threads=1',
            '-p', contig['id'],
            contigFragmentBPath,
            contigFragmentAPath
        ],
        cwd = config['tmp'],
        env = config['env'],
        stdout = sp.DEVNULL,
        stderr = sp.STDOUT
    )

    contig['is_circular'] = False
    hasMatch = False
    with open( config['tmp'] + '/' + contig['id'] + '.delta', 'r' ) as fh:
        for line in fh:
            line = line.rstrip()
            if( line[0] == '>' ):
                hasMatch = True
            elif( hasMatch ):
                cols = line.split( ' ' )
                if( len(cols) == 7 ):
                    bStart = int(cols[0])
                    bEnd   = int(cols[1])
                    aStart = int(cols[2])
                    aEnd   = int(cols[3])
                    mismatches = int(cols[4])
                    aAlignment = aEnd - aStart + 1
                    bAlignment = bEnd - bStart + 1
                    if( aAlignment == bAlignment  and  aAlignment > 50
                        and (mismatches/aAlignment) < 0.05
                        and bEnd == len(contigFragmentBSequence)
                        and aStart == 1):
                            contig['is_circular'] = True
                            contig['circular_link'] = {
                                'length': aAlignment,
                                'mismatches': mismatches,
                                'prime5Start': 1,
                                'prime5End': aAlignment,
                                'prime3Start': contig['length'] - bAlignment + 1,
                                'prime3End': contig['length'],
                            }
                            break
    return




# search for inc type signatures
def search_inc_type( config, contig ):

    contigPath = config['tmp'] + '/' + contig['id'] + '.fasta'
    outPath = config['tmp'] + '/' + contig['id'] + '.inc.blast.out'
    sp.check_call( [ 'blastn',
            '-query', config['db'] + '/inc-types.fasta',
            '-subject', contigPath,
            '-num_threads', '1',
            '-perc_identity', '80',
            '-culling_limit', '1',
            '-outfmt', '6 qseqid sstart send sstrand pident qcovs',
            '-out', outPath
        ],
        cwd = config['tmp'],
        env = config['env'],
        stdout = sp.DEVNULL,
        stderr = sp.STDOUT
    )

    with open( outPath, 'r' ) as fh:
        for line in fh:
            cols = line.rstrip().split( '\t' )
            hit = {
                'type': cols[0],
                'start': int(cols[1]),
                'end': int(cols[2]),
                'strand': '+' if cols[3] == 'plus' else '-',
                'identity': float(cols[4]) / 100,
                'coverage': float(cols[5]) / 100
            }
            if( hit['coverage'] >= 0.9 ):
                contig['inc_types'].append( hit )
    return




# search for rRNAs
def search_rrnas( config, contig ):

    contigPath = config['tmp'] + '/' + contig['id'] + '.fasta'
    outPath = config['tmp'] + '/' + contig['id'] + '.rrna.cmscan.tsv'
    sp.check_call( [ 'cmscan',
            '--noali',
            '--cut_tc',
            '--cpu', '1',
            '--tblout', outPath,
            config['db'] + '/rRNA',
            contigPath,
        ],
        cwd = config['tmp'],
        env = config['env'],
        stdout = sp.DEVNULL,
        stderr = sp.STDOUT
    )

    with open( outPath, 'r' ) as fh:
        for line in fh:
            if( line[0] != '#' ):
                cols = line.rstrip().split()
                hit = {
                    'type': cols[0],
                    'start': int(cols[7]),
                    'end': int(cols[8]),
                    'strand': cols[9],
                    'bitscore': float(cols[14]),
                    'evalue': float(cols[15])
                }
                contig['rrnas'].append( hit )
    return




# search for amr genes
def search_amr_genes( config, contigs, filteredProteinsPath ):

    outPath = config['tmp'] + '/amr.hmm.out'
    sp.check_call( [ 'hmmsearch',
            '--noali',
            '--cpu', '1',
            '--cut_tc',
            '--tblout', outPath,
            config['db'] + '/ncbifam-amr',
            filteredProteinsPath
        ],
        cwd = config['tmp'],
        env = config['env'],
        stdout = sp.DEVNULL,
        stderr = sp.STDOUT
    )

    hits = set()
    with open( outPath, 'r' ) as fh:
        for line in fh:
            if( line[0] != '#' ):
                cols = line.rstrip().split()
                if( not cols[0] in hits ):
                    tmp = cols[0].rsplit('_', 1 )
                    contigId = tmp[0]
                    contig = contigs[ contigId ]
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
                    hits.add( cols[0] )
                    contig['amr_hits'].append(hit)

    return




# search for reference plasmid hits
def search_reference_plasmids( config, contig ):

    # reduce blastn word size to overcome segmentation faults due to too many
    # HSPs. As filtered contigs are at least 1k bp long, word size cannot be
    # smaller than 10.
    blastWordSize = int( contig['length'] / 100 )

    contigPath = config['tmp'] + '/' + contig['id'] + '.fasta'
    outPath = config['tmp'] + '/' + contig['id'] + '.refplas.blast.out'
    sp.check_call( [ 'blastn',
            '-query', contigPath,
            '-db', config['db'] + '/refseq-plasmids',
            '-num_threads', '1',
            '-culling_limit', '1',
            '-perc_identity', '80',
            '-word_size', str(blastWordSize),
            '-outfmt', '6 sseqid qstart qend sstart send slen length nident',
            '-out', outPath
        ],
        cwd = config['tmp'],
        env = config['env'],
        stdout = sp.DEVNULL,
        stderr = sp.STDOUT
    )

    with open( outPath, 'r' ) as fh:
        for line in fh:
            line = line.rstrip()
            cols = line.split( '\t' )
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
            if( hit['coverage'] >= 0.8 ):
                contig['plasmid_hits'].append( hit )
    return




def filter_contig( contig ):

    # include all circular contigs
    if( contig['is_circular'] ):
        return True

    # include all contigs containing Inc type signatures
    if( len(contig['inc_types']) > 0 ):
        return True

    # include all contigs containing replication genes
    if( len(contig['replication_hits']) > 0 ):
        return True

    # include all contigs containing mobilization genes
    if( len(contig['mobilization_hits']) > 0 ):
        return True

    # include all contigs with high confidence protein scores
    if( contig['protein_score'] > PROTEIN_SCORE_TRUSTED_THRESHOLD ):
        return True

    # include all contigs with mediocre protein scores but additional blast hit evidence without rRNAs
    if( contig['protein_score'] > PROTEIN_SCORE_CONSERVATIVE_THRESHOLD  and  len(contig['plasmid_hits']) > 0  and  len(contig['rrnas']) == 0 ):
        return True

    return False




# predict ORFs
def predict_orfs( config, contigs, filteredDraftGenomePath ):
    proteinsPath = config['tmp'] + '/proteins.faa'
    gffPath = config['tmp'] + '/prodigal.gff'
    sp.check_call( [ 'prodigal',
            '-i', filteredDraftGenomePath,
            '-a', proteinsPath,
            '-c', # closed ends
            '-f', 'gff', # GFF output
            '-o', gffPath # prodigal output
        ],
        cwd = config['tmp'],
        env = config['env'],
        stdout = sp.DEVNULL,
        stderr = sp.STDOUT
    )

    # parse orfs
    with open( gffPath, 'rU' ) as fh:
        for line in fh:
            if( line[0] != '#' ):
                cols = line.split( '\t' )
                orfId = cols[8].split(';')[0].split('=')[1].split('_')[1]
                orf = {
                    'start': cols[3],
                    'end': cols[4],
                    'strand': cols[6],
                    'id': orfId
                }
                contig = contigs.get(cols[0], None)
                if( contig is not None ):
                    contig['orfs'][ orfId ] = orf
    return proteinsPath




# search for replication genes
def search_replication_genes( config, contigs, filteredProteinsPath ):

    outPath = config['tmp'] + '/rep.hmm.out'
    sp.check_call( [ 'hmmsearch',
            '--noali',
            '--cpu', '1',
            '-E', '1E-100',
            '--tblout', outPath,
            config['db'] + '/replication',
            filteredProteinsPath
        ],
        cwd = config['tmp'],
        env = config['env'],
        stdout = sp.DEVNULL,
        stderr = sp.STDOUT
    )

    hits = set()
    with open( outPath, 'r' ) as fh:
        for line in fh:
            if( line[0] != '#' ):
                cols = line.rstrip().split()
                if( not cols[0] in hits ):
                    tmp = cols[0].rsplit('_', 1 )
                    contigId = tmp[0]
                    contig = contigs[ contigId ]
                    orf = contig['orfs'][tmp[1]]
                    hit = {
                        'type': cols[2],
                        'start': orf['start'],
                        'end': orf['end'],
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add( cols[0] )
                    contig['replication_hits'].append(hit)

    return




# search for mobilization genes
def search_mobilization_genes( config, contigs, filteredProteinsPath ):

    outPath = config['tmp'] + '/mob.hmm.out'
    sp.check_call( [ 'hmmsearch',
            '--noali',
            '--cpu', '1',
            '-E', '1E-100',
            '--tblout', outPath,
            config['db'] + '/mobilization',
            filteredProteinsPath
        ],
        cwd = config['tmp'],
        env = config['env'],
        stdout = sp.DEVNULL,
        stderr = sp.STDOUT
    )

    hits = set()
    with open( outPath, 'r' ) as fh:
        for line in fh:
            if( line[0] != '#' ):
                cols = line.rstrip().split()
                if( not cols[0] in hits ):
                    tmp = cols[0].rsplit('_', 1 )
                    contigId = tmp[0]
                    contig = contigs[ contigId ]
                    orf = contig['orfs'][tmp[1]]
                    hit = {
                        'type': cols[2],
                        'start': orf['start'],
                        'end': orf['end'],
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add( cols[0] )
                    contig['mobilization_hits'].append(hit)

    return




# search for conjugation genes
def search_conjugation_genes( config, contigs, filteredProteinsPath ):

    outPath = config['tmp'] + '/conj.hmm.out'
    sp.check_call( [ 'hmmsearch',
            '--noali',
            '--cpu', '1',
            '-E', '1E-100',
            '--tblout', outPath,
            config['db'] + '/conjugation',
            filteredProteinsPath
        ],
        cwd = config['tmp'],
        env = config['env'],
        stdout = sp.DEVNULL,
        stderr = sp.STDOUT
    )

    hits = set()
    with open( outPath, 'r' ) as fh:
        for line in fh:
            if( line[0] != '#' ):
                cols = line.rstrip().split()
                if( not cols[0] in hits ):
                    tmp = cols[0].rsplit('_', 1 )
                    contigId = tmp[0]
                    contig = contigs[ contigId ]
                    orf = contig['orfs'][tmp[1]]
                    hit = {
                        'type': cols[2],
                        'start': orf['start'],
                        'end': orf['end'],
                        'strand': orf['strand'],
                        'bitscore': float(cols[5]),
                        'evalue': float(cols[4])
                    }
                    hits.add( cols[0] )
                    contig['conjugation_hits'].append(hit)

    return

def setup_configuration():

    config = {
        'env': os.environ.copy(),
        'tmp': tempfile.mkdtemp(),
        'bundled-binaries': False
    }
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
#    print( '__file__=' + __file__ )
#    print( 'base_dir=' + base_dir )
    share_dir = os.path.join( base_dir, 'share' )
    if( os.access( share_dir, os.R_OK & os.X_OK ) ):
        config['env']["PATH"] = share_dir + ':' + config['env']["PATH"]
        config['bundled-binaries'] = True

    db_path = os.path.join( base_dir, 'db' )
    if( os.access( db_path, os.R_OK & os.X_OK ) ):
        config['db'] = db_path

    return config

def test_binaries():

    # test prodigal
    try:
        sp.check_call( [ 'prodigal',
                '-v'
            ],
            stdout = sp.DEVNULL,
            stderr = sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit( 'ERROR: \'prodigal\' not executable!' )
    except:
        pass

    # test ghostz
    try:
        sp.check_call( ['ghostz'],
            stdout = sp.DEVNULL,
            stderr = sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit( 'ERROR: \'ghostz\' not executable!' )
    except:
        pass

    # test blastn
    try:
        sp.check_call( [ 'blastn',
                '-version'
            ],
            stdout = sp.DEVNULL,
            stderr = sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit( 'ERROR: \'blastn\' not executable!' )
    except:
        pass

    # test hmmsearch
    try:
        sp.check_call( ['hmmsearch'],
            stdout = sp.DEVNULL,
            stderr = sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit( 'ERROR: \'hmmsearch\' not executable!' )
    except:
        pass

    # test nucmer
    try:
        sp.check_call( ['nucmer',
                '-V'
            ],
            stdout = sp.DEVNULL,
            stderr = sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit( 'ERROR: \'nucmer\' not executable!' )
    except:
        pass

    # test cmscan
    try:
        sp.check_call( ['cmscan'],
            stdout = sp.DEVNULL,
            stderr = sp.STDOUT
        )
    except FileNotFoundError:
        sys.exit( 'ERROR: \'cmscan\' not executable!' )
    except:
        pass

    return
