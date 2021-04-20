
import json
import logging
import os
import sys

import platon
import platon.config as cfg

log = logging.getLogger('DB')


def check():
    """Check if database directory exists, is accessible and contains necessary files."""

    if(cfg.db_path is None):
        log.error('directory not provided nor detected!')
        sys.exit('ERROR: database directory not provided nor detected! Please provide a valid path to the database directory.')

    if(not os.access(str(cfg.db_path), os.R_OK & os.X_OK)):
        log.error('directory (%s) not readable/accessible!', cfg.db_path)
        sys.exit(f'ERROR: database directory ({cfg.db_path}) not readable/accessible!')

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
        path = cfg.db_path.joinpath(file_name)
        if(not os.access(str(path), os.R_OK) or not path.is_file()):
            log.error('file not readable! file=%s', file_name)
            sys.exit(f'ERROR: database file ({file_name}) not readable!')
