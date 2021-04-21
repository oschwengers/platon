import logging
import os
import sys
import tempfile
from pathlib import Path

log = logging.getLogger('CONFIG')

# runtime configurations
env = os.environ.copy()
threads = None
verbose = None

# input / output configuration
db_path = None
genome_path = None
prefix = None
tmp_path = None
output_path = None

# workflow configuration
mode = None
characterize = None


def setup(args):
    """Test environment and build a runtime configuration."""
    # runtime configurations
    global env, threads, verbose
    threads = args.threads
    log.info('threads=%i', threads)
    verbose = args.verbose
    log.info('verbose=%s', verbose)

    # input / output path configurations
    global db_path, genome_path, prefix, tmp_path, output_path
    if(args.db):
        db_dir = args.db
        log.debug('test parameter db: db_tmp=%s', db_dir)
        try:
            db_tmp_path = Path(db_dir).resolve()
            if(db_tmp_path.is_dir()):
                db_path = db_tmp_path
                log.info('database detected: type=parameter, path=%s', db_path)
            else:
                log.error('unvalid database path: type=parameter, path=%s', db_tmp_path)
                raise IOError()
        except:
            sys.exit(f'ERROR: wrong database path! --db={db_dir}')
    elif('PLATON_DB' in env):
        db_dir = env['PLATON_DB']
        log.debug('test env db: db_tmp=%s', db_dir)
        try:
            db_tmp_path = Path(db_dir).resolve()
            if(db_tmp_path.is_dir()):
                db_path = db_tmp_path
                log.info('database detected: type=environment, path=%s', db_path)
            else:
                log.error('unvalid database path: type=environment, path=%s', db_tmp_path)
                raise IOError()
        except:
            sys.exit(f'ERROR: wrong database path! PLATON_DB={db_dir}')
    else:
        base_dir = Path(__file__).parent.parent
        db_tmp_path = base_dir.joinpath('db')
        log.debug('test base_dir db: db_tmp=%s', db_tmp_path)
        if(db_tmp_path.is_dir()):
            db_path = db_tmp_path
            log.info('database detected: type=base-dir, path=%s', db_path)
        else:
            log.error('unvalid database path: type=base-dir, path=%s', db_tmp_path)
            sys.exit('ERROR: database neither auto-detected nor provided!\nPlease, download the mandatory db and provide it either via the --db parameter, via a PLATON_DB environment variable or copy it into the Platon base directory.\nFor further information please read the readme.md')

    try:
        genome_path = Path(args.genome).resolve()
        if(not os.access(str(genome_path), os.R_OK)):
            log.error('genome file not readable! path=%s', genome_path)
            sys.exit(f'ERROR: genome file ({genome_path}) not readable!')
        if(genome_path.stat().st_size == 0):
            log.error('empty genome file! path=%s', genome_path)
            sys.exit(f'ERROR: genome file ({genome_path}) is empty!')
    except:
        log.error('provided genome file not valid! path=%s', args.genome)
        sys.exit(f'ERROR: genome file ({args.genome}) not valid!')
    log.info('genome-path=%s', genome_path)

    tmp_path = Path(tempfile.mkdtemp())
    log.info('tmp-path=%s', tmp_path)
    log.info('output-path=%s', output_path)

    # workflow configurations
    global mode, characterize
    mode = args.mode
    log.info('mode=%s', mode)
    characterize = args.characterize
    log.info('characterize=%s', characterize)
