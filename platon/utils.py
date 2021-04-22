import argparse
import collections
import logging
import multiprocessing as mp
import os
import re
import sys
import subprocess as sp

import platon
import platon.constants as pc


log = logging.getLogger('UTILS')
Version = collections.namedtuple('Version', ['major', 'minor', 'patch'], defaults=[0, 0]) # named tuple for version checking, defaults are zero for missing minor/patch

def print_version(self):
    return f'v{self.major}.{self.minor}.{self.patch}'

Version.__str__ = print_version

VERSION_MIN_DIGIT = -1
VERSION_MAX_DIGIT = 1000000000000
VERSION_REGEX = re.compile(r'(\d+)\.(\d+)(?:\.(\d+))?')  # regex to search for version number in tool output. Takes missing patch version into consideration.
DEPENDENCIES = [  # List of dependencies: tuples for: min version, max version, tool name & command line parameter, dependency check exclusion options
    (Version(2,6,3), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('prodigal', '-v')),
    (Version(2,0,4), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('diamond', 'help')),
    (Version(2,10,1), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('blastn', '-version')),
    (Version(3,3,1), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('hmmsearch', '-h')),
    (Version(4,0,0), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('nucmer', '-V')),
    (Version(1,1,2), Version(VERSION_MAX_DIGIT, VERSION_MAX_DIGIT, VERSION_MAX_DIGIT), VERSION_REGEX, ('cmscan', '-h'))
]


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='platon',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Identification and characterization of bacterial plasmid contigs from short-read draft assemblies.',
        epilog=f'Citation:\n{pc.CITATION}\n\nGitHub:\nhttps://github.com/oschwengers/platon',
        add_help=False
    )
    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('genome', metavar='<genome>', help='draft genome in fasta format')
    arg_group_io.add_argument('--db', '-d', action='store', help='database path (default = <platon_path>/db)')
    arg_group_io.add_argument('--prefix', '-p', action='store', default=None, help='Prefix for output files')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
    
    arg_group_workflow = parser.add_argument_group('Workflow')
    arg_group_workflow.add_argument('--mode', '-m', action='store', type=str, choices=['sensitivity', 'accuracy', 'specificity'], default='accuracy', help='applied filter mode: sensitivity: RDS only (>= 95%% sensitivity); specificity: RDS only (>=99.9%% specificity); accuracy: RDS & characterization heuristics (highest accuracy) (default = accuracy)')
    arg_group_workflow.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize all contigs')
    
    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)')
    arg_group_general.add_argument('--version', action='version', version=f'%(prog)s {platon.__version__}')
    return parser.parse_args()


def read_tool_output(dependency):
        """Method for reading tool version with regex. Input: regex expression, tool command. Retursn: version number."""
        version_regex = dependency[2]
        command = dependency[3]
        try:
            tool_output = str(sp.check_output(command, stderr=sp.STDOUT)) # stderr must be added in case the tool output is not piped into stdout
        except FileNotFoundError:
            log.exception('dependency not found! tool=%s', command[0])
            sys.exit(f'ERROR: {command[0]} not found or not executable! Please make sure {command[0]} is installed and executable.')
        except sp.CalledProcessError:
            log.exception('dependency check failed! tool=%s', command[0])
            sys.exit(f'ERROR: {command[0]} could not be executed! Please make sure {command[0]} is installed and executable.')
        version_match = re.search(version_regex, tool_output)
        
        try:
            if version_match is None:
                log.error('no dependency version detected! no regex hit in dependency output: regex=%s, command=%s', version_regex, command)
                sys.exit(f'ERROR: Could not detect/read {command[0]} version!')

            major = version_match.group(1)
            minor = version_match.group(2)
            patch = version_match.group(3)
            if major is None:
                log.error('no dependency version detected! no regex hit in dependency output: regex=%s, command=%s', version_regex, command)
                sys.exit(f'ERROR: Could not detect/read {command[0]} version!')
            elif minor is None:
                version_output = Version(int(major))
            elif patch is None:
                version_output = Version(int(major), int(minor))
            else:
                version_output = Version(int(major), int(minor), int(patch))
            return version_output
        except:
            log.error('no dependency version detected! no regex hit in dependency output: regex=%s, command=%s', version_regex, command)
            sys.exit(f'ERROR: Could not detect/read {command[0]} version!')


def check_version(tool, min, max):
    """Method for checking tool versions with required version. Input: tool version, minimum and maximum version. Returns: boolean value for positive or negative check."""
    if tool.major < min.major or tool.major > max.major:
        return False
    else:
        if tool.major == min.major or tool.major == max.major:
            if tool.minor < min.minor and tool.major == min.major:
                return False
            elif tool.minor > max.minor and tool.major == max.major:
                return False
            else:
                if tool.minor == min.minor or tool.minor == max.minor:
                    if tool.patch < min.patch and tool.minor == min.minor:
                        return False
                    elif tool.patch > max.patch and tool.minor == max.minor and tool.major == max.major:
                        return False
                    else:
                        return True
                else:
                    return True
        else:
            return True


def test_dependencies():
    """Test the proper installation of necessary 3rd party executables."""
    for dependency in DEPENDENCIES:
        version = read_tool_output(dependency)
        check_result = check_version(version, dependency[0], dependency[1])
        if (check_result == False):
            log.error('wrong dependency version for %s: installed=%s, minimum=%s', dependency[3][0], version, dependency[0])
            sys.exit(f'ERROR: Wrong {dependency[3][0]} version installed. Please, install {dependency[3][0]} version {dependency[0]}!')
        else:
            log.info('dependency check: tool=%s, version=%s', dependency[3][0], version)


