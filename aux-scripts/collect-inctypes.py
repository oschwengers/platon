import argparse
import os
import json
import sys

from pathlib import Path
import platon.constants as pc

parser = argparse.ArgumentParser(
    prog='collect-inctypes',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='Collect and print detected inc types as a present/absence matrix.',
    epilog=f'Citation:\n{pc.CITATION}\n\nGitHub:\nhttps://github.com/oschwengers/platon'
)
parser.add_argument('results', metavar='<results>', nargs='+', help='Platon result JSON files')
args = parser.parse_args()

global_inc_types = set()
sample_inc_types = {}
for result_file in args.results:
    try:
        result_path = Path(result_file).resolve()
        if(not os.access(str(result_path), os.R_OK)):
            sys.exit(f'ERROR: result file ({result_path}) not readable!')
        sample = result_path.stem

        inc_types = []
        with result_path.open() as fh:
            content = json.load(fh)
        for contig_id, value in content.items():
            # for inc_type in value.get('inc_types', []):
            for inc_type in value.get('plasmid_hits', []):
                
                # [{"type":"IncFIB(K)_1_JN233704","start":1,"end":559,"strand":"+","identity":0.98929,"coverage":1.0,"bitscore":1002}]
                # inc_type = inc_type['type']

                # [{"contig_start":1,"contig_end":1475,"plasmid_start":263545,"plasmid_end":265019,"plasmid":{"id":"NC_012556.1","length":324503},"coverage":1.0,"identity":1.0}]
                inc_type = inc_type['plasmid']['id']
                inc_types.append(inc_type)
                global_inc_types.add(inc_type)
        sample_inc_types[sample] = inc_types
    except:
        print(f'could not read file: {result_file}')

global_inc_types = sorted(list(global_inc_types))
print(f"Sample\t{'\t'.join(global_inc_types)}")
for sample, inc_types in sample_inc_types.items():
    inc_type_presence = ['1' if i in inc_types else '0' for i in global_inc_types]
    print(f"{sample}\t{'\t'.join(inc_type_presence)}")
