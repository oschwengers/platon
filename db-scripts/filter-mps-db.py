
import argparse
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description='Calculate RDS and extract found MPS.'
)
parser.add_argument('--mps-raw', action='store', dest='mps_raw', help='Path to raw mps file.')
parser.add_argument('--cluster-seqs', action='store', dest='cluster_seqs', help='Path to UniRef90 fasta file.')
args = parser.parse_args()

mps_raw_path = Path(args.mps_raw).resolve()
cluster_seqs_path = Path(args.cluster_seqs).resolve()

mpss = {}
with mps_raw_path.open() as fh:
    for line in fh:
        (id, product, length, rds) = line.split('\t')
        length = int(length)
        rds = float(rds)
        if(rds < 0.0 or rds > 0.0):
            mpss[id] = {
                'id': id,
                'product': product,
                'length': length,
                'rds': rds
            }

print('write filtered MPS fasta and tsv files...')
mps_fasta_path = Path('mps.faa')
mps_tsv_path = Path('mps.tsv')
with cluster_seqs_path.open() as cluster_in_fh, mps_fasta_path.open('wt') as mps_fasta_fh, mps_tsv_path.open('wt') as mps_tsv_fh:
    for record in SeqIO.parse(cluster_in_fh, 'fasta'):
        id = record.id
        if(id in mpss):
            mps = mpss[id]
            mps_fasta_fh.write(f'>{id}\n{record.seq}\n')
            mps_tsv_fh.write(f"{id}\t{mps['product']}\t{mps['length']}\t{mps['rds']}\n")
print('\t...done')
