
import argparse
from pathlib import Path

from Bio import SeqIO
from scipy import stats


def calc_rds_pval(abs_diff_mean, no_plasmids, plasmid_hits, no_chromosomes, chromosome_hits):
    p_hit_frequency = plasmid_hits / no_plasmids
    c_hit_frequency = chromosome_hits / no_chromosomes
    hit_frequency_ratio = ((p_hit_frequency / (p_hit_frequency + c_hit_frequency)) - 0.5) * 2
    odds_ratio, p_value = stats.fisher_exact([[plasmid_hits, chromosome_hits], [no_plasmids, no_chromosomes]], alternative='two-sided')
    abs_diff_norm = abs(p_hit_frequency - c_hit_frequency) / abs_diff_mean
    rds = hit_frequency_ratio * (1 - p_value) * abs_diff_norm
    return rds, p_value


parser = argparse.ArgumentParser(
    description='Calculate RDS and extract found MPS.'
)
parser.add_argument('--plasmids', action='store', help='Path to plasmid fasta file.')
parser.add_argument('--chromosomes', action='store', help='Path to chromosome fasta file.')
parser.add_argument('--counts', action='store', help='Path to protein hit counts tsv file.')
parser.add_argument('--cluster-seqs', action='store', dest='cluster_seqs', help='Path to UniRef90 fasta file.')
parser.add_argument('--cluster-info', action='store', dest='cluster_info', help='Path to UniRef90 tsv file.')
args = parser.parse_args()

counts_path = Path(args.counts).resolve()
plasmids_path = Path(args.plasmids).resolve()
chromosomes_path = Path(args.chromosomes).resolve()
cluster_seqs_path = Path(args.cluster_seqs).resolve()
cluster_info_path = Path(args.cluster_info).resolve()

# count plasmid records
print('count plasmids...')
no_plasmids = 0
with plasmids_path.open() as fh:
    for record in SeqIO.parse(fh, 'fasta'):
        no_plasmids += 1
print(f'\tplasmids found: {no_plasmids}')

# count chromosome records
print('count chromosomes...')
no_chromosomes = 0
with chromosomes_path.open() as fh:
    for record in SeqIO.parse(fh, 'fasta'):
        no_chromosomes += 1
print(f'\tchromosomes found: {no_chromosomes}')

# store UniRef90 cluster info
print('store UniRef90 cluster info...')
clusters = {}
with cluster_info_path.open() as fh:
    for line in fh:
        (id, product, length) = line.split('\t')
        clusters[id] = {
            'id': id,
            'product': product,
            'length': int(length)
        }
print('\t...done')

# compute mean absolute hit frequency differences
print('compute mean absolute hit frequency differences...')
no_mps = 0
abs_hit_freq_diff_sum = 0.
with counts_path.open() as fh:
    for line in fh:
        (id, plasmid_hits, chromosome_hits) = line.split('\t')
        plasmid_hits = int(plasmid_hits)
        chromosome_hits = int(chromosome_hits)
        if(plasmid_hits > 0 or chromosome_hits > 0):
            abs_hit_freq_diff_sum += abs((plasmid_hits / no_plasmids) - (chromosome_hits / no_chromosomes))
            no_mps += 1
mahfd = abs_hit_freq_diff_sum / no_mps
print(f'\tmean abs hit freq diff: {mahfd}')

# calculate RDS for all MPS
print('calculate RDS for all MPS...')
i = 0
i_hit = 0
mpss = {}
with counts_path.open() as fh:
    for line in fh:
        i += 1
        (id, plasmid_hits, chromosome_hits) = line.split('\t')
        plasmid_hits = int(plasmid_hits)
        chromosome_hits = int(chromosome_hits)
        if(plasmid_hits > 0 or chromosome_hits > 0):
            i_hit += 1
            mps = clusters[id]
            rds, p_value = calc_rds_pval(mahfd, no_plasmids, plasmid_hits, no_chromosomes, chromosome_hits)
            mps['rds'] = rds
            mps['p_value'] = p_value
            mpss[id] = mps
        if((i % 1000000) == 0):
            print(f'\t... {i} processed')
clusters.clear()

# write MPS fasta and tsv files
print('write filtered MPS fasta and tsv files...')
mpss = {k: v for k, v in mpss.items() if v['rds'] != 0.0}
mps_fasta_path = Path('mps.faa')
mps_tsv_path = Path('mps.tsv')
mps_full_tsv_path = Path('mps.raw.tsv')
with cluster_seqs_path.open() as cluster_in_fh, mps_fasta_path.open('wt') as mps_fasta_fh, mps_tsv_path.open('wt') as mps_tsv_fh, mps_full_tsv_path.open('wt') as mps_full_tsv_fh:
    for record in SeqIO.parse(cluster_in_fh, 'fasta'):
        id = record.id
        if(id in mpss):
            mps = mpss[id]
            mps_fasta_fh.write(f'>{id}\n{record.seq}\n')
            mps_tsv_fh.write(f"{id}\t{mps['product']}\t{mps['length']}\t{mps['rds']}\n" )
            mps_full_tsv_fh.write(f"{id}\t{mps['product']}\t{mps['length']}\t{mps['p_value']}\t{mps['rds']}\n")
print('\t...done')
