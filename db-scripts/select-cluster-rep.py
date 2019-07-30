#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument( '--clusters', '-c', help='NCBI PLCA NRP cluster proteins' )
parser.add_argument( '--proteins', '-p', help='NCBI NRP proteins' )
args = parser.parse_args()


nrpc = set()
nrpc_proteins = {}
with open(os.path.realpath(args.clusters), 'r') as fh:
    for line in fh:
        if(line[0] is not '#'):
            cols = line.split('\t')
            nrpc.add(cols[0])
            protein = {
                'cluster': cols[0],
                'wp': cols[1],
                'length': int(cols[5]),
            }
            if(cols[0] in nrpc_proteins):
                prot_list = nrpc_proteins[cols[0]]
                prot_list.append(protein)
            else:
                nrpc_proteins[cols[0]] = [protein]

nrp = set()
for record in SeqIO.parse(os.path.realpath(args.proteins), 'fasta'):
    nrp.add(str(record.id))


for cluster in nrpc:
    prot_list = nrpc_proteins[cluster]
    valid_proteins = []
    for protein in prot_list:
        if(protein['wp'] in nrp):
            valid_proteins.append(protein)
    if(len(valid_proteins) == 1 or len(valid_proteins) == 2):
        print(valid_proteins[0]['wp'])
    elif(len(valid_proteins) > 2):
        valid_proteins = sorted(valid_proteins, key=lambda k: k['length'])
        protein = valid_proteins[int(len(valid_proteins) / 2)]
        print(protein['wp'])
    else:
        #print( cluster )
        pass