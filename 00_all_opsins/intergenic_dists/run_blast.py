#!/usr/bin/env python3
import os
import glob
import argparse
from subprocess import call as unix
from joblib import Parallel, delayed

'''
Code used to run blastp search for each species
proteome against itself. Assumes you have species
proteomes downloaded in the pwd and have creeated
blast databases for each one (and placed in blastdb/)
'''

#Usage: python run_blast.py --threads N

parse = argparse.ArgumentParser()

parse.add_argument("--threads",type=str, help="number of BLAST searches to run in parallel",required=True)

args = parse.parse_args()

threads = int(args.threads)

os.makedirs('blastout', exist_ok=True

def run_blast(fa):
    species = fa.split('-')[0]
    print(species)
    unix('blastp -query ' + fa + ' -db blastdb/' + species + ' -evalue 1e-10 -max_target_seqs 2 -num_threads 4 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out blastout/' + species + '_blastoutput.tsv', shell=True)

fa_list = []
for fa in glob.glob('*fa'):
    fa_list.append(fa)

Parallel(n_jobs=threads)(delayed(run_blast)(fas) for fas in fa_list)    
