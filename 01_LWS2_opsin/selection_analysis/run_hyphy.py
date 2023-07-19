#!/usr/bin/env python3
import os
import sys
import json
import glob
import shutil
import argparse
import subprocess
from Bio import SeqIO
from datetime import datetime
from joblib import Parallel, delayed
from subprocess import call as unix
from collections import Counter, defaultdict

# Author: Peter Mulhair
# Date: 15/11/2021
# Usage: python3 run_hyphy.py --gene LWS

'''
Pipeline to run selective pressure analysis
using Hyphy. Uses FitMG94 to assess rate of 
selection on genes - pulls out dN, dS, and 
dN/dS for each gene.
'''

if shutil.which('mafft') is None:
    sys.exit('ERROR: mafft not installed locally')
elif shutil.which('pal2nal.pl') is None:
    sys.exit('ERROR: pal2nal not installed locally')
elif shutil.which('vespasian') is None:
    sys.exit('ERROR: vespasian not installed locally')
elif shutil.which('hyphy') is None:
    sys.exit('ERROR: hyphy not installed locally')

parse = argparse.ArgumentParser()

parse.add_argument("--gene", type=str, help="name of opsin gene to analyse eg. LWS, Blue, UV, Copsin, Rh7", required = True)
parse.add_argument("--genomes", type=str, help="path to genome fasta files", required = True)

args = parse.parse_args()

date = datetime.today().strftime('%Y_%m_%d')

#Get species IDs that are in Lepidoptera
lep_sp = []
for fasta in glob.glob(args.genomes + '*fasta'):
    species = fasta.split('/')[-1].split('_')[2].split('.')[0][2:-1]
    order = fasta.split('/')[-1].split('_')[2].split('.')[0][:2]
    if order == 'il':
        lep_sp.append(species)


#Find genes which are present in more than 1 copy
sp_IDs = defaultdict(list)
sp_copy = []
multi_copy_gene = defaultdict(list)
for opsin_fa in glob.glob('data/*_opsins.fas'):
    with open(opsin_fa) as f:
        for record in SeqIO.parse(f, 'fasta'):
            ID = record.id
            sp = ID.split('|')[0]
            if (ID not in sp_IDs[sp]) and (args.gene in ID) and (sp in lep_sp):
                sp_IDs[sp].append(ID)
                sp_copy.append(sp)
                
sp_count = Counter(sp_copy)
for species, num in sp_count.items():
    if num != 1:
        genes = sp_IDs[species]
        multi_copy_gene[species] = genes

#Make sure there are genes in multicopy in gene family
if len(multi_copy_gene) == 0:
    sys.exit('No multi copy genes present in ' + args.gene)

lep_nuc_opsin = {}
with open('data/' + args.gene + '_nucl_lepi.fas') as f:
    for record in SeqIO.parse(f, 'fasta'):
        ID = record.id
        seq = str(record.seq)
        lep_nuc_opsin[ID] = seq

lep_SC_prot_opsin = {}
lep_SC_nucl_opsin = {}
for species, num in sp_count.items():
    if num == 1:
        with open('data/' + species + '_opsins.fas') as f:
            for record in SeqIO.parse(f, 'fasta'):
                ID = record.id
                if args.gene in ID:
                    seq = str(record.seq)
                    lep_SC_prot_opsin[species] = seq
        
        with open('data/' + args.gene + '_nucl_lepi.fas') as f:
            for record in SeqIO.parse(f, 'fasta'):
                ID = record.id
                if species in ID:
                    nuc = str(record.seq)
                    lep_SC_nucl_opsin[species] = nuc
                    
#Get intron counts for genes
gene_intron_count = {}
with open('data/opsin_intron_count.tsv') as f:
    for line in f:
        lines = line.split('\t')
        sp_gene = lines[0]
        intron_count = lines[1].strip()
        gene_intron_count[sp_gene] = intron_count
        

def hyphy(species):
    os.makedirs('results', exist_ok=True)
    os.makedirs('results/FitMG94/' + args.gene, exist_ok=True)
    os.mkdir('results/FitMG94/' + args.gene + '/' + species)
    os.mkdir('results/FitMG94/' + args.gene + '/' + species + '/prot_align')
    os.mkdir('results/FitMG94/' + args.gene + '/' + species + '/codon_align')
    
    gene_list = multi_copy_gene[species]

    #Get opsin protein sequences for species in question
    sp_prot = {}
    with open('data/' + species + '_opsins.fas') as f:
        for record in SeqIO.parse(f, 'fasta'):
            ID = record.id
            if ID in gene_list:
                seq = str(record.seq)
                sp_prot[ID] = seq
    
    for gene, seqs in sp_prot.items():
        geneID = gene.split('|')[1]
        intron_count = gene_intron_count[species  + '_' + geneID]
        
        with open('results/FitMG94/' + args.gene + '/' + species + '/prot_align/' + species + '_lep_' + geneID + '.fas', 'w') as outF:
            for lep_sp, lep_gene in lep_SC_prot_opsin.items():
                outF.write('>' + lep_sp + '\n' + lep_gene + '\n')
            outF.write('>' + gene + '\n' +  seqs + '\n')

        unix('mafft --quiet results/FitMG94/' + args.gene + '/' + species + '/prot_align/' + species + '_lep_' + geneID + '.fas > results/FitMG94/' + args.gene + '/' + species + '/prot_align/' + species + '_lep_' + geneID + '_aln.fas', shell=True)
            
        nuc_seq_data = lep_nuc_opsin[gene]
        with open('results/FitMG94/' + args.gene + '/' + species + '/' + species + '_lep_' + geneID + '_nuc.fas', 'w') as outF1:
            for lep_sp_nuc, lep_gene_nuc in lep_SC_nucl_opsin.items():
                outF1.write('>' + lep_sp_nuc + '\n' + lep_gene_nuc + '\n')
            outF1.write('>' + gene + '\n' +  nuc_seq_data + '\n')

        unix('pal2nal.pl results/FitMG94/' + args.gene + '/' + species + '/prot_align/' + species + '_lep_' + geneID + '_aln.fas results/FitMG94/' + args.gene + '/' + species + '/' + species + '_lep_' + geneID + '_nuc.fas >> results/FitMG94/' + args.gene + '/' + species + '/codon_align/' + species + '_lep_' + geneID + '_aln.fasta -output fasta', shell=True)

        ##Get gene trees with vespasian
        subprocess.run("bash -c 'source activate vespasian && vespasian infer-gene-trees --warnings --progress results/FitMG94/" + args.gene + "/" + species + "/codon_align/ raw/lep_tree.nwk -o results/FitMG94/" + args.gene + "/" + species + "/gene-trees/ && source deactivate'", shell=True)
        
        ##Run hyphy with FitMG94 model
        subprocess.run("bash -c 'source activate hyphy && hyphy ~/bin/hyphy-analyses/FitMG94/FitMG94.bf --alignment results/FitMG94/" + args.gene + "/" + species + "/codon_align/" + species + "_lep_" + geneID + "_aln.fasta --tree results/FitMG94/" + args.gene + "/" + species + "/gene-trees/" + species + "_lep_" + geneID + "_aln.nwk --type local && source deactivate'", shell=True)

        ##Parse hyphy results
        #Check to make sure hyphy output is not empty
        if os.stat('results/FitMG94/' + args.gene + '/' + species + '/codon_align/' + species + '_lep_' + geneID + '_aln.fasta.FITTER.json').st_size == 0:
            with open('fail_hyphy.txt', 'a+') as outF3:
                outF3.write(geneID + '\n')
        else:
            with open('results/FitMG94/' + args.gene + '/' + species + '/codon_align/' + species + '_lep_' + geneID + '_aln.fasta.FITTER.json') as f, open('results/FitMG94/' + args.gene + '/' + species + '/' + species + '_' + geneID[:-1] + '_output.tsv', 'a+') as outF2:
                res_dict = json.load(f)
                branch_results = res_dict['branch attributes']
                branch_results = branch_results['0']
                sp_res = branch_results[species + '_' + geneID]
                CI = sp_res['Confidence Intervals']
                omega = CI['MLE']
                dN = sp_res['dN']
                dS = sp_res['dS']
                
                outF2.write(species + '_' + geneID + '\t' + str(omega) + '\t' + str(dN) + '\t' + str(dS) + '\t' + intron_count + '\n')

Parallel(n_jobs=40)(delayed(hyphy)(sp) for sp in multi_copy_gene.keys())
