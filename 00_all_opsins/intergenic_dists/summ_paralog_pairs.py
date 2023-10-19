#!/usr/bin/env python3
import os
import glob
from Bio import SeqIO
import statistics
from collections import Counter

'''
Code to parse blast output and measure
intergenic distances between paralog
gene pairs
'''

#Usage: python summ_paralog_pairs.py

os.chdir('blastout')

spshort_fa = {}#Create dictionary of shortened species name to proteome fasta file (assumes proteome fasta file looks like Aricia_agestis-GCA_905147365.1-2021_11-pep.fa)
for fa in glob.glob('../*fa'):
    sp = fa.split('/')[-1].split('-')[0]
    spshort = sp.split('_')[0][:3] + sp.split('_')[1][:4].title()
    spshort_fa[sp] = fa

mean_para_list = []
median_para_list = []
for blastout in glob.glob('*_blastoutput.tsv'):#Parse every species blastoutput file
    paralog_count = 0
    paralog_pair_dists = []
    sp = blastout.split('_blastout')[0]
    print(sp)
    raw = spshort_fa[sp]
    gene_loc_mid = {}
    gene_chrm = {}
    with open(raw) as f:
        for record in SeqIO.parse(f, 'fasta'):
            header = record.description
            geneID = header.split(' ')[3].split(':')[1]
            chrm = header.split(' ')[2].split(':')[1]
            gene_loc = (int(header.split(' ')[2].split(':')[2]) + int(header.split(' ')[2].split(':')[3]))/2
            gene_loc_mid[geneID] = gene_loc
            gene_chrm[geneID] = chrm

    #Create dictionary of all (non self hit) paralog gene pairs 
    paralog_dict = {}
    with open(blastout) as f:
        for line in f:
            lines = line.split('\t')
            query = lines[0]
            subject = lines[1]
            pident = lines[3]
            if (query != subject):#Remove self hits
                paralog_dict[query] = subject

    with open(blastout) as f:
        for line in f:
            line = line.strip('\n')
            lines = line.split('\t')
            query = lines[0]
            subject = lines[1]
            pident = lines[3]

            if subject in paralog_dict:
                subject_reciphit = paralog_dict[subject]

                query_start = lines[5]
                query_end = lines[6]
                query_hit_len = int(query_end) - int(query_start)
                query_len = lines[7]
                subject_start = lines[8]
                subject_end = lines[9]
                subject_hit_len = int(subject_end) - int(subject_start)
                subject_len = lines[10]
                query_cov = int((int(query_hit_len)/int(query_len))*100)
                subject_cov = int((int(subject_hit_len)/int(subject_len))*100)

                if (query != subject) and (float(pident) >= 65) and (subject_reciphit == query) and (query_cov >= 70) and (subject_cov >= 70):#If blast hit is >=65% identity and the hit reciprocally hits the query as its top hit and hit coverage is >= 70% for query and subject
                    query_loc = gene_loc_mid[query]
                    query_chrm = gene_chrm[query]
                    subject_loc = gene_loc_mid[subject]
                    subject_chrm = gene_chrm[subject]
                    paralog_count+=1
                    if query_chrm == subject_chrm:#If gene pairs are located on same chromosome
                        if query_loc > subject_loc:
                            paralog_dist = query_loc - subject_loc#Get intergenic paralog distances 
                        else:
                            paralog_dist = subject_loc - query_loc
                        paralog_pair_dists.append(paralog_dist)

    #Prints intergenic statistics per species
    print(int(statistics.mean(paralog_pair_dists)), 'mean dist between paralog pairs')
    print(int(statistics.median(paralog_pair_dists)), 'median dist between paralog pairs')
    print(paralog_count)
    print('\n')
    mean_para_dist = statistics.mean(paralog_pair_dists)
    mean_para_list+=paralog_pair_dists

#Prints overall intergenic statistics    
print('\n')
print(int(statistics.mean(mean_para_list)))
print(int(statistics.median(mean_para_list)))
