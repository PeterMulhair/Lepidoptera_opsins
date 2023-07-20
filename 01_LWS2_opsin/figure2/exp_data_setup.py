import glob
import argparse
from Bio import SeqIO
from subprocess import call as unix

'''
Code to take bam coverage file and
genome fasta file to set up files
to plot coverage of certain genes 
with gene tracks using plot_tracks_exp.R
'''

parse = argparse.ArgumentParser()

parse.add_argument("-g", "--genome",type=str, help="path to genomes file input",required=True)
parse.add_argument("-b", "--bam",type=str, help="path to bam file input",required=True)
parse.add_argument("-n", "--name",type=str, help="name of sample to process e.g. OrgAnti_1instar or OrgAnti_male_head",required=True)

args = parse.parse_args()

#For a given genome create index file
unix('cp ' + args.genome + ' .', shell=True)
genome_fas = args.genome.split('/')[-1]
unix('samtools faidx ' + genome_fas, shell=True)
unix('mv ' + genome_fas + '.fai chrom.sizes', shell=True)
unix('rm ' + genome_fas, shell=True)
unix('cut -f1,2 chrom.sizes >> chrom1.sizes', shell=True)
unix('rm chrom.sizes', shell=True)
unix('mv chrom1.sizes chrom.sizes', shell=True)

#Make dictionary of chromosome ID to chromosome number
ID_chrm = {}
with open(args.genome) as f:
    for record in SeqIO.parse(f, 'fasta'):
        header = record.description
        ID = record.id
        if 'chromosome' in header:
            chrm = header.split(': ')[-1].strip('\n')
            ID_chrm[ID] = chrm


with open('chrom.sizes') as f, open('chrom1.sizes', 'w') as outF:
    for line in f:
        lines = line.split('\t')
        chrm_no = lines[0]
        info1 = lines[1].strip('\n')
        if chrm_no in ID_chrm:
            chrmID = ID_chrm[chrm_no]
        else:
            chrmID = chrm_no

        outF.write(chrmID + '\t' + info1 + '\n')

unix('rm chrom.sizes', shell=True)
unix('mv chrom1.sizes chrom.sizes', shell=True)

#Use bedtools to make tsv file of RNA coverage across genome
unix('bedtools genomecov -ibam ' + args.bam + ' -bga >> ' + args.name + '_genome_cov.tsv', shell=True)

with open(args.name + '_genome_cov.tsv') as f, open(args.name + '_genome_cov_chrm.tsv', 'w') as outF:
    for line in f:
        lines = line.split('\t')
        chrm_no = lines[0]
        info1 = lines[1]
        info2 = lines[2]
        info3 = lines[3].strip('\n') 
        if chrm_no in ID_chrm:
            chrmID = ID_chrm[chrm_no]
        else:
            chrmID = chrm_no

        outF.write(chrmID + '\t' + info1 + '\t' + info2 + '\t' + info3 + '\n')

#Create sorted bigwig files to use as input to plotting code
unix('sort -k1,1 -k2,2n ' + args.name + '_genome_cov_chrm.tsv >> ' + args.name + '_genome_cov_chrm_sorted.tsv', shell=True)
unix('bedGraphToBigWig ' + args.name + '_genome_cov_chrm_sorted.tsv chrom.sizes ' + args.name + '_genome_cov_chrm_sorted.bw', shell=True)
