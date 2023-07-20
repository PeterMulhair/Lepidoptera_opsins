import os
import json
import glob
import argparse
from Bio import SeqIO
from collections import defaultdict
from subprocess import call as unix

parse = argparse.ArgumentParser()

parse.add_argument("--gtf",type=str, help="path to full gtf files per species to parse",required=True)
parse.add_argument("--gff",type=str, help="path to opsin gff files i.e. ~/Lepidoptera_opsins/00_all_opsins/gff/ in this repo",required=True)

args = parse.parse_args()

os.makedirs('LWS1_intron', exist_ok=True)
os.makedirs('LWS2_intronless', exist_ok=True)

noctuoidea_sp = []
noctuoidea_sp_short = []
with open('Noctuoidea_sp.txt') as f:
    for line in f:
        sp = line.strip()
        sp_short = sp.split('_')[0][:3] +  sp.split('_')[1][:4].title()
        noctuoidea_sp.append(sp)
        noctuoidea_sp_short.append(sp_short)
            
gtf_dict = {}
for gtf in glob.glob(args.gtf + '*.gtf'):
    sp_name = gtf.split('/')[-1].split('-')[0]
    sp_short = sp_name.split('_')[0][:3] +  sp_name.split('_')[1][:4].title() 
    gtf_dict[sp_short] = gtf
    

name_change_dict = {'PieBras': 'PieBrab', 'HerTars': 'HerTari', 'InaIo': 'AglIoxx', 'ParAege': 'ParAegt', 'YpoSede': 'YpoSedl', 'NemSwam': 'NemSwae', 'EnnQuer': 'EnnQuei', 'SpiLute': 'SpiLutu', 'XesC-Ni': 'XesCnig', 'NymUrti': 'AglUrti', 'AphHype': 'AphHyp', 'HemFuci': 'HemFuc', 'PtiCapu': 'PtiCapc', 'EilDepr': 'EilDepe'}

with open('species_chromosome_dict.json') as f:
    sp_chrm_dict = json.load(f)

for sp_short in noctuoidea_sp_short:
    gff = glob.glob(args.gff + '*' + sp_short + '*gff')
    gff = gff[0]
    if sp_short in gtf_dict:
        intronLWS_count = 0
        intronlessLWS_count = 0
        chromosome_dict = sp_chrm_dict[sp_short]
        gtf_output = gtf_dict[sp_short]
        with open(gff) as f:
            lines = f.read()
            hits = lines.split('# --- START OF GFF DUMP ---')
            for hit in hits:
                gene_intron_lens = defaultdict(list)
                intron_lens = []
                exon_lens = []
                exon_count = 0
                hit_info = hit.split('\n')
                for info in hit_info:
                    if info.startswith('#'):
                        continue
                    else:
                        gene_info = info.split('\t')
                        if len(gene_info) > 1:
                            if (gene_info[2] == 'gene') and ('|LWS' in gene_info[-1]):
                                if 'intron' in hit:
                                    intronLWS_count+=1
                                    chrm = gene_info[0]
                                    chrm_num = chromosome_dict[chrm]
                                    print(sp_short, chrm_num, intronLWS_count)
                                    gene_start = gene_info[3]
                                    gene_end = gene_info[4]
                                    gene_centre = (int(gene_start) + int(gene_end))/2
                                    gene_centre = int(round(gene_centre))
                                    gene_region_start = gene_centre - 50000
                                    gene_region_end = gene_centre + 50000
                                    with open('LWS1_intron/' + sp_short + '_' + str(intronLWS_count)  + '_gtf_tracks.ini', 'w') as outF:
                                        outF.write('[genes 1]\nfile = /home/zoo/zool2500/tree_of_life/raw/updated_DToL_data/opsin_genes/expression_analysis/plot_tracks/GenomeTracks/' + gtf_output + '\nheight = 1.4\nstyle = tssarrow\nmax_labels = 1\n\n[spacer]\nheight = 3\n\n[x-axis]\nwhere = bottom\nfontsize = 10\n')
                                    with open('genome_plots.sh', 'a+') as outF1:
                                        outF1.write('pyGenomeTracks --tracks LWS1_intron/' + sp_short + '_' + str(intronLWS_count)  + '_gtf_tracks.ini  --region ' + str(chrm_num) + ':' + str(gene_region_start) + '-' + str(gene_region_end) + ' --trackLabelFraction 0.2 --width 25 --dpi 130 -o LWS1_intron/' + sp_short + '_' + str(intronLWS_count) + '_LWS1.pdf\n')
                                else:
                                    intronlessLWS_count+=1
                                    chrm = gene_info[0]
                                    chrm_num = chromosome_dict[chrm]
                                    gene_start = gene_info[3]
                                    gene_end = gene_info[4]
                                    gene_centre = (int(gene_start) + int(gene_end))/2
                                    gene_centre = int(round(gene_centre))
                                    gene_region_start = gene_centre - 50000
                                    gene_region_end = gene_centre + 50000
                                    with open('LWS2_intronless/' + sp_short + '_' + str(intronlessLWS_count)  + '_gtf_tracks.ini', 'w') as outF:
                                        outF.write('[genes 1]\nfile = /home/zoo/zool2500/tree_of_life/raw/updated_DToL_data/opsin_genes/expression_analysis/plot_tracks/GenomeTracks/' + gtf_output + '\nheight = 4\nstyle = tssarrow\nmax_labels = 1\n\n[spacer]\nheight = 0.5\n\n[x-axis]\nwhere = bottom\nfontsize = 10\n')

                                    with open('genome_plots.sh', 'a+') as outF1:
                                        outF1.write('pyGenomeTracks --tracks LWS2_intronless/' + sp_short + '_' + str(intronlessLWS_count)  + '_gtf_tracks.ini  --region ' + str(chrm_num) + ':' + str(gene_region_start) + '-' + str(gene_region_end) + ' --trackLabelFraction 0.2 --width 25 --dpi 130 -o LWS2_intronless/' + sp_short + '_' + str(intronlessLWS_count) + '_LWS2.pdf\n')
