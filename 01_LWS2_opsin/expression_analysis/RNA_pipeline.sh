##Trim and map RNA reads to genome
sudo java -jar trimmomatic-0.39.jar PE -threads 20 -phred33 V4_1.fq V4_2.fq V4_1_paired_trim.fastq V4_1_unpaired_trim.fastq V4_2_paired_trim.fastq V4_2_unpaired_trim.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:10:30 MINLEN:35

sudo bowtie2-build --threads 40 -f GCA_916999025.1_ilOrgAnti1.1_genomic.fasta OrgAnti

bowtie2 --threads 20 -x OrgAnti -1 V4_1_paired_trim.fastq -2 V4_2_paired_trim.fastq | samtools view -S -h -F4 - > OrgAnti_mapped.sam

#Sort mapped read output
samtools sort OrgAnti_mapped.sam -o OrgAnti_sorted_instar.bam
samtools index OrgAnti_sorted_instar.bam

#Reformat fastq headers if required (necessary to run Trinity)
awk '{{print (NR%4 == 1) ? "@1_" ++i "/1": $0}}' SRR11196268_pass_1_paired_trim.fastq > SRR11196268_pass_1_paired_trim_reform.fastq
awk '{{print (NR%4 == 1) ? "@1_" ++i "/2": $0}}' SRR11196268_pass_2_paired_trim.fastq > SRR11196268_pass_2_paired_trim_reform.fastq


##Make transcriptome
Trinity --seqType fq --max_memory 50G --left V4_1_paired_trim.fastq  --right V4_2_paired_trim.fastq --CPU 6 --output trinity_out_dir_OrgAnti

##Do kallisto expression quantification
kallisto index OrgAnti_transcriptome.fasta -i OrgAnti_transcriptome

kallisto quant -i OrgAnti_transcriptome -o OrgAnti_output --threads 20 V4_1_paired_trim.fastq V4_2_paired_trim.fastq
