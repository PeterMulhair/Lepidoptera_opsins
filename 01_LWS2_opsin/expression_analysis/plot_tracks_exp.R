source("https://github.com/PoisonAlien/trackplot/blob/master/R/trackplot.R?raw=true")

bigWigs = c("OrgAnti_1instar_genome_cov_chrm_sorted.bw", "OrgAnti_male_genome_cov_chrm_sorted.bw", "OrgAnti_female_genome_cov_chrm_sorted.bw")

#Change chromosome number and location to what is required
track_data = track_extract(bigWigs = bigWigs, loci = "2:19363924-19381861")

pdf('OrgAnti_1instar_male_female_LWS2.pdf', height = 5, width = 15)
track_plot(summary_list = track_data, draw_gene_track = TRUE, gene_model = "Orgyia_antiqua-GCA_916999025.1-2022_02-genes.gtf.gz", isGTF = TRUE)
dev.off()

track_data = track_extract(bigWigs = bigWigs, loci = "1:12416716-12431743")

pdf('OrgAnti_1instar_male_female_LWS1.pdf', height = 5, width = 15)
track_plot(summary_list = track_data, draw_gene_track = TRUE, gene_model = "Orgyia_antiqua-GCA_916999025.1-2022_02-genes.gtf.gz", isGTF = TRUE)
dev.off()