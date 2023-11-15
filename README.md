# Lepidoptera opsin evolution
This repository contains the datasets and scripts needed to reproduce the results in [Opsin gene duplication in Lepidoptera: retrotransposition, sex linkage, and gene expression](https://academic.oup.com/mbe/article/40/11/msad241/7341929)

Each folder contains data and code required to recreate results and figures in each section of the manuscript.

If you use any code or data from this repository please cite: [Mulhair et al. 2023](https://academic.oup.com/mbe/article/40/11/msad241/7341929)

If you have any questions about any of the code or data please get in touch at peter.mulhair[at]biology.ox.ac.uk

## Instructions

* `00_all_opsins/` contains files on the opsin genes present in each of the species. This include fasta files of the opsin protein sequences per species, as well as gff files per species with the location of each opsin in the genome annotated. It also contains code and data required to reproduce plots from Figure 1, including species tree used and opsin copy number per species.

* `01_LWS2_opsin/` contains code and data to run selective pressure analyses and expression analyses on the LWS2 opsin gene in Noctuoidea. This includes protein and nucleotide fasta files for opsin genes and all code required to reproduce plots from Figure 2. It also contains files on LWS2 synteny in each family and code used to plot synteny blocks (using [gggenes](https://github.com/wilkox/gggenes/tree/master)).

* `02_selective_pressure_analyses/` contains code and data required to run selective pressure analyses on all opsin genes and compare between day-flying and night-flying lineages. This includes all alignemnts and phylogenies required to run the analysis, as well as the code used. All analyses were carried out with [hyphy](https://github.com/veg/hyphy) or [vespasian](https://github.com/bede/vespasian).

---

Genomes used in this analysis from the [Darwin Tree of Life project](https://www.darwintreeoflife.org/) can be downloaded by using code from [here](https://github.com/PeterMulhair/DToL_insects)

---

<div align="center">
<p align="center">
<img src="https://github.com/PeterMulhair/Lepidoptera_opsins/blob/main/opsin_logo.png" width="450" height="540">
</p>
</div>