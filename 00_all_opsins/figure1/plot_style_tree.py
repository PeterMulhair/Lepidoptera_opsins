import glob
import argparse
import toytree
import toyplot
import toyplot.svg
import toyplot.pdf
import numpy as np
from ete3 import Tree
from subprocess import call as unix

parse = argparse.ArgumentParser()

parse.add_argument("--genomes",type=str, help="path to genome fasta files",required=True)

args = parse.parse_args()


diurnal_sp = ['PapMach', 'PyrMalv', 'EryTage', 'ThySylv', 'HesComm', 'OchSylv', 'LepSina', 'ColCroc', 'AntCard', 'ApoCrat', 'PieRapa', 'PieBras', 'PieNapi', 'LycPhla', 'CelArgi', 'GlaAlex', 'PleArgu', 'CyaSemi', 'AriAges', 'LysBell', 'LysCori', 'ParAege', 'MelGala', 'EreLige', 'ManJurt', 'AphHype', 'LimCami', 'BolSele', 'FabAdip', 'MelCinx', 'MelAtha', 'VanCard', 'VanAtal', 'NymPoly', 'NymUrti', 'InaIo', 'ZygFili', 'NeoCorn', 'LimLuna', 'LimMarm', 'SesApif', 'SynVesp', 'AutGamm']

drop_sp = ['LimLuna', 'LimMarm']

change_sp_name = {'NemSwam': 'NemSwae', 'YpoSede': 'YpoSedl', 'HemFuci': 'HemFuc','EnnQuer': 'EnnQuei', 'PtiCapu': 'PtiCapc', 'EucMi': 'EucMixx', 'HerTars': 'HerTari', 'SpiLute': 'SpiLutu', 'EilDepr': 'EilDepe', 'AcrPsi': 'AcrPsix', 'XesC-Ni': 'XesCnig', 'AgrGeni': 'AgoAren'}


#Get list of species to plot ordered by the phylogeny
prune_list = []
with open('final_species_ordered.txt') as f:
    for line in f:
        lines = line.strip()
        if lines not in drop_sp:
            prune_list.append(lines)

#Parse rooted species tree
tree = Tree('rooted_sp_tree.nwk', format = 1)
tree.prune(prune_list)
tree.write(format=1, outfile="pruned_tree.nwk")

unix('cp pruned_tree.nwk pruned_tree_rename.nwk', shell=True)

sp_name_dict = {}
for genome in glob.glob(args.genomes + '*fasta'):
    if '_il' in genome:
        with open(genome) as f:
            first_line = f.readline()
            species = first_line.split(' ')[1] + '_' + first_line.split(' ')[2]
            sp = first_line.split(' ')[1][:3] + first_line.split(' ')[2][:4].title()
            if sp in change_sp_name:
                sp = change_sp_name[sp]
            sp_name_dict[sp] = species

#Convert short names to long names in new newick file
for species, sp in sp_name_dict.items():
    unix("sed -i 's/" + species + "/" + sp + "/g' pruned_tree_rename.nwk", shell=True)

#Get ultrametric tree using OrthoFinders make_ultrametric.py script (root age set to 300mya)
unix("python ~/bin/OrthoFinder/tools/make_ultrametric.py -r 300 pruned_tree_rename.nwk", shell=True)

#Import opsin gene counts per species (which has been ordered to match the phylogeny species list)
spdata = np.genfromtxt('opsin_gene_counts_reblast.csv', delimiter=',')
#Import newly made ultrametric tree
tree = toytree.tree("pruned_tree_rename.nwk.ultrametric.tre")

#Scale the plot to fit tree and data
ctree = tree.mod.node_scale_root_height(spdata.shape[1] * 2)

#Create colour variable for branches which are diurnal species
ecolors = ctree.get_edge_values_mapped({("Micropterix_aruncella", "Neomicropteryx_cornuta"): "#fec44f", ("Incurvaria_masculella", "Nematopogon_swammerdamellus"): "#fec44f", ("Zygaena_filipendulae"): "#fec44f", ("Sesia_apiformis", "Sesia_bembeciformis", "Bembecia_ichneumoniformis", "Synanthedon_andrenaeformis", "Synanthedon_myopaeformis", "Synanthedon_formicaeformis", "Synanthedon_vespiformis"): "#fec44f", ("Papilio_machaon", "Erynnis_tages", "Pyrgus_malvae", "Carterocephalus_palaemon", "Thymelicus_sylvestris", "Hesperia_comma", "Ochlodes_sylvanus", "Leptidea_sinapis", "Colias_croceus", "Anthocharis_cardamines", "Aporia_crataegi", "Pieris_rapae", "Pieris_brassicae", "Pieris_napi", "Lycaena_phlaeas", "Celastrina_argiolus", "Glaucopsyche_alexis", "Plebejus_argus", "Cyaniris_semiargus", "Aricia_agestis", "Aricia_artaxerxes", "Polyommatus_icarus", "Lysandra_bellargus", "Lysandra_coridon", "Bicyclus_anynana", "Lasiommata_megera", "Pararge_aegeria", "Hipparchia_semele", "Melanargia_galathea", "Aphantopus_hyperantus", "Maniola_jurtina", "Erebia_aethiops", "Erebia_ligea", "Limenitis_camilla", "Boloria_selene", "Fabriciana_adippe", "Mellicta_athalia", "Melitaea_cinxia", "Vanessa_atalanta", "Vanessa_cardui", "Nymphalis_polychloros", "Inachis_io", "Nymphalis_urticae"): "#fec44f", ("Euclidia_mi"): "#fec44f", ("Saturnia_pavonia"): "#fec44f", ("Apoda_limacodes"): "#7fc97f", ("Watsonalla_binaria"): "#7fc97f", ("Gymnoscelis_rufifasciata"): "#7fc97f", ("Orgyia_antiqua"): "#7fc97f", ("Catocala_fraxini"): "#7fc97f", ("Phragmatobia_fuliginosa"): "#7fc97f", ("Miltochrista_miniata"): "#7fc97f", ("Autographa_gamma"): "#7fc97f", ("Amphipoea_oculea"): "#7fc97f", ("Mesoligia_furuncula"): "#7fc97f"})

canvas = toyplot.Canvas(width=700, height=1000)
axes = canvas.cartesian()
axes.show = False

#Draw tree
ctree.draw(width=700, height=1000, edge_colors=ecolors, tip_labels_align=False, tip_labels=False, axes=axes);
colourlist = ['#bdbdbd', '#74a9cf', '#74c476', '#fc8d59', '#8c96c6'] 

#Draw columns of opsin copy number
ncols = 5
xoffset = 1.5
for col in range(5):
    data = spdata[:, col]
    axes.scatterplot(np.repeat(col, tree.ntips) + xoffset, np.arange(tree.ntips), marker='s', size=4, mstyle={'stroke': 'black, "stroke-width": 0.25'}, color=colourlist[col], opacity=0.1 + data[::-1] / data.max(), title=data,);

#Add rectangles to delimit Lep families
axes.rectangle(8, 8.2, 0, 51, color="#1c9099",)
axes.rectangle(8, 8.2, 52, 68, color="#a6bddb",)
axes.rectangle(8, 8.2, 69, 70, color="#1c9099",)
axes.rectangle(8, 8.2, 71, 78, color="#a6bddb",)
axes.rectangle(8, 8.2, 79, 112, color="#1c9099",)
axes.rectangle(8, 8.2, 113, 116, color="#a6bddb",)
axes.rectangle(8, 8.2, 117, 120, color="#1c9099",)
axes.rectangle(8, 8.2, 121, 130, color="#a6bddb",)
axes.rectangle(8, 8.2, 131, 132, color="#1c9099",)
axes.rectangle(8, 8.2, 137, 138, color="#a6bddb",)
axes.rectangle(8, 8.2, 139, 157, color="#1c9099",)
axes.rectangle(8, 8.2, 158, 167, color="#a6bddb",)
axes.rectangle(8, 8.2, 168, 174, color="#1c9099",)
axes.rectangle(8, 8.2, 175, 180, color="#a6bddb",)
axes.rectangle(8, 8.2, 182, 184, color="#1c9099",)
axes.rectangle(8, 8.2, 185, 191, color="#a6bddb",)
axes.rectangle(8, 8.2, 195, 209, color="#1c9099",)
axes.rectangle(8, 8.2, 214, 215, color="#a6bddb",)
axes.rectangle(8, 8.2, 216, 217, color="#1c9099",)
axes.rectangle(8, 8.2, 218, 219, color="#a6bddb",)

axes.x.domain.max = 20

#Plot tree
toyplot.svg.render(canvas, 'opsin_style_tree.svg')



####################################################################
#Code to plot Supplementary Figure 1
#################################################################### 

##Replot tree this time just with species names and diurnal/nocturnal branches
tree = toytree.tree("pruned_tree_rename.nwk.ultrametric.tre")
ctree = tree.mod.node_scale_root_height(1 * 1)

ecolors = ctree.get_edge_values_mapped({("Micropterix_aruncella", "Neomicropteryx_cornuta"): "#fec44f", ("Incurvaria_masculella", "Nematopogon_swammerdamellus"): "#fec44f", ("Zygaena_filipendulae"): "#fec44f", ("Sesia_apiformis", "Sesia_bembeciformis", "Bembecia_ichneumoniformis", "Synanthedon_andrenaeformis", "Synanthedon_myopaeformis", "Synanthedon_formicaeformis", "Synanthedon_vespiformis"): "#fec44f", ("Papilio_machaon", "Erynnis_tages", "Pyrgus_malvae", "Carterocephalus_palaemon", "Thymelicus_sylvestris", "Hesperia_comma", "Ochlodes_sylvanus", "Leptidea_sinapis", "Colias_croceus", "Anthocharis_cardamines", "Aporia_crataegi", "Pieris_rapae", "Pieris_brassicae", "Pieris_napi", "Lycaena_phlaeas", "Celastrina_argiolus", "Glaucopsyche_alexis", "Plebejus_argus", "Cyaniris_semiargus", "Aricia_agestis", "Aricia_artaxerxes", "Polyommatus_icarus", "Lysandra_bellargus", "Lysandra_coridon", "Bicyclus_anynana", "Lasiommata_megera", "Pararge_aegeria", "Hipparchia_semele", "Melanargia_galathea", "Aphantopus_hyperantus", "Maniola_jurtina", "Erebia_aethiops", "Erebia_ligea", "Limenitis_camilla", "Boloria_selene", "Fabriciana_adippe", "Mellicta_athalia", "Melitaea_cinxia", "Vanessa_atalanta", "Vanessa_cardui", "Nymphalis_polychloros", "Inachis_io", "Nymphalis_urticae"): "#fec44f", ("Euclidia_mi"): "#fec44f", ("Saturnia_pavonia"): "#fec44f", ("Apoda_limacodes"): "#7fc97f", ("Watsonalla_binaria"): "#7fc97f", ("Gymnoscelis_rufifasciata"): "#7fc97f", ("Orgyia_antiqua"): "#7fc97f", ("Catocala_fraxini"): "#7fc97f", ("Phragmatobia_fuliginosa"): "#7fc97f", ("Miltochrista_miniata"): "#7fc97f", ("Autographa_gamma"): "#7fc97f", ("Amphipoea_oculea"): "#7fc97f", ("Mesoligia_furuncula"): "#7fc97f"})

canvas = toyplot.Canvas(width=800, height=1000)
axes = canvas.cartesian()
axes.show = False

ctree.draw(width=800, height=1000, scalebar=True, edge_colors=ecolors, tip_labels_align=False, tip_labels=True, axes=axes, tip_labels_style={"font-size": "5px"});

toyplot.pdf.render(canvas, 'opsin_style_tree_names.pdf')
