from ete3 import Tree
import toytree
import toyplot
import toyplot.pdf
import toyplot.svg

#List of species in Noctuoidea superfamily
sp_list = ['MegAlbu', 'NycReva', 'LymMona', 'EupSimi', 'OrgAnti', 'CatFrax', 'EucMixx', 'SchCost', 'LasFlex', 'TriEmor', 'HypProb', 'HerTari', 'PhrFuli', 'SpiLubr', 'SpiLutu', 'MilMini', 'CybMeso', 'EilDepe', 'EilSoro', 'AbrTril', 'AbrTrip', 'DiaChry', 'AutGamm', 'AutPulc', 'ProPyga', 'AllOxya', 'CraLigu', 'AcrPsix', 'AcrAcer', 'AcrLepo', 'XylAreo', 'AmpBerb', 'AmpTrag', 'SpoExig', 'SpoFrug', 'CarClav', 'AnoMund', 'ThoDeci', 'HecDyso', 'MamBras', 'MytImpu', 'MytAlbi', 'MytFerr', 'AgrPuta', 'OchPlec', 'DiaRubi', 'NocPron', 'NocFimb', 'NocJant', 'XesCnig', 'XesSexs', 'XesXant', 'EupLuci', 'PhlMeti', 'AteCent', 'CosPyra', 'CosTrap', 'EupTran', 'OmpLuno', 'AgrCirc', 'AgrMaci', 'BraVimi', 'DryErem', 'GriApri', 'ApoLuen', 'ApaMono', 'ApaSord', 'AmpOcul', 'HydMica', 'LupTest', 'MesFuru']

#List of species in Erebidae and Noctuidae families
erebid_list = ['LymMona', 'EupSimi', 'OrgAnti', 'CatFrax', 'EucMixx', 'SchCost', 'LasFlex', 'TriEmor', 'HypProb', 'HerTari', 'PhrFuli', 'SpiLubr', 'SpiLutu', 'MilMini', 'CybMeso', 'EilDepe', 'EilSoro']
noctuid_list = ['AbrTril', 'AbrTrip', 'DiaChry', 'AutGamm', 'AutPulc', 'ProPyga', 'AllOxya', 'CraLigu', 'AcrPsix', 'AcrAcer', 'AcrLepo', 'XylAreo', 'AmpBerb', 'AmpTrag', 'SpoExig', 'SpoFrug', 'CarClav', 'AnoMund', 'ThoDeci', 'HecDyso', 'MamBras', 'MytImpu', 'MytAlbi', 'MytFerr', 'AgrPuta', 'OchPlec', 'DiaRubi', 'NocPron', 'NocFimb', 'NocJant', 'XesCnig', 'XesSexs', 'XesXant', 'EupLuci', 'PhlMeti', 'AteCent', 'CosPyra', 'CosTrap', 'EupTran', 'OmpLuno', 'AgrCirc', 'AgrMaci', 'BraVimi', 'DryErem', 'GriApri', 'ApoLuen', 'ApaMono', 'ApaSord', 'AmpOcul', 'HydMica', 'LupTest', 'MesFuru']


gene_prune = ['TinTrin_LWS1']
tree = Tree('lep_opsins.mft.contree.rooted')
for node in tree.traverse("postorder"):
    if node.is_leaf():
        sp = node.name.split('_')[0]
        if (sp in sp_list) and ('LWS' in node.name):
            gene_prune.append(node.name)

#Prune full opsin gene tree to just include Noctuoidea LW opsin genes
tree.prune(gene_prune)
tree.set_outgroup('TinTrin_LWS1')
tree.write(format=1, outfile="Noctuoidea_LWS_tree.nwk")

#Draw Noctuoidea LW tree
rtre = toytree.tree('Noctuoidea_LWS_tree.nwk')

#Colour branches by family
colorlist = ["#02818a" if tip.split('_')[0] in erebid_list else ("#9dc6c5" if tip.split('_')[0] in noctuid_list else ("black" if tip.split('_')[0] == "TinTrin" else '#d4e8e7') for tip in rtre.get_tip_labels()]

canvas, axes, mark  = rtre.draw(tip_labels_colors = colorlist, tip_labels_align=True, width=500, height=800, tip_labels_style={"font-size": "8px"});

#Plot gene tree
toyplot.pdf.render(canvas, 'Noctuoidea_LWS_tree.pdf')
