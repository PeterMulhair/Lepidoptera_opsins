library(ape)
library(ggtree)

tree<- read.tree('rooted_sp_tree.nwk')

tip <- c("LimLuna", "LimMarm")
tree_trim<- drop.tip(tree, tip, trim.internal = TRUE)

p<-ggtree(tree_trim, layout = "circular") + geom_tiplab(size = 0, align=TRUE, linetype='dashed', linesize=.3) + geom_strip('LymMona', 'EilSoro', barsize=2, color='black', label="Erebidae", offset = 0.02, offset.text=.2) +
geom_strip('MicArun', 'NeoCorn', barsize=2, color='black', label="Micropterigidae", offset = 0.02, offset.text=.2) +
geom_strip('MegAlbu', 'NycReva', barsize=2, color='black', label="Nolidae", offset = 0.02, offset.text=.2) +
geom_strip('LepSina', 'PieNapi', barsize=2, color='black', label="Pieridae", offset = 0.02, offset.text=.2) +
geom_strip('AbrTril', 'MesFuru', barsize=2, color='black', label="Noctuidae", offset = 0.02, offset.text=.2) +
geom_strip('CloCurt', 'NotZicz', barsize=2, color='black', label="Notodontidae", offset = 0.02, offset.text=.2) +
geom_strip('LigAdus', 'TheBrit', barsize=2, color='black', label="Geometridae", offset = 0.02, offset.text=.2) +
geom_strip('DeiPorc', 'MimTili', barsize=2, color='black', label="Sphingidae", offset = 0.02, offset.text=.2) +
geom_strip('AcrSuav', 'AgrTris', barsize=2, color='black', label="Pyralidae", offset = 0.02, offset.text=.2) +
geom_strip('BicAnyn', 'NymUrti', barsize=2, color='black', label="Nymphalidae", offset = 0.02, offset.text=.2) +
geom_strip('LysCori', 'LycPhla', barsize=2, color='black', label="Lycaenidae", offset = 0.02, offset.text=.2) +
geom_strip('OchSylv', 'EryTage', barsize=2, color='black', label="Hesperiidae", offset = 0.02, offset.text=.2) +
geom_strip('EpiDema', 'PanCinn', barsize=2, color='black', label="Tortricidae", offset = 0.02, offset.text=.2) +
geom_strip('TinTrin', 'TinSemi', barsize=2, color='black', label="Tineidae", offset = 0.02, offset.text=.2) +
geom_strip('SesApif', 'SynVesp', barsize=2, color='black', label="Sesiidae", offset = 0.02, offset.text=.2) +
geom_strip('SteBipu', 'EmmMono', barsize=2, color='black', label="Pterophoridae", offset = 0.02, offset.text=.2) +
geom_strip('DreFalc', 'ThyBati', barsize=2, color='black', label="Drepanidae", offset = 0.02, offset.text=.2) +
geom_strip('YpsSequ', 'YpoSedl', barsize=2, color='black', offset = 0.02) +
geom_strip('NemSwae', 'IncMasc', barsize=2, color='black', offset = 0.02) +
geom_strip('ApoLima', 'ZygFili', barsize=2, color='black', offset = 0.02) +
geom_strip('AgoAren', 'BlaLact', barsize=2, color='black', offset = 0.02)


pdf('lep_tree.pdf', height = 10, width = 10)
print(p)
dev.off()

