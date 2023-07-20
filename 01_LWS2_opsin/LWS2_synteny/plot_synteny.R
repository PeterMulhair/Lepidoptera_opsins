library(ggplot2)
library(gggenes)

#Code to plot the LWS2 synteny cluster in Nolidae, Noctuidae, and Erebidae (Supplementary Figure S6).

LWS2_synteny_ereb <- read.delim("LWS2_synteny_ereb.tsv")

dummies <- make_alignment_dummies(LWS2_synteny_ereb,aes(xmin = start, xmax = end, y = Molecule, id = gene),on = "RabGAP-TBC_gene")

p<-ggplot(LWS2_synteny_ereb, aes(xmin = start, xmax = end, y = Molecule, fill = gene)) +
  facet_wrap(~ Molecule, scales = "free", ncol = 1) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  geom_feature(data = LWS2_synteny_ereb, aes(x = from, y = Molecule, forward = direction)) +
  geom_feature_label(data = LWS2_synteny_ereb, aes(x = from, y = Molecule, label = "LWS2", forward = direction)) +
  geom_subgene_arrow(data = LWS2_synteny_ereb,
    aes(xmin = start, xmax = end, y = Molecule, xsubmin = from, xsubmax = to), fill = "white", color="black", alpha=.9) +
  theme_genes() +
  scale_fill_manual(values = c("Fis1_TPR" = "#a6cee3", "Kiaa" = "#1f78b4", "hypothetical" = "#b2df8a", "RabGAP-TBC_gene" = "#33a02c", "Metallophos" = "#fb9a99", "Cyclin" = "#e31a1c", "HMG_box" = "#fdbf6f")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size = 20))


pdf("Erebidae_LWS2_synteny.pdf", height = 5, width = 7)
print(p)
dev.off()


LWS2_synteny <- read.delim("LWS2_synteny_noct.tsv")

dummies <- make_alignment_dummies(LWS2_synteny,aes(xmin = start, xmax = end, y = Molecule, id = gene),on = "LWS2")

p<-ggplot(LWS2_synteny, aes(xmin = start, xmax = end, y = Molecule, fill = gene)) +
  facet_wrap(~ Molecule, scales = "free", ncol = 1) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  geom_subgene_arrow(data = LWS2_synteny,
    aes(xmin = start, xmax = end, y = Molecule, xsubmin = from, xsubmax = to), fill = "white", color="black", alpha=.9) +
  theme_genes() +
  scale_fill_manual(values = c("hypothetical" = "#e0e0e0", "zf_Hakai_gene" = "#8dd3c7", "YjeF_PCNA_gene" = "#ffffb3", "PCNA_gene" = "#bebada", "INTS5" = "#fb8072", "Septin" = "#80b1d3", "ATP_bind_1_gene" = "#fdb462", "LWS2" = "#b3de69", "MFS_1" = "#fccde5")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size = 20))



pdf("Noctuidae_LWS2_synteny.pdf", height = 5, width = 7)
print(p)
dev.off()




LWS2_synteny_nolid <- read.delim("LWS2_synteny_nolid.tsv")

dummies <- make_alignment_dummies(LWS2_synteny_nolid,aes(xmin = start, xmax = end, y = Molecule, id = gene),on = "LWS2")

LWS2_synteny_nolid$Molecule <- factor(LWS2_synteny_nolid$Molecule, levels = unique(LWS2_synteny_nolid$Molecule))
LWS2_synteny_nolid$gene <- factor(LWS2_synteny_nolid$gene, levels = unique(LWS2_synteny_nolid$gene))

p<-ggplot(LWS2_synteny_nolid, aes(xmin = start, xmax = end, y = factor(Molecule, levels = unique(Molecule)), fill = gene)) +
  facet_wrap(~ Molecule, scales = "free", ncol = 1) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  geom_subgene_arrow(data = LWS2_synteny_nolid,
    aes(xmin = start, xmax = end, y = Molecule, xsubmin = from, xsubmax = to), fill = "white", color="black", alpha=.9) +
  theme_genes() +
  scale_fill_manual(values = c("hypothetical" = "#e0e0e0", "LRR_8" = "#cb4f42", "WW_domain" = "#c98443", "LWS2" = "#9a963f", "AMP_synthetase" = "#64ac48", "C-type_lectin" = "#4bad90", "Willebrand_factor" = "#6e7ecb", "Fork_head_domain" = "#b45ac2")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size = 20))



pdf("Nolidae_LWS2_synteny.pdf", height = 5, width = 7)
print(p)
dev.off()

