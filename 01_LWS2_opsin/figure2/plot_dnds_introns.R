library(ggplot2)

hyphy_results_LWS1 <- read.delim("hyphy_results_LWS_Noct.tsv", header=FALSE)

plot<-ggplot(hyphy_results_LWS1, aes(V6, V3))+
  geom_point(aes(colour = V6)) +
  geom_path(aes(group=V1), size=0.05) + 
  scale_colour_gradient(low = "#FCCF31", high = "#b30000", name = "Intron count") +
  ylim(0,0.2) +
  ggtitle("LWS multicopy genes") +
  scale_x_continuous(breaks=seq(0, 7, 7)) +
  xlab("Intron count") +
  ylab("dN/dS") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 20), plot.title = element_text(size = 20, face = "bold"))


pdf("LWS_dnds_Noct.pdf", height = 5, width = 6)
print(plot)
dev.off()

