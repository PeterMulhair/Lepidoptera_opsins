library(ggplot2)
library(patchwork)

lep_UV_FitMG94_output <- read.delim("lep_UV_FitMG94_output.tsv", header=FALSE)

p<-ggplot(lep_UV_FitMG94_output, aes(y=V4, x=V1, fill=V1)) + 
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.1, alpha=0.5) +
  scale_fill_manual(values = c('Diurnal' = '#fec44f', 'Nocturnal' = '#4d004b', 'Both' = '#fcc5c0')) +
  theme_classic() + 
  ylab('dN') +
  scale_x_discrete() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 20), legend.position = "none")

p1<-ggplot(lep_UV_FitMG94_output, aes(y=V3, x=V1, fill = V1)) + 
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.1, alpha=0.5) +
  scale_fill_manual(values = c('Diurnal' = '#fec44f', 'Nocturnal' = '#4d004b', 'Both' = '#fcc5c0')) +
  theme_classic() + 
  ylab('dN/dS') +
  scale_x_discrete() +
  theme(axis.text.x = element_text(size = 15), axis.title.x = element_blank(), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 20))


plot<-p/p1

pdf('DnDs_photic_niche.pdf', height = 6, width = 7)
print(plot)
dev.off()


#Plot rainplots
p<-ggplot(lep_UV_FitMG94_output,aes(x=V1,y=V3, fill = V1, colour = V1))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE, alpha = 0.7)+
  geom_point(position = position_jitter(width = .1), size = .25)+
  geom_boxplot(aes(x = V1, y = V3),position = position_nudge(x = .25, y = 0),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  ylab('dN/dS')+xlab(NULL)+coord_flip()+theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = c('Diurnal' = '#fec44f', 'Nocturnal' = '#4d004b', 'Both' = '#7fc97f')) +
  scale_fill_manual(values = c('Diurnal' = '#fec44f', 'Nocturnal' = '#4d004b', 'Both' = '#7fc97f')) 


p1<-ggplot(lep_UV_FitMG94_output,aes(x=V1,y=V4, fill = V1, colour = V1))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE, alpha = 0.7)+
  geom_point(position = position_jitter(width = .1), size = .25)+
  geom_boxplot(aes(x = V1, y = V4),position = position_nudge(x = .25, y = 0),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  ylab('dN')+xlab(NULL)+coord_flip()+theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = c('Diurnal' = '#fec44f', 'Nocturnal' = '#4d004b', 'Both' = '#7fc97f')) +
  scale_fill_manual(values = c('Diurnal' = '#fec44f', 'Nocturnal' = '#4d004b', 'Both' = '#7fc97f')) 


plot<-(p1 + ggtitle("UV"))/p

pdf('DnDs_photic_niche_rainplot.pdf', height = 6, width = 5)
print(plot)
dev.off()




#Stats
X<-split(lep_UV_FitMG94_output, lep_UV_FitMG94_output$V1)
Both<-X$Both
Nocturnal<-X$Nocturnal
Diurnal<-X$Diurnal

Both_Nocturnal <- wilcox.test(Both$V3, Nocturnal$V3) 
Both_Diurnal <- wilcox.test(Both$V3, Diurnal$V3) 
Nocturnal_Diurnal <- wilcox.test(Nocturnal$V3, Diurnal$V3) 
