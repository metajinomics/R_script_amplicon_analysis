# beta-diversity
setwd("~/")
##########   Load library   ############
library(ape)
library(ggplot2)
library(grid)
library(plyr)
library(reshape)
library(reshape2)

library(gridExtra)
library(corrplot)

library(vegan)


library(phyloseq)
library(DESeq2)

##########  USE saved object ############
physeq <- readRDS("water-project-object_physeq_normalized.rds")


##########  PCOA plot ############
qiime.ord <- ordinate(physeq, "PCoA","unifrac",weighted = TRUE)
p1 = plot_ordination(physeq, qiime.ord, color="new_group")
pdf("pcoa.pdf", width=6, height=6)
print(p1 + stat_ellipse())
dev.off()

cor = qiime.ord$vectors
met = sample_data(physeq)[,c(4,5,6,7,8,9,10,11,12,13,14)]
datEF = envfit(cor, met)
datEF.df = as.data.frame(datEF$vectors$arrows*sqrt(datEF$vectors$r))
datEF.df$species <- rownames(datEF.df)

pdf("pcoa_with_env.pdf", width=6, height=6)
p1+ geom_segment(data = datEF.df, aes(x=0, xend=Axis.1, y=0, yend=Axis.2), arrow = arrow(length=unit(0.2, "cm")), colour="grey")+geom_text(data= datEF.df, aes(x= Axis.1, y=Axis.2, label=species), size=5,colour="black") + stat_ellipse()
dev.off()


#NMDS
GP.ord <- ordinate(physeq, "NMDS", "bray")
p1 = plot_ordination(physeq, GP.ord, color="new_group")

pdf("NMDS.pdf", width=6, height=6)
print(p1+ stat_ellipse())
dev.off()

# some stastical numbers
d = phyloseq::distance(physeq, method="bray")
df <- data.frame(sample_data(physeq))
ado <- adonis(d ~ new_group, data=df)
capture.output(ado, file="adonis.txt")

b_result <- betadisper(d, df$new_group, "median")
ano_re <- anova(b_result)
capture.output(ano_re, file="anova.txt")

