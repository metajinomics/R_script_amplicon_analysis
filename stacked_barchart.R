#install phyloseq
#source("http://bioconductor.org/biocLite.R") 
#biocLite('phyloseq') 

library(phyloseq)
library(ggplot2)
#read data from mothur output
mothur_data <- import_mothur(mothur_shared_file = "stability.opti_mcc.shared",mothur_constaxonomy_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy")
meta <- read.table("meta.txt", header=T)
sample_data(mothur_data) <- meta


##rarefy the samples to get rids of thoese with low reads
set.seed(1)
rarefied <- rarefy_even_depth(mothur_data, sample.size = 9000,rngseed=TRUE)


glom <- tax_glom(rarefied , taxrank="Rank2")
psmelt <- psmelt(glom)
ggplot(psmelt, aes(x=Sample, y=Abundance, fill=Rank2))+geom_bar(stat="identity", position = "fill")+facet_grid(.~new_group, scales= "free")+theme(axis.text.x = element_text(angle=90))
