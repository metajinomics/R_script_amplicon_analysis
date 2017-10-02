#Find core microbiome

setwd("~/Box Sync/2017/1January/water/water64/")


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

data <- as.data.frame(otu_table(physeq))
dataPA = (data>0)*1
ubiq = subset(dataPA, rowSums(dataPA) > 32) #64 x 50% make 426 taxa
names = rownames(ubiq)
physeq= prune_taxa(names, physeq)
#down to 423 taxa, 64 sample

#count(taxa_sums(physeq) < 286)
summary( sample_sums(physeq))
#summary( sample_sums(physeq))
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 7724   10600   15570   29340   22900  286300 
  
summary(taxa_sums(physeq))
# summary(taxa_sums(physeq))
#   Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
 #   48.22    378.10    985.00   4440.00   2714.00 227700.00 
 
#physeq = prune_taxa(taxa_sums(physeq) > 286, physeq)
#down to 347 taxa

#physeq = prune_taxa(taxa_sums(physeq) > 4000, physeq)

physeq = prune_taxa(taxa_sums(physeq) > 2714, physeq) #3rd qu
#106 taxa
saveRDS(physeq, "water-project-object_physeq_core.rds")