# alpha-diversity
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


fin_table  = data.frame()
for (group in unique(sample_data(physeq)$new_group)){
	temp_physeq = prune_samples(sample_data(physeq)$new_group == group, physeq)

	new_table <- data.frame(taxa_sums(temp_physeq), tax_table(temp_physeq)[,2])
	colnames(new_table) = c("abundance", "phylum")
	sum_table = data.frame()
	for (x in unique(new_table$phylum) ){
		temp = subset(new_table, phylum==x)
		su = sum(temp$abundance)
		sum_table = rbind(sum_table, data.frame(su, temp[1,2]))
	}
	colnames(sum_table) = c("abundance", "phylum")
	perc_table = sum_table
	for (i in 1:nrow(sum_table)){
			perc_table[i,1] = sum_table[i,1] / sum(sum_table[,1])
	}
	
	other_table = data.frame()
	other = c()
	for (i in 1:nrow(perc_table)){
		if(perc_table[i,1] > 0.01) {
			other_table = rbind(other_table, perc_table[i,])
		}else{
			other = c(other, perc_table[i,1])
		}
		
	}
	
	sum(other)
	tep = data.frame(sum(other), "other")
	colnames(tep) = c("abundance", "phylum")
	tfin = rbind(other_table, tep)
	ttfin = cbind(tfin,group)
	fin_table = rbind(fin_table, ttfin)
	phy = unique(fin_table$phylum)
}

fin_table$group = factor(fin_table$group, levels = c("NlimHC", "NlimLC","PlimHC","PlimLC"))
pdf("alpha_diversity.pdf", width=6, height=6)
ggplot(fin_table, aes(x=group,y=abundance, fill=phylum))+geom_bar(stat="identity")+labs(y="relative abundance")
dev.off()

phy[1:7]


## memo below
#get_taxa_unique(temp_physeq, "Rank2")
#sample_sums(physeq)
#taxa_sums(temp_physeq)
#get_taxa(physeq, names)
#get_sample()



physeq <- readRDS("water-project-object_physeq_normalized.rds")
#relative abundance
rephyseq= transform_sample_counts(physeq, function(x) x / sum(x) )

# p__Proteobacteria  p__Actinobacteria  p__Cyanobacteria   p__Verrucomicrobia p__Chlorobi   p__Planctomycetes  p__Bacteroidetes 
fin_table  = data.frame()


phylo= "p__Bacteroidetes"
subphy <- subset_taxa(rephyseq, Rank2 == "p__Bacteroidetes" )


for (group in unique(sample_data(subphy)$new_group)){
	temp_physeq = prune_samples(sample_data(subphy)$new_group == group, subphy)

	new_table <- data.frame(taxa_sums(temp_physeq), tax_table(temp_physeq)[,6])
	colnames(new_table) = c("abundance", "genus")
	sum_table = data.frame()
	for (x in unique(new_table$genus) ){
		temp = subset(new_table, genus==x)
		su = sum(temp$abundance)
		sum_table = rbind(sum_table, data.frame(su, temp[1,2]))
	}
	colnames(sum_table) = c("abundance", "genus")
	perc_table = sum_table
	#for (i in 1:nrow(sum_table)){
	#		perc_table[i,1] = sum_table[i,1] / sum(sum_table[,1])
	#}
	
	other_table = data.frame()
	other = c()
	for (i in 1:nrow(perc_table)){
		if(perc_table[i,1] > 0.01) {
			other_table = rbind(other_table, perc_table[i,])
		}else{
			other = c(other, perc_table[i,1])
		}
		
	}
	
	sum(other)
	tep = data.frame(sum(other), "other")
	colnames(tep) = c("abundance", "genus")
	tfin = rbind(other_table, tep)
	ttfin = cbind(tfin,group,phylo)
	fin_table = rbind(fin_table, ttfin)
}

fin_table$group = factor(fin_table$group, levels = c("NlimHC", "NlimLC","PlimHC","PlimLC"))
filename = paste0(phylo,"_alpha_diversity.pdf" )
pdf(filename, width=6, height=6)
ggplot(fin_table, aes(x=group,y=abundance, fill=genus))+geom_bar(stat="identity")+labs(y="relative abundance")
dev.off()

saveRDS(fin_table, file="genus_table.rds")
fin_table = readRDS("genus_table.rds")
pdf("alpha_all_phyum.pdf", width=15, height=6)
ggplot(fin_table, aes(x=group,y=abundance, fill=genus))+geom_bar(stat="identity")+labs(y="relative abundance")+facet_grid(~phylo)
dev.off()

