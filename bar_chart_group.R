library(phyloseq)
library(ggplot2)

# You need to change "new_group" into your name of treatment
group_bar_chart <- function(physeq){
fin_table  = data.frame()
for (group in unique(sample_data(physeq)$new_group) ){
	temp_physeq = prune_samples( (sample_data(physeq)$new_group == group), physeq)

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
return(fin_table)
}

#read phyloseq object
physeq <- readRDS("water-project-object_physeq_normalized.rds")


fin_table <- group_bar_chart(physeq)
fin_table$group = factor(fin_table$group, levels = c("NlimHC", "NlimLC","PlimHC","PlimLC"))
ggplot(fin_table, aes(x=group,y=abundance, fill=phylum))+geom_bar(stat="identity")+labs(y="relative abundance")
