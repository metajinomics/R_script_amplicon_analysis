library(phyloseq)
library(ggplot2)

otu_dada2 <- readRDS("seq_table.RDS")
tax_dada2 <- readRDS("tax_table.RDS")

ps <- phyloseq(otu_table(otu_dada2, taxa_are_rows=TRUE), tax_table(tax_dada2))


#relative abundance
rel <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))

#top 20 
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, ps.top20)
glom <- tax_glom(ps.top20, taxrank="phylum")
p <- plot_bar(glom, fill="phylum")
