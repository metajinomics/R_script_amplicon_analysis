# This script read result of QIIME to phylosq object


##########   Load library   ############
library(ape)
library(phyloseq)

##########  Read data  ############
biom_file = "otu_table_mc2_w_tax.json.biom"
biom_otu_tax <- import_biom(biom_file)
tree = read.tree("rep_set.tre")
meta  <- read.csv("meta_w_new_group.new64.csv",header=T, row.names=1)
sampledata = sample_data(meta)
physeq = merge_phyloseq(sampledata, biom_otu_tax, tree)

########## Save physeq object ##########
saveRDS(physeq, "water-project-object_physeq_before_normalization.rds")


