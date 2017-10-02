#save data table for co-occurence

setwd("~/Box Sync/2017/1January/water/water64/")
##########   Load library   ############
library(phyloseq)

##########  USE saved object ############
physeq <- readRDS("water-project-object_physeq_normalized.rds")
#physeq <- readRDS("water-project-object_physeq_core.rds")
meta <- as.data.frame(sample_data(physeq))
count <- as.data.frame(otu_table(physeq))
ftable <- cbind(meta, t(count))

#write.csv(ftable, "table_for_co_core64.csv")
write.csv(ftable, "table_for_co_water_full.csv")