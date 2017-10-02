

##########   Load library   ############
library(phyloseq)
library(DESeq2)

##########  Read data  ############
physeq <- readRDS("water-project-object_physeq_before_normalization.rds")

##########  Normalizatoin   ############
countData <- as.data.frame(as.matrix(otu_table(physeq)))
colData <- as.data.frame(sample_data(physeq))

## Create a DESeqDataSet object
deseq2_deseqds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ new_group)

deseq2_de <- DESeq(deseq2_deseqds, fitType = "local")
pdf("dirpersion_plot.pdf", width=6, height=6)
plotDispEsts(deseq2_de, main="Dispersion plot")
dev.off()

########## Save physeq object with normalized value ##########
count <- counts(deseq2_de, normalized=TRUE)
otu_table(physeq) <- otu_table(count, taxa_are_rows=TRUE)
saveRDS(physeq, "water-project-object_physeq_normalized.rds")

