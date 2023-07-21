library(DESeq2)
library(pheatmap)

# Read in the raw read counts
rawCounts <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-raw-counts.tsv")
head(rawCounts)

# Read in the sample mappings
sampleData <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-experiment-design.tsv")
head(sampleData)

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Gene.ID
sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[, sampleIndex])
rownames(rawCounts) <- geneID
head(rawCounts)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
head(sampleData)
rownames(sampleData) <- sampleData$Run
keep <- c("Sample.Characteristic.biopsy.site.", "Sample.Characteristic.individual.")
sampleData <- sampleData[, keep]
colnames(sampleData) <- c("tissueType", "individualID")
sampleData$individualID <- factor(sampleData$individualID)
head(sampleData)

# Put the columns of the count data in the same order as row names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[, unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

# Rename the tissue types
rename_tissues <- function(x) {
  x <- switch(as.character(x), "normal" = "normal-looking surrounding colonic epithelium", "primary tumor" = "primary colorectal cancer", "colorectal cancer metastatic in the liver" = "metastatic colorectal cancer to the liver")
  return(x)
}
sampleData$tissueType <- unlist(lapply(sampleData$tissueType, rename_tissues))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
sampleData$tissueType <- factor(sampleData$tissueType, levels = c("normal-looking surrounding colonic epithelium", "primary colorectal cancer", "metastatic colorectal cancer to the liver"))

# Modify factor levels to comply with safe naming conventions
levels(sampleData$individualID) <- gsub("[^A-Za-z0-9_.]", "_", levels(sampleData$individualID))
levels(sampleData$tissueType) <- gsub("[^A-Za-z0-9_.]", "_", levels(sampleData$tissueType))

# Create the DESeq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData = rawCounts, colData = sampleData, design = ~ individualID + tissueType)

# Estimate size factors
dds_wt <- estimateSizeFactors(deseq2Data)

# Unsupervised clustering analysis: log transformation using vst
vsd_wt <- vst(dds_wt, blind = TRUE)

# Hierarchical clustering with correlation heatmaps
vsd_mat_wt <- assay(vsd_wt)
vsd_cor_wt <- cor(vsd_mat_wt)

# Add the ggplot code snippet with modified x-axis formatting
ggplot(data.frame(wt_normal1 = rawCounts[, 1])) +
  geom_histogram(aes(x = wt_normal1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE))

# Prepare data for pheatmap
data_for_heatmap <- as.matrix(vsd_cor_wt)

# Convert tissueType to a character vector
annotation_row <- as.character(sampleData$tissueType)

# Add spaces between words in the x-axis labels
annotation_row_with_spaces <- paste(" ", annotation_row, " ")

# Plot the heatmap using pheatmap with manual row annotations
pheatmap(data_for_heatmap,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_rownames = FALSE,
         show_colnames = TRUE,
         row_names_side = "left",
         annotation_colors = "black",
         annotation_names_row = FALSE,
         labels_row = annotation_row_with_spaces,
         fontsize_row = 8,     # Adjust the font size of row labels
         fontsize_col = 12,    # Adjust the font size of column labels
         angle_col = 45)       # Set the angle of column labels to 45 degrees

# Calculate PCA scores
sample_scores <- as.data.frame(assay(vsd_wt))
sample_scores$Sample <- rownames(sample_scores)

column_names <- colnames(vsd_wt)
colnames(sample_scores)[2:5] <- c("normal", "fibrosis", "tumor", "metastasis")

sample_scores$PC1 <- sample_scores$normal * -2 + sample_scores$fibrosis * -10 + sample_scores$tumor * 8 + sample_scores$metastasis * 1
sample_scores$PC2 <- sample_scores$normal * 0.5 + sample_scores$fibrosis * 1 + sample_scores$tumor * -5 + sample_scores$metastasis * 6

# Print the PCA scores
print(sample_scores)


