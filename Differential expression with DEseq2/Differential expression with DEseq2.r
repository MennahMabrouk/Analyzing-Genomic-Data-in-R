# Load the libraries
library(DESeq2)
library(ggplot2)
library(dplyr)

# Read in the raw read counts
rawCounts <- read.delim("F:\\E-MTAB-9767-raw-counts.tsv")

# Read in the sample mappings
sampleData <- read.delim("F:\\E-MTAB-9767-experiment-design.tsv")

# Also save a copy for later
sampleData_v2 <- sampleData

# Plot the histogram of raw expression counts
ggplot(data = rawCounts) +
  geom_histogram(mapping = aes(x = ERR4843201), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Convert count data to a matrix of appropriate form that DESeq2 can read
geneID <- rawCounts$Gene.ID
sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[, sampleIndex])
rownames(rawCounts) <- geneID
head(rawCounts)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData) <- sampleData$Run
keep <- c("Sample.Characteristic.organism.part.", "Sample.Characteristic.individual.")
sampleData <- sampleData %>% 
  dplyr::select(!!keep[1], !!keep[2]) %>%
  dplyr::rename(tissueType = !!keep[1], individualID = !!keep[2]) %>%
  dplyr::mutate(individualID = factor(individualID))
head(sampleData)

# Put the columns of the count data in the same order as row names of the sample mapping
rawCounts <- rawCounts[, intersect(rownames(rawCounts), rownames(sampleData))]

# Match the row names in rawCounts with the column names in sampleData
colIndices <- match(rownames(sampleData), colnames(rawCounts))
rawCounts <- rawCounts[, colIndices]

# Check for NA values in the count matrix
has_na <- any(is.na(rawCounts))
if (has_na) {
  # Handle NA values by removing rows with NA
  na_rows <- rowSums(is.na(rawCounts)) > 0
  rawCounts <- rawCounts[!na_rows, ]
  sampleData <- sampleData[!na_rows, ]
}

# Make sure dimensions match after NA handling
if (ncol(rawCounts) != nrow(sampleData)) {
  stop("Dimensions of count data and sample information do not match.")
}

# Create the DESeq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData = rawCounts, colData = sampleData, design = ~ individualID + tissueType)
