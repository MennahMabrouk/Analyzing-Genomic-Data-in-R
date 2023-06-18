"""
IRanges
As you've learned in the previous chapters, 
you can store sequences with their own alphabets, order, 
and focus on certain intervals of these sequences. 
To extract sequence intervals, you use ranges. 
The Bioconductor package IRanges comes in handy with its function IRanges(), 
which creates a vector representation of a sequence, 
used to facilitate subsetting and annotation.

The IRanges package has been already loaded. 
To help you with this exercise you can check the documentation of IRanges().

Fill in the blank:

IRanges objects can be defined using two data types: ___ or ___ vectors.

"""
#logical , numeric

###########################################################

#Interacting with IRanges
#One
# Create the first sequence seq_1
seq_1 <- IRanges(start = 10, end = 37)

# Create the second sequence seq_2
seq_2 <- IRanges(start = c(5, 35, 50),
                 end = c(12, 39, 61),
                 names = LETTERS[1:3])

# Check the width of seq_1 and seq_2
width(seq_1)
width(seq_2)

#Two
# Create the first sequence seq_1
seq_1 <- IRanges(start = 10, end = 37)

# Create the second sequence seq_2
seq_2 <- IRanges(start = c(5, 35, 50),
                 end = c(12, 39, 61),
                 names = LETTERS[1:3])

# Check the width of seq_1 and seq_2
lengths(seq_1)
lengths(seq_2)

########################################################

#From tabular data to Genomic Ranges
# Load GenomicRanges package
library(GenomicRanges)

# Create seq_intervals
print(seq_intervals)


# Convert seq_intervals to GRanges using as()
myGR <- as(seq_intervals, "GRanges")

# Print myGR
print(myGR)

######################################################

#GenomicRanges accessors
# Load GenomicRanges
library(GenomicRanges)

# Print the seqinfo of myGR
seqinfo(myGR)

# Check the metadata
mcols(myGR)

####################################################

"""
ABCD1 mutation
You have just learned about the gene ABCD1. 
It encodes the protein in charge of the normal transport of 
fats that keep brain and lung cells functioning normally. 
When these groups of fats are not broken down, they build up in the body 
and become toxic. This affects the adrenal glands 
(small glands on top of each kidney) and the insulation 
(myelin) that surrounds neurons, causing hormonal problems 
and deteriorating movement, vision, and hearing. 
More than 650 mutations in the ABCD1 gene have been found to 
cause X-linked adrenoleukodystrophy, a rare genetic disease.

Since you are going to be studying this gene in the coming exercises, 
it is important to remember where it is located. If you are unsure, 
check the gene ABCD1 and its location using the Ensembl genome browser.

Where is ABCD1 located?

"""
# Chromosome X

##################################################################
#Human genome chromosome X
# Load human reference genome hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Assign hg38 to hg, then print it
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg

#Two
# Extract all the genes in chromosome X as hg_chrXg, then print it
hg_chrXg <- genes(hg, filter = list(tx_chrom = c("chrX")))
hg_chrXg

#Three
# Extract all positive stranded genes in chromosome X, assign to #hg_chrXgp, then sort it

hg_chrXgp <- genes(hg, filter = list(tx_chrom = c("chrX"), tx_strand = "+"))
sort(hg_chrXgp)

##########################################################

"""
A sequence window
To temporarily partition sections of your sequences, 
you will use the concept of windows of given widths that can move in steps. 
As you have seen in the video, GRanges provides the slidingWindows() function, 
with arguments such as width and step.

slidingWindows(x, width, step = 1L)
The GRanges object called ABCD1, with gene id 215 and length 19895, 
has been pre-loaded for this exercise. Use the ranges() 
function to see its structure.

ranges(ABCD1)
If you needed exactly 2 windows using step = 1L, 
what will be the maximum number allowed for the width of each window?
"""
#GRanges Objects -1 

#######################################
#More about ABCD1
#One
# Store the overlapping range in rangefound
rangefound <- subsetByOverlaps(hg_chrX, ABCD1)

# Print names of rangefound
names(rangefound)

#Two
# Print the gene of interest 
print(ABCD1)

# Print rangefound
print(rangefound)

########################################

#How many transcripts?
# Load the human transcripts DB to hg
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Prefilter chromosome X "chrX" using seqlevels()
seqlevels(hg) <- c("chrX")

# Get all transcripts by gene and print it
hg_chrXt <- transcriptsBy(hg, by = "gene")
hg_chrXt

# Select gene `215` from the hg_chrXt
hg_chrXt$`215`

#############################################
#From GRangesList object into a GRanges object
#One
# Unlist hg_ChrX and save result as myGR
myGR <- unlist(hg_ChrX)

# Compare classes of hg_ChrX and myGR
lengths(hg_ChrX)
lengths(myGR)

#Two
# Unlist hg_ChrX and save result as myGR
myGR <- unlist(hg_ChrX)

# Compare classes of hg_ChrX and myGR
class(hg_ChrX)
class(myGR)

# Compare length of hg_ChrX and myGR
lengths(hg_ChrX)
lengths(myGR)

#Three
"""
Question
Fill in the blank:

You would expect the GRanges object to have a length ___ the GRangesList.
"""
# Longer