"""
Why fastq?
So far, you have been introduced to the ShortRead package with useful tools
for the initial stages of short-read sequence analysis. 
As a summary, the main functionalities of this package are:

Data input
Quality assessment
Data transformation
Access to downstream analysis
In this context, you have been exposed to two different file types: 
fasta and fastq.

What do you think is the reason to use fastq formatted files?
"""
# ALL

############################################

"""
Reading in files
From the video, you've learned the difference between fasta and fastq files,
and what information can be stored in those files.

You have also seen examples of reading fasta and fastq files with their 
respective functions. You have learned that readFasta() and readFastq() 
need both a location and a file pattern to read one or many files.

Which two arguments are particular to the function readFastq() 
from the package ShortRead?
"""
# Qualitytype and Filter 

#############################################

#Exploring a fastq file
#One
# Load ShortRead
library(ShortRead)

# Print fqsample
print(fqsample)

#Two
# Load ShortRead
library(ShortRead)

# Print fqsample
fqsample

# Check class of fqsample
class(fqsample)

# Check class sread fqsample
class(sread(fqsample))

#Three
# Check ids of fqsample
id(fqsample)

#######################################
#Extract a sample from a fastq file
# Load ShortRead
library(ShortRead)

# Set a seed for sampling
set.seed(1234)

# Use FastqSampler with f and select 100 reads
fs <- FastqSampler(con = f, n = 100)

# Generate new sample yield
my_sample <- yield(fs)

# Print my_sample
print(my_sample)

##########################################

#Exploring sequence quality
#One
# load ShortRead
library(ShortRead)

# Check quality
quality(fqsample)

#Two
# Check encoding of quality
encoding(quality(fqsample))

#Three
# load ShortRead
library(ShortRead)

# Check quality
quality(fqsample)

# Check encoding of quality
encoding(quality(fqsample))

# Check baseQuality
qaSummary[["baseQuality"]]

#####################################

"""
Base quality plot
This figure has been created for you from the base 
quality encodings of a complete fastq file. 
This file uses Illumina encoding. 
The maximum encoding is I, or a score of 40. 
Good quality encodings have scores of 33 or above, or in other words, 
are B or above. baseQuality.png

You can see encoding details for Illumina here.

Which scores are represented most in the figure?

"""
# 37, 38 and 39

#####################################

# Glimpse nucByCycle
glimpse(nucByCycle)

# Create a line plot of cycle vs. count
nucByCycle %>% 
  # Gather the nucleotide letters in alphabet and get a new count column
  pivot_longer(-cycle, names_to = "alphabet", values_to = "count") %>% 
  ggplot(aes(x = cycle, y =  count, color = alphabet)) +
  geom_line(size = 0.5 ) +
  labs(y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())

#######################################

#Filtering reads on the go!
#One
# Load package ShortRead
library(ShortRead)

# Check class of fqsample
class(fqsample)

# Filter reads into selectedReads using myStartFilter
selectedReads <- fqsample[myStartFilter(fqsample)]

# Check class of selectedReads
class(selectedReads)

#Two
# Check detail of selectedReads
detail(selectedReads)

#######################################

"""
It is always a good practice to check that your sequence 
reads don't contain too many duplicates.

# Sample with duplicates of class: ShortReadQ
dfqsample

# Get the reads from dfqsample
mydReads <- sread(dfqsample)

# Counting duplicates
table(srduplicated(mydReads))
How would you go about removing duplicated reads in a file? 
Pay attention to what the condition should be in this filter.

"""
# == False

#################################################

#More filtering!
#One
# Check reads of fqsample
sread(fqsample)

# Create myFil using polynFilter
myFil  <- polynFilter(threshold= 3 , nuc = c("A"))

# Check myFil
myFil

#Two
# Check reads of fqsample
sread(fqsample)

# Create myFil using polynFilter
myFil <- polynFilter(threshold = 3, nuc = c("A"))

# Apply your filter to fqsample
filterCondition <- myFil(fqsample)

# Use myFil with fqsample
filteredSequences <- fqsample[filterCondition]

# Check reads of filteredSequences
sread(filteredSequences)

#Three
"""
How many reads with a maximum of 3 consecutive A's per read did you find?

""" 
# 3

##################################################

#Plotting cycle average quality
#One
# Load package Rqc
library(Rqc)

# Average per cycle quality plot
rqcCycleAverageQualityPlot(qa)

#Two
# Average per cycle quality plot with white background
rqcCycleAverageQualityPlot(qa) + theme_minimal()

#Three
# Read quality plot with white background
rqcReadQualityPlot(qa) + theme_minimal()


