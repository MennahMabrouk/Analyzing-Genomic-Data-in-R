#Chapter One
# Bioconductor version
# Load the BiocManager package
library(BiocManager)

# Explicitly check the Bioconductor version
version()


#################################################
#BiocManager to install packages
# Load GenomicRanges
library(GenomicRanges)
# Check versions for reproducibility
library(BSgenome)
sessionInfo()


#################################################
"""
S4 class definition
We will use the class BSgenome, which is already loaded for you.

Let's check the formal definition of this class by using 
the function showClass("className"). 
Check the BSgenome class results and find its parent classes (Extends) 
and the classes that inherit from it (Subclasses).

"""
#Annotated and MaskedBSgenome


#################################################
#Interaction with classes
# Investigate the a_genome using show()
.S4methods(class = "GenomeDescription")

showMethods(classes = "GenomeDescription", where = search())

show(a_genome)

#Two
# Investigate the a_genome using show()
show(a_genome)

# Investigate some other accesors
organism(a_genome)
provider(a_genome)
seqinfo(a_genome)

#################################################


#Discovering the yeast genome
#One
# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

#Two
# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the head of seqnames and tail of seqlengths for yeastGenome
head(seqnames(yeastGenome))
tail(seqlengths(yeastGenome))

#Three
# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the head of seqnames and tail of seqlengths for yeastGenome
head(seqnames(yeastGenome))
tail(seqlengths(yeastGenome))

# Print chromosome M, alias chrM
print(yeastGenome$chrM)

#################################################



#Partitioning the yeast genome
# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the first 30 bases of chrM
getSeq(yeastGenome,name="chrM",end=30)


#################################################



"""
Available genomes
As a recap, the BSgenome package contains various public genomes. 
If you want to explore the available genomes from this package, you can use:

available.genomes()
The list of names will appear in the following format: 
BSgenome.speciesName.provider.version.

After running this function, which of the following 
is a major provider of available genomes?

"""
#UCSC








