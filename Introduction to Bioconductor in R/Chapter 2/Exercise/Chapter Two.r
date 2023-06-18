#Chapter Two
#One
#Exploring the Zika virus sequence
# Load packages
library(Biostrings)

# Check the alphabet of the zikaVirus
alphabet(zikaVirus)

#Two
# Load packages
library(Biostrings)

# Check the alphabet of the zikaVirus
alphabet(zikaVirus)

# Check the alphabetFrequency of the zikaVirus
alphabetFrequency(zikaVirus)

#Three
# Load packages
library(Biostrings)

# Check the alphabet of the zikaVirus
alphabet(zikaVirus)

# Check the alphabetFrequency of the zikaVirus
alphabetFrequency(zikaVirus)

# Check alphabet of the zikaVirus using baseOnly = TRUE
alphabet(zikaVirus, baseOnly = TRUE)

#########################

"""
Biostrings containers
Now that you know how to check the alphabet of a sequence, 
can you check what the container class of the zikaVirus object is?

"""
#DNAStringset

##############################################

#Manipulating Biostrings
#One
# Unlist the set, select the first 21 letters, and assign to dna_seq
dna_seq <- subseq(unlist(zikaVirus), end = 21)
dna_seq

# Transcribe dna_seq into an RNAString object and print it
rna_seq <- RNAString(dna_seq) 
rna_seq

#Two
# Unlist the set, select the first 21 letters, and assign to dna_seq
dna_seq <- subseq(unlist(zikaVirus), end = 21)
dna_seq

# Transcribe dna_seq into an RNAString object and print it
rna_seq <- RNAString(dna_seq) 
rna_seq

# Translate rna_seq into an AAString object and print it
aa_seq <- translate(rna_seq)
aa_seq

#Three
# Unlist the set, select the first 21 letters, and assign to dna_seq
dna_seq <- subseq(unlist(zikaVirus), end = 21)
dna_seq

# Transcribe and translate dna_seq into an AAString object and print it
aa_seq <- translate(dna_seq)
aa_seq

###################################################

#From a set to a single sequence
#One
# Create zikv with one collated sequence using `zikaVirus`
zikv <- unlist(zikaVirus)

# Check the length of zikaVirus and zikv
length(zikaVirus)
length(zikv)

#Two
# Create zikv with one collated sequence using zikaVirus
zikv <- unlist(zikaVirus)

# Check the length of zikaVirus and zikv
length(zikaVirus)
length(zikv)

# Check the width of zikaVirus
width(zikaVirus)

#Three
# Create zikv with one collated sequence using zikaVirus
zikv <- unlist(zikaVirus)

# Check the length of zikaVirus and zikv
length(zikaVirus)
length(zikv)

# Check the width of zikaVirus
width(zikaVirus)

# Subset zikv to only the first 30 bases
subZikv <- subseq(zikv, end = 30)
subZikv

#################################################

"""
Subsetting a set
In the previous exercise, you used subseq() to subset a single sequence.
Here, you can try subseq() using a set with 3 sequences. 
The arguments are the same as before: object, start, and end. 
The last two should be in the form of a vector, 
for each of the sequences on the set.

subseq(zikaSet, 
        start = c(20, 40, 2), 
        end = c(50, 45, 22)
     )
What is the width of the sequences after subsetting?

"""
#31,6,21

###############################################
#Common sequence manipulation functions
# Reverse the zikv sequence
reverse(zikv)

# Complement the zikv sequence
complement(zikv)

# Reverse complement the zikv sequence
reverseComplement(zikv)

# Translate the zikv sequence
translate(zikv)

###############################################

"""
Searching for a pattern
Let's find some occurrences of a pattern in the zikaVirus set using
vmatchPattern(). Then, let's try the same pattern search using matchPattern() 
with a single sequence, zikv.

# For Sets
vmatchPattern(pattern = "ACATGGGCCTACCATGGGAG", 
              subject = zikaVirus, max.mismatch = 1)
# For single sequences
matchPattern(pattern = "ACATGGGCCTACCATGGGAG", 
              subject = zikv, max.mismatch = 1)
Both functions should find the same number of occurrences, 
but you will notice a different output. 
How many matches do we get when running each pattern search individually?

"""
# 1 match

####################################

#Finding Palindromes
# Find palindromes in zikv
findPalindromes(zikv)

######################################
#Finding a conserved region within six frames
#One
# Print rnaframesZikaSet 
print(rnaframesZikaSet)

# Translate rnaframesZikaSet 
AAzika6F <- translate(rnaframesZikaSet)
AAzika6F

#Two
# Print rnaframesZikaSet
rnaframesZikaSet

# Translate rnaframesZikaSet
AAzika6F <- translate(rnaframesZikaSet)
AAzika6F

# Count NS5 protein matches in AAzika6F, allowing 15 mismatches
vcountPattern(pattern = NS5, subject = AAzika6F, max.mismatch = 15)

#Three
# Print rnaframesZikaSet
rnaframesZikaSet

# Translate rnaframesZikaSet
AAzika6F <- translate(rnaframesZikaSet)
AAzika6F

# Count NS5 protein matches in AAzika6F, allowing 15 mismatches
vcountPattern(pattern = NS5, subject = AAzika6F, max.mismatch = 15)

# Subset the frame that contains the match from AAzika6F
selectedSet <- AAzika6F[unlist(vcountPattern(pattern = NS5, subject = AAzika6F, max.mismatch = 15) > 0)] 

# Convert selectedSet into a single sequence
selectedSeq <- unlist(selectedSet)

####################################################################

#Looking for a match
#One
# Use vmatchPattern() with the set
vmatchPattern(pattern = ns5, subject = selectedSet, max.mismatch = 15)

#Two
# Use matchPattern() with the single sequence
matchPattern(pattern= ns5, subject=selectedSeq, max.mismatch = 15)

# Take your time to see the similarities/differences in the result.



