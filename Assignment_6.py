#Keerthana Vishwanath
#BNFO 301: Introduction to Bioinformatics
#Assignment #6: Divides fasta file into unique k-mers and analyzes them--each analysis includes the following:
# 1. Length of the sequence
# 2. The observed nucleotide frequencies (as percentages) of the sequence
# 3. The total number of unique k-mers observed and how many could have been observed
# 4. Sort k-mers by number of observations and by chi-square value  and report in a table
    # The table consists of the following: The sequence of a k-mer, the observed count of that k-mer in the sequence
    # the expected count of that k-mer in the sequence, based on the observed nucleotide frequencies (Rounded to nearest integer)
    # a chi-square value comparing the observed count to the expected count (chi-sqare = (observed - expected)^2/expected
import sys

# class to calculate and store observed kmer values
# allows for sorting on observed values (Frequency)
class KMERDIST:
    def __init__(self, kmer, frequency, probability, expectedVal, chiSqVal):
        self.kmer = kmer
        self.frequency = frequency
        self.probability = probability
        self.expectedVal = expectedVal
        self.chiSqVal = chiSqVal
    def __getitem__(self, key):
        return key

    def __lt__(self, input):
        return input.frequency < self.frequency

# class to calculate and store probability, expected values, and chi-squared values
# allows for sorting on chi-squared
class CHISQKMER:
    def __init__(self, kmer, frequency, probability, expectedVal, chiSqVal):
        self.kmer = kmer
        self.frequency = frequency
        self.probability = probability
        self.expectedVal = expectedVal
        self.chiSqVal = chiSqVal

    def __getitem__(self, key):
        return key

    def __lt__(self, input):
        return input.chiSqVal < self.chiSqVal

#accepts fasta file and k-mer length from user
inputfile = input("Enter a fasta file name: ")
print("Input file name is: ", inputfile)
inputkmerlength = input("Enter a k-mer length: ")
print("Input k-mer length is: ", inputkmerlength)
fh = None
kmerlength = int(inputkmerlength)
kmers = []
data = None
try:
    fh = open(inputfile, "r")
    data = fh.readlines()
except:
    print("File Not Found")
finally:
    if fh != None:
        fh.close()
if data == None:
    sys.exit("No data in file.")
#print(data)
sequence = ""
i = 0
# Each line in the data list represents a line in the input fasta file
# Each line has a newline character at the end.
# Extract all the characters in each line except the newline
# Concatenate the extracted characters to sequence string
for i in range(1,(len(data) -1)):
    sequence += data[i][0:len(data[i]) -1]


sequenceLength = len(sequence) #length of sequence

# Calculate the percent of each nucleotide in the sequence
print("\nSequence length: " + str(sequenceLength))
nucApercentage = round(float(float(sequence.count("A"))/sequenceLength),3)
nucTpercentage = round(float(float(sequence.count("T"))/sequenceLength),3)
nucCpercentage = round(float(float(sequence.count("C"))/sequenceLength),3)
nucGpercentage = round(float(float(sequence.count("G"))/sequenceLength),3)
print("\nNucleotide Frequencies: \n")
print("A         " + str(nucApercentage))
print("C         " + str(nucCpercentage))
print("G         " + str(nucGpercentage))
print("T         " + str(nucTpercentage))
print("\n")

#Separating each of the sequences into the k-mer length chosen via user input
for i in range((sequenceLength-1)-kmerlength):
    kmers.append(sequence[i:i+kmerlength])



kmers.sort() #sorts the kmers


#for kmerstring in sortedkmers:
#previous variable, current variable (if current is same as previous then increment the counter
#if current is not same as previous then populate the class and add it to the list and
# reset the counter to one and initialize previous kmer with current kmer
#if the for loop terminates the last kmer and its count should added to the list
previousKmer = ""
kcounter = 0
currentKmer = ""
kmerdistribution = []
for currentKmer in kmers:
    if currentKmer == previousKmer or len(previousKmer) == 0:
        previousKmer = currentKmer
        kcounter += 1
    else:
        kmerdistribution.append(KMERDIST(previousKmer, kcounter, 0.0, 0, 0.0))
        previousKmer = currentKmer
        kcounter = 1

if kcounter != 0:
    kmerdistribution.append(KMERDIST(previousKmer, kcounter, 0.0, 0, 0.0))
sumFrequency = 0
for distribution in kmerdistribution:
    sumFrequency += distribution.frequency

uniqueKmers = set()
for kmer in kmers:
    uniqueKmers.add(kmer)
print(str(len(uniqueKmers)) + " of " + str(len(kmerdistribution)) + " possible " + inputkmerlength + "-mers found.\n")

#Sorted the kmerdistribution on observed values
kmerdistribution.sort()

kmerChiSqDistribution = []
for kmer in kmerdistribution:
    probKmer = 1.0
    probKmerA = 1.0
    probKmerT = 1.0
    probKmerC = 1.0
    probKmerG = 1.0
    for nuc in kmer.kmer:
        if nuc == "A":
            probKmerA *= nucApercentage
        if nuc == "T":
            probKmerT *= nucTpercentage
        if nuc == "C":
            probKmerC *= nucCpercentage
        if nuc == "G":
            probKmerG *= nucGpercentage
    probKmer = probKmerA * probKmerT * probKmerC * probKmerG
    kmerChiSqDistribution.append(CHISQKMER(kmer.kmer, kmer.frequency, probKmer, round(probKmer * sequenceLength),
                                    round((kmer.frequency - round(probKmer * sequenceLength,4))**2/round(probKmer * sequenceLength,4),1)))

#stores the top 10 kmers
topTenKmers = []

count = 0
for kmer in kmerChiSqDistribution:
    if count < 10:
        topTenKmers.append(kmer)
        count += 1
    else:
        break


print("\nTop 10 Highest Frequency\n")
print("Seq	        Obs	        Exp	        Chi-Square")
print("--------------------------------------------------")
for kmer in topTenKmers:
    print(kmer.kmer + "       " + str(kmer.frequency) + "         " + str(kmer.expectedVal) + "        " + str(kmer.chiSqVal))

#sorts the kmers by chi-squared values
kmerChiSqDistribution.sort()
print("\nTop 10 Highest Chi-Square\n")
print("Seq	        Obs	        Exp	        Chi-Square")
print("--------------------------------------------------")
count = 0
for kmer in kmerChiSqDistribution:
    if count < 10:
        print(kmer.kmer + "       " + str(kmer.frequency) + "         "  + str(kmer.expectedVal) + "        " + str(kmer.chiSqVal))
        count += 1
    else:
        break