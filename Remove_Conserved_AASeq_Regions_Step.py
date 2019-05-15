#Keerthana Vishwanath
# How to find and removes conserved regions of amino acid sequences (this is the 3rd step in the overall process)
# Our group's work looks at the HA and NA segments across all Human Influenza A strains from 1968 to 2019
# All sequences passed through this program must be aligned beforehand
import sys
import os

def load_fasta(data, seqDict, seqCountDict): # loads file into dictionary
    prevSeqKey = ""
    sequence = ""

    for line in data:
        if line[0:1] == ">":
            seqKey = line[0:len(line) - 1]  # current sequence key

            if prevSeqKey != "":
                seqDict[prevSeqKey] = sequence
                try:
                    seqCountDict[prevSeqKey] += 1
                except:
                    seqCountDict[prevSeqKey] = 1
                sequence = ""
            prevSeqKey = seqKey

        else:
            sequence += line[0:len(line)-1]
    if sequence != "":
        seqDict[prevSeqKey] = sequence
        try:
            seqCountDict[prevSeqKey] += 1
        except:
            seqCountDict[prevSeqKey] = 1
    print(len(seqDict))
    return seqDict, seqCountDict

seqDict = {}  # dictionary of UNIQUE sequences
seqCountDict = {}  # tracks number of repeating sequences (duplicates)

# takes in a directory with all relevant fasta files (should ONLY contain fasta files)
fastafiles = []
inputdir = raw_input("Enter a fasta directory name: ")

print ("Input fasta directory name is: ", inputdir)

fh = None
data = None
try:
    fastafiles = os.listdir(inputdir)  # gives files in the directory
    for fastafile in fastafiles:
        fh = open(inputdir + "/" + fastafile, "r")
        data = fh.readlines()
        load_fasta(data, seqDict, seqCountDict)  # loads all files from directory into dictionary
except:
    print "Directory Not Found"
finally:
    if fh != None:
        fh.close()
if data == None:
    sys.exit("No data in directory.")
outCountFile = open("count.out", "w")
for seqKey, count in seqCountDict.items():
    outCountFile.write(seqKey + "--" + str(count))
    outCountFile.write("\n")
outCountFile.close()

sequence = ""
i = 0
j = 0
prevSeqChar = ""
curSeqChar = ""
modDict = {}
processDict = {}
print(len(seqDict))

# read the first character of the first value
# compare the character of the first value across all other values
outDebugFile = open("debug.out", "w") # checking to see whether the modified fasta file has the same number of sequences
for seqKey, sequence in seqDict.items():
    outDebugFile.write(seqKey + "--" + str(len(sequence)))
    outDebugFile.write("\n")
outDebugFile.close()

curSeqLen = 0
seqLen = 0
for seqKey, sequence in seqDict.items():  # walks through dictionary using both key and value
    curSeqLen = len(sequence)  # determines longest sequence
    if curSeqLen > seqLen:  # current sequence length is longer than the existing sequence length
        seqLen = curSeqLen  # sets sequence length to current sequence length--length of longest sequence
for j in range(seqLen - 1):  # takes max sequence length as the range, loops through using sequence character position (j)
    charMatches = True
    for seqKey, sequence in seqDict.items(): # for each character position--goes through entire sequence dictionary

        if j < len(sequence): # process ONLY if character position is less than the length of sequence, otherwise move onto next sequence
            curSeqChar = sequence[j]
            if prevSeqChar == "":
                prevSeqChar = curSeqChar
                processDict[seqKey] = curSeqChar
                continue
            if curSeqChar != prevSeqChar:
                charMatches = False

            processDict[seqKey] = curSeqChar
    prevSeqChar = ""
    # if character in position across all sequences does NOT match, add character into dictionary for that particular sequence, else drop the character from all sequences
    if charMatches != True:

        for processKey, processSeq in processDict.items():
            try:
                modDict[processKey] += processSeq
            except:
                modDict[processKey] = processSeq

outFile = open("modifiedSequence.fasta", "w") #writes modified sequences to fasta file
# WILL overwrite the existing file when re-run with different data
for seqKey, sequence in modDict.items():
    outFile.write(seqKey + "\n")
    outFile.write(sequence + "\n")

outFile.close()









