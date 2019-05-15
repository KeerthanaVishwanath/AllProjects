#Keerthana Vishwanath

import sys, re


def load_fasta(filepath):

    fh = open(filepath, "r")
    data = fh.readlines()
    fh.close() #reads in file
    label = ""
    fullSequence = ""
    firstTime = True

    for x in data: #looks at each line in file
        if x[0] == '>': #checks for label beginning marker
            if firstTime == True: #boolean created to store label marker on first iteration
                label = x
                firstTime = False
            else:
                checkSign = does_regex_match(fullSequence) #runs does_regex_match function on the currently stored sequence
                if checkSign == True: #if the regex exists within sequence, prints sequence label
                    print(label[1:])
                fullSequence = "" #resets fullSequence
                label = x #gets label of next sequence
        else:
            fullSequence = fullSequence + x #appends each line of sequence to a string
    checkSign = does_regex_match(fullSequence) #runs does_regex_match function on the last sequence
    if checkSign == True:
        print(label[1:]) #if the regex exists within sequence, prints sequence label
    print("Regex pattern for Alcohol Dehydrogenase Activity is: GHE[ARNDCQEGHILKMFPOSUTWYVBZXJ][^EL]G[^AP][ARNDCQEGHILKMFPOSUTWYVBZXJ]{4}[GA][ARNDCQEGHILKMFPOSUTWYVBZXJ]{2}[IVSAC]")
    #prints the Regex pattern used to find the Alcohol Dehydrogenase Activity


def does_regex_match(sequence):
    #compiles regex
    pattern = re.compile('GHE[ARNDCQEGHILKMFPOSUTWYVBZXJ][^EL]G[^AP][ARNDCQEGHILKMFPOSUTWYVBZXJ]{4}[GA][ARNDCQEGHILKMFPOSUTWYVBZXJ]{2}[IVSAC]')

    matches = re.findall(pattern, sequence) #finds all occurrences of pattern in sequence
    if len(matches) > 0:
        return True #if pattern exists--returns true
    else:
        return False #if pattern does not exist--returns false




if __name__ == '__main__':
    filepath = sys.argv[1] #gets filepath name from commandline call
    load_fasta(filepath) #runs load_fasta method to read data from file



    pass