# Keerthana Vishwanath
# Must use command line prompt to call program
# Use a sample .pdb (protein database) file when entering arguments
import math, re, sys


# EuclideanDistance()
# Input: two lists of equal size containing floats
# Output: a float giving the distance between two multidimensional points
def EuclideanDistance(vec1, vec2):
    total = 0
    for i, j in zip(vec1, vec2):
        total += (i - j) ** 2
    return math.sqrt(total)


# A simple class for storing and retrieving data associated with atoms in the .pdb
class ATOM:
    def __init__(self, chainID, residue_number, x, y, z):
        self.chainID = chainID
        self.residue_number = residue_number
        self.x = x
        self.y = y
        self.z = z

    def coords(self):
        return [self.x, self.y, self.z]


# LoadPDB()
# Input: path to .pdb file containing docking frames
# Output: list of dictionaries, each dictionary mapping chain names (A, B, etc) to lists of alpha-carbon ATOM instances
def LoadPDB(filepath):
    models = []
    with open(filepath, "r") as fh:
        chains = {}

        for line in fh:
            line = re.split("\s+", line)

            # if a given line beginning with MODEL, we are at a new model
            # consequently, we should dump the current dictionary into the models list and create an empty working dictionary
            if line[0] == "MODEL":
                models.append({key: chains[key] for key in chains.keys()})
                chains = {}

            # if given an ATOM line representing an alpha-carbon, we need to store information in the line
            # for legibility's sake we store the data in variables, add a key to the chains dictionary if the chain is not already present,
            # and store the newly-instanced ATOM into the appropriate list in the chains dictionary
            elif line[0] == "ATOM" and line[2] == "CA":
                chainID = line[4]
                residue_number = int(line[5]) - 1
                x = float(line[6])
                y = float(line[7])
                z = float(line[8])

                if chainID not in chains.keys():
                    chains[chainID] = []

                chains[chainID].append(ATOM(chainID, residue_number, x, y, z))

            else:
                continue

    del models[0]
    return models


# GenerateIPCM()
# Input: two lists of ATOMs, a threshold value (in angstroms)
# Output: a 2D list containing an intermolecular protein contact map
def GenerateIPCM(chainA, chainB, threshold):
    IPCM_matrix = [[0 for i in range(len(chainA))] for j in range(len(chainB))]

    for i in range(len(chainA)):
        for j in range(len(chainB)):
            if EuclideanDistance((chainA[i].coords()), (chainB[j].coords())) < threshold:
                IPCM_matrix[i][j] = 1

    return IPCM_matrix


# GenerateConsensusIPCM()
# Input: a list of 2D lists, each 2D list containing an intermolecular protein contact map
# Output: a 2D list containing a consensus intermolecular protein contact map
def GenerateConsensusIPCM(maps):
    IPCM = maps[0]
    numRow = len(IPCM)
    numCol = len(IPCM[0])
    numIPCM = len(maps)  # should be 15 using example.pdb file

    print(numRow)
    print(numCol)
    print(numIPCM)
    consensusIPCM = [[0 for x in range(numRow)] for y in range(numCol)]
    for IPCM in maps:
        for row in range(numRow):
            for col in range(numCol):
                consensusIPCM[row][col] += IPCM[row][col]  # creation of the consensus IPCM
    for row in range(numRow):
        for col in range(numCol):
            consensusIPCM[row][col] = consensusIPCM[row][col] / float(
                numIPCM)  # modification of the consensus IPCM to have final values
    return consensusIPCM


# WritePCM()
# Input: a filepath to write to, a 2D list
# writes a 2D list to a file in the form of a matrix
def WritePCM(filepath, map):
    fh = open(filepath, "w")
    for i in range(len(map)):
        temp = ""
        for j in range(len(map[i])):
            temp += str(map[i][j]) + "\t"
        fh.write(temp.strip() + "\n")
    fh.close()


models = LoadPDB(sys.argv[1])

IPCMs = []
for i in range(len(models)):

    print("Generating IPCM for model " + str(i))

    IPCMs.append(GenerateIPCM(models[i]["A"], models[i]["B"], 5.0))

    WritePCM("model" + str(i) + ".txt", IPCMs[-1])

consensus = GenerateConsensusIPCM(IPCMs)

WritePCM("consensus.txt", consensus)
