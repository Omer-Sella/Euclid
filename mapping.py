
import numpy as np
import re
import matplotlib.pyplot as plt
import math
import copy
import matplotlib.cm as cm
import matplotlib.colors

def binaryPlusOne(operand):
    """
    Adds 1 to the binary string operand.
    No safety,i.e.: operand is assumed to be a string, with only binary content 9i.e.: 0s or 1s),
            however a string may be returened only if it has the same length as operand.
    """
    
    #print(str(operand))
    length = len(operand)
    k = len(operand) - 1
    newSymbol = operand
    notDone = True
    
    while (k >= 0 ) and notDone:
        if operand[k] == '0':
            
            notDone = False
            
            if k == 0:
                newSymbol = "1" + newSymbol[1:]
            if k == length - 1:
                newSymbol = newSymbol[0 : k ] + "1"
            else:
                newSymbol = newSymbol[0 : k ] + "1" + newSymbol[k + 1 : ]
        else:
            if k == 0:
                newSymbol = "0" + newSymbol[1:]
            if k == length - 1:
                newSymbol = newSymbol[0 : k] + "0"
            else:
                newSymbol = newSymbol[0 : k] + "0" + newSymbol[k + 1:]
            k = k - 1
        assert(len(newSymbol) == length)
    return newSymbol

def byteToBases(inputByte, numOfBits = 8):
    if type(inputByte) != str:
        print("Input stream is not string ... ")
        bitStream = (bin(int(inputByte.hex(), 16))[2:]).zfill(numOfBits)
    else:
        bitStream = inputByte
    nuceleotideStream = ""
    for i in range(numOfBits // 2):
        twoBits = bitStream[2 * i : 2 * i + 2]
        #print(twoBits)
        if twoBits == "00":
            nuceleotide = 'A'
        elif twoBits == "01":
            nuceleotide = 'C'
        elif twoBits == "11":
            nuceleotide = 'T'
        elif twoBits == "10":
            nuceleotide = 'G'
        else:
            nuceleotide = 'bad'
        if nuceleotide != 'bad':
            nuceleotideStream = nuceleotideStream + nuceleotide
        else:
            nuceleotideStream = 'bad'
        
    return nuceleotideStream

 
def byte2binary(byte, numOfBits = 8):
    bitStream = (bin(int(byte.hex(), 16))[2:]).zfill(numOfBits)
    return bitStream


def binaryStreamToBases(bitStream, numOfBits = 8):
    nuceleotideStream = ""
    if type(bitStream) != str:
        bitStream = (bin(int(bitStream.hex(), 16))[2:]).zfill(numOfBits)
    numOfBits = len(bitStream)
    for i in range(numOfBits // 2):
        twoBits = bitStream[2 * i : 2 * i + 2]
        #print(twoBits)
        if twoBits == "00":
            nuceleotide = 'A'
        elif twoBits == "01":
            nuceleotide = 'C'
        elif twoBits == "11":
            nuceleotide = 'T'
        elif twoBits == "10":
            nuceleotide = 'G'
        else:
            nuceleotide = 'bad'
        if nuceleotide != 'bad':
            nuceleotideStream = nuceleotideStream + nuceleotide
        else:
            nuceleotideStream = 'bad'
                
    return nuceleotideStream

def binaryToDna(inputFilename, lineLength, outputFilename = None):
    if lineLength % 4 != 0:
        actualLineLenghth = lineLength + (4 - lineLength % 4)
    else:
        actualLineLenghth = lineLength
    with open(outputFilename, "w") as outFile:
        with open(inputFilename, "rb") as inFile:
            line = ""
            lineCounter = 0
            byte = inFile.read(1)
            while byte != b"":
                newBases = byteToBases(byte)
                if newBases == 'bad':
                    print('***Stopped due to badness ...')
                    break
                line = line + newBases
                
                lineCounter = lineCounter + 4
                if lineCounter == actualLineLenghth:
                    line = line + "\n"
                    outFile.write(line)
                    line = ""
                    lineCounter = 0
                # Read new byte
                byte = inFile.read(1)

def produceStatsFromDnaFile(inputFilename, userSubSequences = None):
    basePairs = ['A','T','C','G']
    subSequences = []
    for c0 in basePairs:
        subSequences = subSequences + [c0]
    for c0 in basePairs:
        for c1 in basePairs:
            subSequences = subSequences + [c0 + c1]
    for c0 in basePairs:
        for c1 in basePairs:
            for c2 in basePairs:
                subSequences = subSequences + [c0 + c1 + c2]
    if userSubSequences != None:
        subSequences = subSequences + userSubSequences
    occurences = np.zeros(len(subSequences), dtype = np.int32)
    print(subSequences)
    with open(inputFilename, "r") as inFile:
        lines = inFile.readlines()
        for line in lines:
            for s in range(len(subSequences)):
                allMatches = re.findall(subSequences[s], line)
                occurences[s] = occurences[s] + len(allMatches)
    print(occurences)
    return subSequences, occurences


def produceStatsFromDnaStream(inputStream, userSubSequences = None):
    basePairs = ['A','T','C','G']
    subSequences = []
    for c0 in basePairs:
        subSequences = subSequences + [c0]
    for c0 in basePairs:
        for c1 in basePairs:
            subSequences = subSequences + [c0 + c1]
    for c0 in basePairs:
        for c1 in basePairs:
            for c2 in basePairs:
                subSequences = subSequences + [c0 + c1 + c2]
    if userSubSequences != None:
        subSequences = subSequences + userSubSequences
    occurences = np.zeros(len(subSequences), dtype = np.int32)
    print(subSequences)
    for s in range(len(subSequences)):
        allMatches = re.findall(subSequences[s], inputStream)
        occurences[s] = occurences[s] + len(allMatches)
    print(occurences)
    return subSequences, occurences

"""
sub = ['A', 'T', 'C', 'G', 'AA', 'AT', 'AC', 'AG', 'TA', 'TT', 'TC', 'TG', 'CA', 'CT', 'CC', 'CG', 'GA', 'GT', 'GC', 'GG', 'AAA', 'AAT', 'AAC', 'AAG', 'ATA', 'ATT', 'ATC', 'ATG', 'ACA', 'ACT', 'ACC', 'ACG', 'AGA', 'AGT', 'AGC', 'AGG', 'TAA', 'TAT', 'TAC', 'TAG', 'TTA', 'TTT', 'TTC', 'TTG', 'TCA', 'TCT', 'TCC', 'TCG', 'TGA', 'TGT', 'TGC', 'TGG', 'CAA', 'CAT', 'CAC', 'CAG', 'CTA', 'CTT', 'CTC', 'CTG', 'CCA', 'CCT', 'CCC', 'CCG', 'CGA', 'CGT', 'CGC', 'CGG', 'GAA', 'GAT', 'GAC', 'GAG', 'GTA', 'GTT', 'GTC', 'GTG', 'GCA', 'GCT', 'GCC', 'GCG', 'GGA', 'GGT', 'GGC', 'GGG']
occ = np.array([1774547, 1629353, 1351267, 1628193,  377500,  419691,  415585,
        431407,  445577,  320004,  361464,  390256,  402150,  348748,
        199148,  360649,  418888,  428938,  334330,  335052,  107800,
        125560,  104975,  101667,   94901,   94823,  109722,  111912,
        106795,   89565,   86488,  125447,  116640,  113617,   73143,
        119172,  111723,   99509,   99673,  126178,  116411,   93890,
         71386,  101906,  131523,   79841,   50381,   92929,  107735,
         87763,   80803,  106567,   95742,  104736,   88050,  107268,
         89757,   97302,   70098,   84150,   68114,   75476,   29666,
         49989,   92005,  107704,   72979,   81066,  124795,   81565,
        116339,   87755,  136078,   97586,  102755,   84946,   88571,
         97048,   56711,   84111,   93774,  112292,  100609,   92962])
"""


def subSequenceGraphics(subSequences, occurances):
    #print(occurances)
    fig, ax = plt.subplots(subplot_kw=dict(polar=True))
    numberOfCircles = 3
    size = 0.3
    singles = copy.copy(occurances[0:4])
    doubles = copy.copy(occurances[4:20])
    triplets = copy.copy(occurances[20:])
    singlesNormed = 2 * np.pi * singles / np.sum(singles)
    singlesLeft = np.cumsum(np.append(0, singlesNormed[:-1]))
    doublesNormed = np.zeros(len(doubles), dtype = np.float32)
    doublesLeft = np.zeros(len(doubles), dtype = np.float32)
    tripletsNormed = np.zeros(len(triplets), dtype = np.float32)
    tripletsLeft = np.zeros(len(triplets), dtype = np.float32)
    #print(singlesNormed)
    #print(singlesLeft)
    
    for i in range(4):
        relevantSum = np.sum(doubles[4 * i : 4 * (i + 1)])
        #print(relevantSum)
        #print(doubles[4 * i : 4 * (i + 1)])
        #print(singlesNormed[i])
        #BUG
        #Cheasy bug fix - in case relevantSum is 0
        vector = singlesNormed[i] * doubles[4 * i : 4 * (i + 1)] / (1 + relevantSum )
        #print(vector)
        doublesNormed[4 * i : 4 * (i + 1)] = vector #singlesNormed[i] * doublesNormed[4 * i : 4 * (i + 1)] / relevantSum
        #print(doublesNormed[4 * i : 4 * (i + 1)])
        doublesLeft[4 * i : 4 * (i + 1)] = np.cumsum(np.append(0, doublesNormed[4 * i : 4 * (i + 1) - 1] )) + singlesLeft[i]
        #print(doublesLeft[4 * i : 4 * (i + 1)])
    
    for i in range(16):
        relevantSum = np.sum(triplets[4 * i : 4 * (i + 1)])
        print(relevantSum)
        #print(doubles[4 * i : 4 * (i + 1)])
        #print(singlesNormed[i])
        #BUG
        #Cheasy bug fix - in case relevantSum is 0
        vector = doublesNormed[i] * triplets[4 * i : 4 * (i + 1)] / (1 + relevantSum ) 
        #print(vector)
        tripletsNormed[4 * i : 4 * (i + 1)] = vector #singlesNormed[i] * doublesNormed[4 * i : 4 * (i + 1)] / relevantSum
        #print(doublesNormed[4 * i : 4 * (i + 1)])
        tripletsLeft[4 * i : 4 * (i + 1)] = np.cumsum(np.append(0, tripletsNormed[4 * i : 4 * (i + 1) - 1] )) + doublesLeft[i]
        #print(doublesLeft[4 * i : 4 * (i + 1)])
    cmap = plt.get_cmap("tab20c")
    singlesColors = cmap(np.arange(4)*4)
    doublesColors = cmap(np.arange(16))
    tripletsColors = cmap(np.arange(4)*4)
    s1 = ax.bar(x=singlesLeft,
           width=singlesNormed, bottom = 0, height = size,
           color = singlesColors, edgecolor='w', linewidth=1, align="edge", label = 'Single base')
    s2 = ax.bar(x=doublesLeft,
           width=doublesNormed, bottom = 0.33, height = size,
           color = doublesColors, edgecolor='w', linewidth=1, align="edge", label = 'Two bases')
    s3 = ax.bar(x=tripletsLeft,
           width=tripletsNormed, bottom = 0.66, height = size,
           color = tripletsColors, edgecolor='w', linewidth=1, align="edge", label = 'Three bases')

    for s in range(len(s1)):
        rect = s1[s]
        height = rect.get_height()
        ax.annotate('{}'.format(subSequences[s]),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 0),  
                    textcoords="offset points",
                    ha='center', va='center', size = 24)
    for s in range(len(s2)):
        rect = s2[s]
        height = rect.get_height()
        ax.annotate('{}'.format(subSequences[4 + s]),
                    xy=(rect.get_x() + rect.get_width() / 2, 2 * height),
                    xytext=(0, 0),  
                    textcoords="offset points",
                    ha='center', va='center', size = 20)
    for s in range(len(s3)):
        rect = s3[s]
        height = rect.get_height()
        ax.annotate('{}'.format(subSequences[4 + 16 + s]),
                    xy=(rect.get_x() + rect.get_width() / 2, 3 * height),
                    xytext=(0, 0),  
                    textcoords="offset points",
                    ha='center', va='center', size = 18)

    ax.set_title("Distribution of subsequences.", fontsize = 18)
    ax.set_axis_off()
    plt.show()


def dnaToBinary(inputStream, outputFile = None):
    """
    sample input ACTGAAACCTGA
    Convert DNA text into binary.
    input is assumed to be a list object so input[i] is a DNA strand, and input[i][j] is the j'th nucleotide of strand i.
    """
    status = 'OK'
    newBinarySequence = ''
    newBinaryArray = np.zeros(2 * len(inputStream))
    idx = 0
    for nucleotide in inputStream:
        if nucleotide == 'A':
            bitPairText = '00'
            bitPairNumerical = np.array([0,0])
        elif nucleotide == 'C':
            bitPairText = '01'
            bitPairNumerical = np.array([0,1])
        elif nucleotide == 'T':
            bitPairText = '11'
            bitPairNumerical = np.array([1,1])
        elif nucleotide == 'G':
            bitPairText = '10'
            bitPairNumerical = np.array([1,0])
        elif nucleotide == '\n':
            bitPairText = '\n'
            bitPairNumerical = np.array([2,2])
        else:
            bitPairText = 'BAD_NUCLEOTIDE_ENCOUNTERED'
            bitPairNumerical = np.array([2,2])
            status = 'FAIL'
        newBinarySequence = newBinarySequence + bitPairText
        newBinaryArray[idx : idx + 2] = bitPairNumerical
        idx = idx + 2
        if outputFile != None:
            with open(outputFile, "a+") as fid:
                fid.write(newBinarySequence + "\n")
    return status, newBinarySequence, newBinaryArray