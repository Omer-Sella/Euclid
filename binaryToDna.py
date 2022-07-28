"""
Created on Wed Jul 22 09:20:04 2020

@author: Omer Sella
"""
import numpy as np
import re
import matplotlib.pyplot as plt
import math
import copy
import matplotlib.cm as cm
import matplotlib.colors
#import xlsxwriter

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
        bitStream = (bin(int(inputByte.hex(), 16))[2:]).zfill(numOfBits)
    else:
        bitStream = inputByte
    fourNuceleotides = ""
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
            fourNuceleotides = fourNuceleotides + nuceleotide
        else:
            fourNuceleotides = 'bad'
    return fourNuceleotides

def binaryToDna(inputFilename, lineLength, outputFilename = None, format = 'IDT'):
    
    if lineLength % 4 != 0:
        actualLineLenghth = lineLength + (4 - lineLength % 4)
    else:
        actualLineLenghth = lineLength
    
    
    
    with open(outputFilename, "w") as outFile:
        with open(inputFilename, "rb") as inFile:
            sequenceNumber = 0
            if format == 'IDT':
                line = "seq" + str(sequenceNumber) + "\t"
            else:
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
                    print(lineCounter)
                    print(len(line))
                    line = line + "\n"
                    outFile.write(line)
                    sequenceNumber = sequenceNumber + 1
                    if format == 'IDT':
                        line = "seq" + str(sequenceNumber) + "\t"
                    else:
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


def produceStatsFromSequence(sequence, userSubSequences = None):
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
    
    for s in range(len(subSequences)):
        allMatches = re.findall(subSequences[s], sequence)
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
    singlesNormed = 2 * np.pi * singles / np.sum(singles)
    singlesLeft = np.cumsum(np.append(0, singlesNormed[:-1]))
    cmap = plt.get_cmap("tab20c")
    singlesColors = cmap(np.arange(4)*4)
    if len(occurances) > 4 :
        doubles = copy.copy(occurances[4:20])
        doublesNormed = np.zeros(len(doubles), dtype = np.float32)
        doublesLeft = np.zeros(len(doubles), dtype = np.float32)    
        doublesColors = cmap(np.arange(16))
        for i in range(4):
            relevantSum = np.sum(doubles[4 * i : 4 * (i + 1)])
            vector = singlesNormed[i] * doubles[4 * i : 4 * (i + 1)] / relevantSum 
            doublesNormed[4 * i : 4 * (i + 1)] = vector #singlesNormed[i] * doublesNormed[4 * i : 4 * (i + 1)] / relevantSum
            doublesLeft[4 * i : 4 * (i + 1)] = np.cumsum(np.append(0, doublesNormed[4 * i : 4 * (i + 1) - 1] )) + singlesLeft[i]
        
        tripletsColors = cmap(np.arange(4)*4)
        if len(occurances) > 20 :
            triplets = copy.copy(occurances[20:])
            tripletsNormed = np.zeros(len(triplets), dtype = np.float32)
            tripletsLeft = np.zeros(len(triplets), dtype = np.float32)
            for i in range(16):
                relevantSum = np.sum(triplets[4 * i : 4 * (i + 1)])
                vector = doublesNormed[i] * triplets[4 * i : 4 * (i + 1)] / relevantSum 
                tripletsNormed[4 * i : 4 * (i + 1)] = vector #singlesNormed[i] * doublesNormed[4 * i : 4 * (i + 1)] / relevantSum
                tripletsLeft[4 * i : 4 * (i + 1)] = np.cumsum(np.append(0, tripletsNormed[4 * i : 4 * (i + 1) - 1] )) + doublesLeft[i]
    
        
    
    
    # for i in range(4):
    #     relevantSum = np.sum(doubles[4 * i : 4 * (i + 1)])
    #     #print(relevantSum)
    #     #print(doubles[4 * i : 4 * (i + 1)])
    #     #print(singlesNormed[i])
    #     vector = singlesNormed[i] * doubles[4 * i : 4 * (i + 1)] / relevantSum 
    #     #print(vector)
    #     doublesNormed[4 * i : 4 * (i + 1)] = vector #singlesNormed[i] * doublesNormed[4 * i : 4 * (i + 1)] / relevantSum
    #     #print(doublesNormed[4 * i : 4 * (i + 1)])
    #     doublesLeft[4 * i : 4 * (i + 1)] = np.cumsum(np.append(0, doublesNormed[4 * i : 4 * (i + 1) - 1] )) + singlesLeft[i]
    #     #print(doublesLeft[4 * i : 4 * (i + 1)])
    
    

    
    
    
    
    s1 = ax.bar(x=singlesLeft,
           width=singlesNormed, bottom = 0, height = size,
           color = singlesColors, edgecolor='w', linewidth=1, align="edge", label = 'Single base')

    for s in range(len(s1)):
        rect = s1[s]
        height = rect.get_height()
        ax.annotate('{}'.format(subSequences[s]),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 0),  
                    textcoords="offset points",
                    ha='center', va='center', size = 36)
    if len(occurances) > 4 :
        s2 = ax.bar(x=doublesLeft,
           width=doublesNormed, bottom = 0.33, height = size,
           color = doublesColors, edgecolor='w', linewidth=1, align="edge", label = 'Two bases')
        for s in range(len(s2)):
            rect = s2[s]
            height = rect.get_height()
            ax.annotate('{}'.format(subSequences[4 + s]),
                        xy=(rect.get_x() + rect.get_width() / 2, 2 * height),
                        xytext=(0, 0),  
                        textcoords="offset points",
                        ha='center', va='center', size = 22)
    if len(occurances) > 20:
        s3 = ax.bar(x=tripletsLeft,
           width=tripletsNormed, bottom = 0.66, height = size,
           color = tripletsColors, edgecolor='w', linewidth=1, align="edge", label = 'Three bases')
        for s in range(len(s3)):
            rect = s3[s]
            height = rect.get_height()
            ax.annotate('{}'.format(subSequences[4 + 16 + s]),
                        xy=(rect.get_x() + rect.get_width() / 2, 3 * height),
                        xytext=(0, 0),  
                        textcoords="offset points",
                        ha='center', va='center', size = 20)

    ax.set_title("Distribution of subsequences.", fontsize = 36)
    ax.set_axis_off()
    #ax.legend() Doesn't work well. gives a blue rectangle.
    plt.show()
    




def checkViolation(symbol1, symbol2, constraintList):
    
    violations = 0
    gcBalanceViolation = False
    #reversedSymbol2 = symbol2[::-1]
    #reversedSymbol1 = symbol1[::-1]
    #concatenatedSymbols = reversedSymbol1 + symbol2
    concatenatedSymbols = symbol1 + symbol2
    #print(concatenatedSymbols)
    if ('runLength' in constraintList):
        runL = constraintList['runLength']
        for l in ['A', 'C', 'T', 'G']:
            constraintList['regexRunL' + l] = l * runL
    
    for i in constraintList:
        if i[0:5] == 'regex':
            subSequence = constraintList[i]
            violations = violations + len(re.findall(subSequence, concatenatedSymbols))
    
    if ('gcMin' in constraintList):
        gcMin = constraintList['gcMin']
    else:
        gcMin = 0.25
        
    if ('gcMax' in constraintList):
        gcMax = constraintList['gcMax']
    else:
        gcMax = 0.65
    
    gcCount = len(re.findall('G', concatenatedSymbols)) + len(re.findall('C', concatenatedSymbols))
    if ((gcCount / len(concatenatedSymbols)) > gcMax) or ( (gcCount / len(concatenatedSymbols))  < gcMin):
        gcBalanceViolation = True
        
    if (violations > 0) or (gcBalanceViolation == True):
        noViolation = 0
    else:
        noViolation = 1
    
    return noViolation, violations, gcBalanceViolation
        
def plotConnectivityMatrix(connectivityMatrix, xLabels = None, yLabels = None, xSize = None, ySize = None, fileName = None, figureSize = None, figParams = None, colour = None, alpha = None, verticalSpacing = None, horizontalSpacing = None):
    
    if xSize is None:
        xFontSize = 28
    else:
        xFontSize = xSize
    
    if ySize is None:
        yFontSize = 28
    else:
        yFontSize = ySize
    
    (verticalDimension, horizontalDimension) = connectivityMatrix.shape
    if figParams is None:
        fig, ax = plt.subplots()
    else:
        fig = figParams[0]
        ax = figParams[1]
        plt.sca(ax)
        
    if colour is None:
        colourMap = 'Greys'
    else:
        colourMap = colour
        
    
    ax.imshow((-1 * connectivityMatrix) + 1, cmap=colourMap,  interpolation = None, alpha = alpha)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    
    if verticalSpacing is None:
        spacingVertical = 5 #verticalDimension // 50
    else:
        spacingVertical = verticalSpacing
    if horizontalSpacing is None:
        spacingHorizontal = 5 #horizontalDimension // 50
    else:
        spacingHorizontal = horizontalSpacing
    
    if xLabels is not None :
        xTickLocations = np.arange(0, horizontalDimension, spacingHorizontal)
        xTickValues = []
        for i in xTickLocations:
            xTickValues.append(xLabels[i])
        plt.xticks(xTickLocations, xTickValues, fontsize = xFontSize, rotation = 90)
    
    if yLabels is not None :
        yTickLocations = np.arange(0, verticalDimension, spacingVertical)
        yTickValues = []
        for i in yTickLocations:
            yTickValues.append(yLabels[i])
    
        plt.yticks(yTickLocations, yTickValues, fontsize = yFontSize)
        
    if fileName is not None:
        if figureSize is not None:
            plt.savefig(fileName, figsize = figureSize)
        else:
            plt.tight_layout()
            plt.savefig(fileName)
            
    return fig, ax
    
    

def connectivityMatrix(symbolSize = 4, consraintList = {}):
    
    sequencesBinary = []
    sequencesBases = []
    k = 0
    symbol = "0" * (2 * symbolSize)
    sequencesBinary.append(symbol)
    sequencesBases.append(byteToBases(symbol))
    
    while (k < 2 ** (2 * symbolSize) - 1):
        symbol = binaryPlusOne(symbol)
        sequencesBinary.append(symbol)
        sequencesBases.append(byteToBases(symbol))
        k = k + 1
    
    #for number, value in enumerate(sequencesBinary):
    #    print(str(number) + "  " + str(value))
    connectivityDictionaryBinary = {}
    connectivityMatrix = np.zeros((2**(symbolSize * 2), 2**(symbolSize * 2)), dtype = np.int8)
    violationsMatrix = np.zeros((2**(symbolSize * 2), 2**(symbolSize * 2)), dtype = np.int8)
    for i in range(2**(symbolSize * 2)):
        connectivityDictionaryBinary[sequencesBinary[i]] = {}
        for j in range(2**(symbolSize * 2)):
            connectivityMatrix[i,j], violations, gcBalance = checkViolation(sequencesBases[i], sequencesBases[j], consraintList)
            violationsMatrix[i,j] = violations + gcBalance
            connectivityDictionaryBinary[sequencesBinary[i]][sequencesBinary[j]] = connectivityMatrix[i,j]
            
    # fig, ax = plt.subplots()
    # ax.imshow((-1 * connectivityMatrix) + 1, cmap='Greys',  interpolation = None)
    # ax.spines['top'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    
    # tickLocations = np.arange(0,256,5)
    # tickValues = []
    # for i in tickLocations:
    #     tickValues.append(sequencesBases[i])
    # plt.yticks(tickLocations, tickValues, fontsize = 16)
    # plt.xticks(tickLocations, tickValues, fontsize = 16, rotation = 90)
    # figV, axV = plt.subplots()
    # im = axV.imshow(violationsMatrix, cmap="YlGn")
    # cbar = axV.figure.colorbar(im)
    return connectivityMatrix, sequencesBinary, sequencesBases, violationsMatrix, connectivityDictionaryBinary
        


    
def capacityDictatedByMatrix(matrix, symbolSize):
    rightEigenValues, rightEigenVectors = np.linalg.eig(matrix)
    maximumPower = 7
    matrix = np.uint64(matrix)
    m,n = matrix.shape
    wann = np.ones((m,1), dtype = np.uint64)
    pathsArray = np.zeros(maximumPower + 1, dtype = np.uint64)
    capacity = np.zeros(maximumPower + 1)
    poweredMatrix = np.eye(m) #copy.copy(matrix) #np.eye(m)
    for i in range(maximumPower):
        poweredMatrix = poweredMatrix.dot(matrix)
        pathsArray[i + 1] = np.sum(poweredMatrix.dot(wann))
        print(pathsArray[i + 1])
        capacity[i + 1] =  (1 / (i + 2)) * math.log(pathsArray[i + 1], 4**symbolSize)
            
    bitCount = pathsArray * symbolSize
    return capacity, bitCount, pathsArray

def constrainedChannelGraphics(matrix, symbolSize):
    rightEigenValues, rightEigenVectors = np.linalg.eig(matrix)
    lambda1 = np.abs(rightEigenValues[0])
    print("***First right eigenvalue:" + str(lambda1))
    print("***Log of first right eigenvalue:" + str(math.log(lambda1, 4**symbolSize)))
    fig, ax = plt.subplots()
    capacity, bitCount, pathsArray = capacityDictatedByMatrix(matrix, symbolSize)
    nAxis = np.arange(1,8,1)
    ax.plot(nAxis, capacity[1:], label = "Mackay capacity", color = "green")
    ax.hlines(math.log(lambda1, 4**symbolSize), nAxis[0], nAxis[1], label = 'Log of first eigenvalue', color = 'blue')
    
    #ax.plot(pathsArray)
    #ax.set_yscale('log')
    ax.legend()
    return fig, ax, capacity[-1]

def generateXsequence(symbolToAdd, sequencesBases):
    
    xSequence = []
    for seq in sequencesBases:
        newSeq = symbolToAdd + seq[1 : ]
        xSequence.append(newSeq)
    
    return xSequence
# def matrixReduce(matrix, axis):
    
#     symbolSize = 4
#     sequencesBinary = []
#     sequencesBases = []
#     k = 0
#     symbol = "0" * (2 * symbolSize)
#     sequencesBinary.append(symbol)
#     sequencesBases.append(byteToBases(symbol))
    
#     while (k < 2 ** (2 * symbolSize) - 1):
#         symbol = binaryPlusOne(symbol)
#         sequencesBinary.append(symbol)
#         sequencesBases.append(byteToBases(symbol))
#         k = k + 1
    
#     #for number, value in enumerate(sequencesBinary):
#     #    print(str(number) + "  " + str(value))
    
    
    
#     tickLocations = np.arange(0,4**symbolSize,5)
#     tickValues = []
#     for i in tickLocations:
#         tickValues.append(sequencesBases[i])
#     #plt.yticks(tickLocations, tickValues, fontsize = 16)
#     #plt.xticks(tickLocations, tickValues, fontsize = 16, rotation = 90)
    
    
    
#     assert (axis == 0 or axis == 1)
#     (m,n) = matrix.shape
    
#     if axis == 0:   
#         newMatrix1 = matrix[0 : m // 2, :]
#         newMatrix2 = matrix[m // 2 : , :]
#     else:
#         newMatrix1 = matrix[:, 0 : n // 2]
#         newMatrix2 = matrix[:, n // 2 : ]
    
#     assert (newMatrix1.shape == newMatrix2.shape)
#     foldedMatrix = np.where(newMatrix1 < newMatrix2, newMatrix1, newMatrix2)
#     if np.all(foldedMatrix == 0):
#         print("*** All zeros")
#     return foldedMatrix, newMatrix1, newMatrix2


"""    
gcAxis = np.array([0.15, 0.35])#np.arange(0.15, 0.35, 0.1)
capacity = np.zeros(gcAxis.shape)
logEigenValues = np.zeros(gcAxis.shape)
symbolSize = 4

i = 0
for v in gcAxis:
    constraintList = {"runLength": 8, 'gcMin': 0.5 - v, 'gcMax':0.5 + v}
    A, seq, V = connectivityMatrix(4, constraintList)
    rightEigenValues, rightEigenVectors = np.linalg.eig(A)
    lambda1 = np.abs(rightEigenValues[0])
    logEigenValues[i] = math.log(lambda1, 4**symbolSize)
    fig, ax, capacity[i] = constrainedChannelGraphics(A, 4)
    i = i + 1
    
    

fig, ax = plt.subplots()
ax.plot(gcAxis, capacity, label = "Mackay capacity", color = "green", marker = 'D', linewidth = 3)
ax.plot(gcAxis, logEigenValues, label = "Base 256 log of first eigenvalue", color = "pink", marker = 'p', linewidth = 3)
ax.legend(fontsize = 18)
ax.tick_params(axis = 'x', labelsize = 18)
ax.tick_params(axis = 'y', labelsize = 18)
ax.set_title('Capacity as a function of G-C balance', fontsize = 18)
ax.set_xlabel('G-C variance allowed, centred around 0.5', fontsize = 18)
ax.set_ylabel('Capacity', fontsize = 18)
"""


gcContentDemoList = {'gcMin': 0.15, 'gcMax': 0.85, 'runLength': 10}#, 'regex1': '[ACTG][ACTG][ACTG]AA[ACTG][ACTG][ACTG]'}
gcContentDemoList2 = {'gcMin': 0.35, 'gcMax': 0.65, 'runLength': 10}#, 'regex1': '[ACTG][ACTG][ACTG]AA[ACTG][ACTG][ACTG]'}
#demoConstraintList = {'gcMin': 0.25, 'gcMax': 0.65, 'runLength': 7}#, 'regex1': '[ACTG][ACTG][ACTG]AA[ACTG][ACTG][ACTG]'}
demoConstraintList = {'gcMin': 0.0, 'gcMax': 1.0, 'runLength': 10, 
                      'regex1': '[ACTG][ACTG][ACTG]AA[ACTG][ACTG][ACTG]',
                      'regex2': '[ACTG][ACTG][ACTG]CC[ACTG][ACTG][ACTG]',
                      'regex3': '[ACTG][ACTG][ACTG]GA[ACTG][ACTG][ACTG]',
                      'regex4': '[ACTG][ACTG][ACTG]CT[ACTG][ACTG][ACTG]'}#,
                      #'regex5': '[ACTG][ACTG][ACTG]CG[ACTG][ACTG][ACTG]'}

#constraintListThatDeniesAAinTheMiddle = {'gcMin': 0.0, 'gcMax': 1.0, 'runLength': 9, 'regex1': '[ACTG][ACTG][ACTG]AA[ACTG][ACTG][ACTG]'}

def produceFigure4():
    #demoConstraintList = gcContentDemoList2#gcContentDemoList
    demoConstraintList = {'gcMin': 0.15, 'gcMax': 0.85, 'runLength': 6, 'regex1': 'AGCCTAG', 'regex2': 'TCGGATC', 'regex3': 'GATCCGA', 'regex4': 'CTAGGCT','regex5': 'AGCAGC', 'regex6': 'A[AGCT][AGCT][AGCT]CTAG'}
    A, seqBinary, seqBases, V, _ = connectivityMatrix(4, demoConstraintList)

    plotConnectivityMatrix(A, seqBases, seqBases, xSize = 28, ySize = 28)

    (m,n) = A.shape
    ALeft = A[:, 0 : n//2]
    ARight = A[:, n//2 : n]
    AUp = A[0 : n//2, :]
    ADown = A[n//2 : n, :]
    
    plotConnectivityMatrix(ALeft, seqBases[0 : n//2], seqBases)
    plotConnectivityMatrix(ARight, seqBases[n//2 : n], seqBases)
    
    
    #plotConnectivityMatrix(AUp, seqBases, seqBases[0 : n//2])
    #plotConnectivityMatrix(ADown, seqBases, seqBases[n//2 : n])
    
    ySequence = generateXsequence('Ry', seqBases[0 : n//2])
    xSequence = generateXsequence('Rx', seqBases[0 : n//2])
    
    reducedALR = np.where(ALeft > ARight, ALeft, ARight)
    reducedAUD = np.where(ADown > AUp, ADown, AUp)
    #finalReduced = np.where(reducedAUD[:, 0 : n//2] > reducedAUD[:, n//2 : n], reducedAUD[:, 0 : n//2], )
    #plotConnectivityMatrix(reducedAUD, seqBases, xSequence)
    plotConnectivityMatrix(reducedALR, xLabels = xSequence, yLabels = seqBases)
    
    reducedALR_left =  reducedALR[: , 0 : n // 4]
    reducedALR_right = reducedALR[: , n // 4 : n // 2]
    xxSequence = generateXsequence('RR', seqBases[0 : n//4])
    plotConnectivityMatrix(reducedALR_left, yLabels = seqBases, xLabels = xxSequence)
    plotConnectivityMatrix(reducedALR_right, yLabels = seqBases, xLabels = xxSequence)
    
    # Omer Sella: Very
    reducedReducedALR = np.where(reducedALR_left > reducedALR_right, reducedALR_left, reducedALR_right)
    plotConnectivityMatrix(reducedReducedALR, yLabels = seqBases, xLabels = xxSequence)
    
    
    
    reducedReducedALR_left =  reducedReducedALR[: , 0 : n // 8]
    reducedReducedALR_right = reducedReducedALR[: , n // 8 : n // 4]
    xxxSequence = generateXsequence('RR', seqBases[0 : n//4])
    plotConnectivityMatrix(reducedReducedALR_left, yLabels = seqBases)
    plotConnectivityMatrix(reducedReducedALR_right, yLabels = seqBases)
    xRRSequence = generateXsequence('RRx', seqBases[0 : n//4])
    # Omer Sella: Very
    RRRSequence = generateXsequence('RRRx', seqBases[0 : n // 8])
    reducedReducedReducedALR = np.where(reducedReducedALR_left > reducedReducedALR_right, reducedReducedALR_left, reducedReducedALR_right)
    plotConnectivityMatrix(reducedReducedReducedALR, yLabels = seqBases, xLabels = RRRSequence)
    
    fig, ax = plt.subplots(1, 4, gridspec_kw={'width_ratios': [2, 2, 2, 1]})
    plotConnectivityMatrix(ALeft, xLabels = seqBases[0 : n//2], yLabels = seqBases, figParams = [fig, ax[0]])
    plotConnectivityMatrix(ARight, xLabels = seqBases[n//2 : ], yLabels = seqBases, figParams = [fig, ax[1]])
    plotConnectivityMatrix(reducedALR, xLabels = xSequence, yLabels = seqBases, figParams = [fig, ax[2]])
    plotConnectivityMatrix(reducedReducedALR, xLabels = xxSequence, yLabels = seqBases, figParams = [fig, ax[3]])
    ax[0].set_title('(a)', size = 40, pad = 15)
    ax[1].set_title('(b)', size = 40, pad = 15)
    ax[2].set_title('(c)', size = 40, pad = 15)
    ax[3].set_title('(d)', size = 40, pad = 15)
    return xxxSequence


def reduceMatrix(inMatrix):
    _, l = inMatrix.shape
    left = inMatrix[:, 0: l//2 ]
    right = inMatrix[:, l//2 : l]
    outMatrix = np.where(left > right, left, right)
    fullyConnected = np.all(outMatrix)
    return outMatrix, fullyConnected

def findBitValue(matrix, maxBits):
    i = 0
    fullyConnected = np.all(matrix)
    while (i < maxBits) and not fullyConnected:
        outMatrix, fullyConnected = reduceMatrix(matrix)
        matrix = outMatrix
        i = i + 1
    return i, fullyConnected

def scanGCValues():
    low = np.arange(0.05, 0.5, 0.05)
    high = np.arange(0.55, 0.95, 0.05)
    results = []
    for l in low:
        for h in high:
            gcContentConstraint = {'gcMin': l, 'gcMax': h, 'runLength': 10}#, 'regex1': '[ACTG][ACTG][ACTG]AA[ACTG][ACTG][ACTG]'}
            matrix, seqBinary, seqBases, V = connectivityMatrix(4, gcContentConstraint)
            numOfBits, isDone = findBitValue(matrix, 7)
            print(numOfBits)
            print(isDone)
            if isDone:
                results.append([l,h,numOfBits])
    return results

def getPareto():
    colors = cm.rainbow(np.linspace(0, 1, 8))
    results = scanGCValues()
    fig, ax = plt.subplots()
    markerSize = 128
    for g in range(8):
        lows = []
        highs = []
        for [l,h,b] in results:
            if b == g:
                lows.append(l)
                highs.append(h)
        if len(lows) > 0:
            text = str(g) + " reserved bits"
            ax.scatter(lows, highs, c=colors[g], label = text, s = np.ones(len(lows))*markerSize)
    ax.set_title("Pareto fronts of G-C content limits as a function of reserved bits", fontsize = 18)
    ax.set_xlabel("G-C content low limit", fontsize = 18)
    ax.set_ylabel("G-C content high limit", fontsize = 18)
    ax.set_ylim(0.92, 0.53)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)
    ax.legend(fontsize = 18)
    
    return fig, ax

def multipleConstraintsGraphics():
    fig, ax = plt.subplots(1, 4, gridspec_kw={'width_ratios': [1, 1, 1, 1]})
    
    
    
    homopolymerList = {'gcMin': 0.0, 'gcMax': 1.0, 'runLength': 6}#, 'regex1': '[ACTG][ACTG][ACTG]AA[ACTG][ACTG][ACTG]'}
    homopolymerA, seqBinary, seqBases, V, _ = connectivityMatrix(4, homopolymerList)
    homopolymerCM = matplotlib.colors.ListedColormap(['white', 'red'])
    
    primerList = {'gcMin': 0.0, 'gcMax': 1.0, 'runLength': 10, 'regex1': 'AGCCTAG', 'regex2': 'TCGGATC', 'regex3': 'GATCCGA', 'regex4': 'CTAGGCT','regex5': 'AGCAGC', 'regex6': 'A[AGCT][AGCT][AGCT]CTAG'}
    primerA, _, _, _ = connectivityMatrix(4, primerList)
    primerCM = matplotlib.colors.ListedColormap(['white', 'steelblue'])
    
    gcList = {'gcMin': 0.15, 'gcMax': 0.85, 'runLength': 10}#, 'regex1': '[ACTG][ACTG][ACTG]AA[ACTG][ACTG][ACTG]'}
    gcA, _, _, _ = connectivityMatrix(4, gcList)
    gcCM = matplotlib.colors.ListedColormap(['white', 'slateblue'])
    
    A = np.minimum(np.minimum(homopolymerA, primerA), gcA)
    
    
    
    plotConnectivityMatrix(homopolymerA, xLabels = seqBases, yLabels = seqBases, figParams = [fig, ax[0]], colour = homopolymerCM, xSize = 20, ySize = 20, verticalSpacing = 7, horizontalSpacing = 7)
    plotConnectivityMatrix(primerA, xLabels = seqBases, yLabels = seqBases, figParams = [fig, ax[1]], colour = primerCM, xSize = 20, ySize = 20, verticalSpacing = 7, horizontalSpacing = 7)
    plotConnectivityMatrix(gcA, xLabels = seqBases, yLabels = seqBases, figParams = [fig, ax[2]], colour = gcCM, xSize = 20, ySize = 20, verticalSpacing = 7, horizontalSpacing = 7)
    plotConnectivityMatrix(A, xLabels = seqBases, yLabels = seqBases, figParams = [fig, ax[3]], xSize = 20, ySize = 20, verticalSpacing = 7, horizontalSpacing = 7)
    
    plotConnectivityMatrix(homopolymerA, xLabels = seqBases, yLabels = seqBases, colour = homopolymerCM, xSize = 28, ySize = 28)
    plotConnectivityMatrix(primerA, xLabels = seqBases, yLabels = seqBases, colour = primerCM, xSize = 28, ySize = 28)
    plotConnectivityMatrix(gcA, xLabels = seqBases, yLabels = seqBases, colour = gcCM, xSize = 28, ySize = 28)
    plotConnectivityMatrix(A, xLabels = seqBases, yLabels = seqBases, xSize = 28, ySize = 28)
    
    #plotConnectivityMatrix(gcA, xLabels = seqBases, yLabels = seqBases, figParams = [fig, ax[3]], colour = gcCM, alpha = 0.9)
    #plotConnectivityMatrix(primerA, xLabels = seqBases, yLabels = seqBases, figParams = [fig, ax[3]], colour = primerCM, alpha = 0.9)
    #plotConnectivityMatrix(homopolymerA, xLabels = seqBases, yLabels = seqBases, figParams = [fig, ax[3]], colour = homopolymerCM, alpha = 0.9)
    ax[0].set_title('(a)', size = 40, pad = 15)
    ax[1].set_title('(b)', size = 40, pad = 15)
    ax[2].set_title('(c)', size = 40, pad = 15)
    ax[3].set_title('(d)', size = 40, pad = 15)
    
    return A, fig, ax