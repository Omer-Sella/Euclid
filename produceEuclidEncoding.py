# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 10:10:30 2022

@author: omers
      
"""
# import os
# import sys
# projectDir = os.environ.get('EUCLID')
# if projectDir == None:
#     projectDir = "E:/Euclid"
# sys.path.insert(1, projectDir)
from binaryToDna import reduceMatrix, connectivityMatrix, plotConnectivityMatrix
import copy
from itertools import combinations_with_replacement
import numpy as np


# def isFullyConnected(matrix):
#     result = (np.sum(matrix) == matrix.size)
#     return result

# def reduce(matrix):
#     (m,n) = matrix.shape
#     matrixLeft = matrix[:, 0 : n//2]
#     matrixRight = matrix[:, n//2 : n]
#     reducedLeftRight = np.where(matrixLeft > matrixRight, matrixLeft, matrixRight)
#A, seqBinary, seqBases, V = connectivityMatrix(4, demoConstraintList)



def generateCandidates(inMatrix, connectivityDictionary, verticalSymbols, horizontalSymbols, alphabet):
    numberOfReservedSymbols = 0
    # I'm assuming here that inMatrix is not all 1s, otherwise there is no point in using this encoding.
    fullyConnected = False
    reducedMatrix = copy.deepcopy(inMatrix)
    while not fullyConnected and (np.mod(reducedMatrix.shape[1], 2) == 0):
        reducedMatrix, fullyConnected = reduceMatrix(reducedMatrix)
        numberOfReservedSymbols = numberOfReservedSymbols + 1
    print("Number of reserved bits is: " + str(numberOfReservedSymbols))
    
    #Generate pool of suffixes:

    prefixes = copy.deepcopy(alphabet)
    i = 1
    while i < numberOfReservedSymbols:
        i = i + 1
        newPrefixes = []
        for prefix in prefixes:
            for aleph in alphabet:
                newPrefix = copy.deepcopy(prefix + aleph)
                newPrefixes.append(newPrefix)
        prefixes = newPrefixes
        
    #Generate unique suffixes, as they appear in the original list
    newHorizontalSymbols = []
    for h in horizontalSymbols:
        suffixH = h[numberOfReservedSymbols :]
        if not suffixH in newHorizontalSymbols:
            newHorizontalSymbols.append(suffixH)
    
    candidatesDictionary = {}
    for verticalSymbol in verticalSymbols:
        candidatesDictionary[verticalSymbol] = {}
        for horizontalSymbol in newHorizontalSymbols:
            tempList = []
            for prefix in prefixes:
                tempSymbol = prefix + horizontalSymbol
                #print(verticalSymbol)
                #print(horizontalSymbol)
                if connectivityDictionary[verticalSymbol][tempSymbol] == 1:
                    tempList.append(prefix)
            candidatesDictionary[verticalSymbol][horizontalSymbol] = tempList
             
        
    return candidatesDictionary, verticalSymbols, newHorizontalSymbols

def euclidCandidatesForTheJoin():
    restrictionSite = 'GAGTC'
    startPrimer = 'GAACCGTGCCGAGTCTGAGC'
    sp0 = 'GAACCGTG'
    sp1 = 'AACCGTGC'
    sp2 = 'ACCGTGCC'
    sp3 = 'CCGTGCCG'
    sp4 = 'CGTGCCGA'
    sp5 = 'GTGCCGAG'
    sp6 = 'TGCCGAGT'
    sp7 = 'GCCGAGTC'
    sp8 = 'CCGAGTCT'
    sp9 = 'CGAGTCTG'
    sp10 = 'GAGTCTGA'
    sp11 = 'AGTCTGAG'
    sp12 = 'GTCTGAGC'
    
    endPrimer = 'CTCAGGACTCGCAACGCTGG'
    ep0 = 'CTCAGGAC'
    ep1 = 'TCAGGACT'
    ep2 = 'CAGGACTC'
    ep3 = 'AGGACTCG'
    ep4 = 'GGACTCGC'
    ep5 = 'GACTCGCA'
    ep6 = 'ACTCGCAA'
    ep7 = 'CTCGCAAC'
    ep8 = 'TCGCAACG'
    ep9 = 'CGCAACGC'
    ep10 ='GCAACGCT'
    ep11 ='CAACGCTG'
    ep12 ='AACGCTGG'
    
    PARTSSUP_TABLE_ID_PRIMER = 'GCGACTGGATGACCTGACGC'
    
    
    partssupPrimer0 = 'GCGACTGG'
    partssupPrimer1 = 'CGACTGGA'
    partssupPrimer2 = 'GACTGGAT'
    partssupPrimer3 = 'ACTGGATG'
    partssupPrimer4 = 'CTGGATGA'
    partssupPrimer5 = 'TGGATGAC'
    partssupPrimer6 = 'GACCTGAC'
    partssupPrimer7 = 'ACCTGACG'
    partssupPrimer8 = 'CCTGACGC'
    
    
    
    PARTS_TABLE_ID_PRIMER = 'GCAGACCGGAGACCTGTCGG'
    
    partsPrimer0 = 'GCAGACCG'
    partsPrimer1 = 'CAGACCGG'
    partsPrimer2 = 'AGACCGGA'
    partsPrimer3 = 'GACCGGAG'
    partsPrimer4 = 'ACCGGAGA'
    partsPrimer5 = 'CCGGAGAC'
    partsPrimer6 = 'CGGAGACC'
    partsPrimer7 = 'GGAGACCT'
    partsPrimer8 = 'GAGACCTG'
    partsPrimer9 = 'AGACCTGT'
    partsPrimer10 = 'GACCTGTC'
    partsPrimer11 = 'ACCTGTCG'
    partsPrimer12 = 'CCTGTCGG'
    
    
    joinConstraintList = {'gcMin': 0.25, 'gcMax': 0.65, 'runLength': 8, 
                      'regex1': 'GAGTC',
                      # Restriction site GAGTC is specifc, i.e.: the enzyme recognizes it exactly.
                      #'regex2': '[ACTG]AGTC',
                      #'regex3': 'G[ACTG]GTC',
                      #'regex4': 'GA[ACTG]TC',
                      #'regex5': 'GAG[ACTG]C',
                      #'regex6': 'GAGT[ACTG]',
                      'regex7': sp0,
                      'regex8': sp1,
                      'regex9': sp2,
                      'regex10': sp3,
                      'regex11': sp4,
                      'regex12': sp5,
                      'regex13': sp6,
                      'regex14': sp7,
                      'regex15': sp8,
                      'regex16': sp9,
                      'regex17': sp10,
                      'regex18': sp11,
                      'regex19': sp12,
                      
                      'regex20': ep0,
                      'regex21': ep1,
                      'regex22': ep2,
                      'regex23': ep3,
                      'regex24': ep4,
                      'regex25': ep5,
                      'regex26': ep6,
                      'regex27': ep7,
                      'regex28': ep8,
                      'regex29': ep9,
                      'regex30': ep10,
                      'regex31': ep11,
                      'regex32': ep12,
                      
                      'regex33': partssupPrimer0,
                      'regex34': partssupPrimer1,
                      'regex35': partssupPrimer2,
                      'regex36': partssupPrimer3,
                      'regex37': partssupPrimer4,
                      'regex38': partssupPrimer5,
                      'regex39': partssupPrimer6,
                      'regex40': partssupPrimer7,
                      'regex41': partssupPrimer8,
                      
                      'regex42': partsPrimer0,
                      'regex43': partsPrimer1,
                      'regex44': partsPrimer2,
                      'regex45': partsPrimer3,
                      'regex46': partsPrimer4,
                      'regex47': partsPrimer5,
                      'regex48': partsPrimer6,
                      'regex49': partsPrimer7,
                      'regex50': partsPrimer8,
                      'regex51': partsPrimer9,
                      'regex52': partsPrimer10,
                      'regex53': partsPrimer11,
                      'regex54': partsPrimer12
                      }#,
    
    matrix, seqBinary, seqBases, V, connectivityDictionary = connectivityMatrix(4, joinConstraintList)

    plotConnectivityMatrix(matrix, seqBases, seqBases, xSize = 28, ySize = 28)
    
    candidates, vSymbols, hSymbols = generateCandidates(matrix, connectivityDictionary, seqBinary, seqBinary, ['0', '1'])
    
    return candidates, connectivityDictionary, vSymbols, hSymbols

def getStats(candidates, verticalSymbols, horizontalSymbols):
    countMatrix = np.zeros((len(verticalSymbols), len(horizontalSymbols)) )
    for key in candidates.keys():
        i = verticalSymbols.index(key)
        for subKey in candidates[key].keys():
            j = horizontalSymbols.index(subKey)
            countMatrix[i,j] = len(candidates[key][subKey])
    return countMatrix

def makeFSM(candidates, verticalSymbols, horizontalSymbols, mechanism):
    #An FSM is made of:
    #states, triggers, outputTable, transitionTable, initialState
    countMatrix = getStats(candidates, verticalSymbols, horizontalSymbols)
    states = verticalSymbols
    triggers = horizontalSymbols
    # Now comes the hard part: we need to choose an output (validPrefix + newHorizontalSymbol) for every (state, trigger) pair.
    outputDictionary = {}
    for vs in verticalSymbols:
        outputDictionary[vs] = {}
        for hs in horizontalSymbols:
            outputDictionary[vs][hs] = mechanism(vs, candidates[vs][hs], hs)
            
    return countMatrix, outputDictionary


def minimiseReservedValue(a,candidatesList,c):
    candidatesAsInt = np.zeros(len(candidatesList))
    for i in range(len(candidatesList)):
        candidatesAsInt[i] = int(candidatesList[i], 2)
    return candidatesList[(np.argmin(candidatesAsInt))]

def minimiseLeftRightMagnitude(a, candidatesList, c):
    aInt = int(a, 2)
    candidatesAndCAsInt = np.zeros(len(candidatesList))
    for i in range(len(candidatesList)):
        candidatesAndCAsInt[i] = int((candidatesList[i] + c), 2)
    difference = np.abs(candidatesAndCAsInt - aInt)
    return candidatesList[np.argmin(difference)]

def testEuclid():
    c, cd, vsym, hsym = euclidCandidatesForTheJoin()
    numberOfPossibleCandidatesCountMatrix, outputDictionary = makeFSM(c, vsym, hsym, minimiseReservedValue)
    