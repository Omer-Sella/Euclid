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
                tempSymbol = prefix + suffixH
                #print(verticalSymbol)
                #print(horizontalSymbol)
                if connectivityDictionary[verticalSymbol][tempSymbol] == 1:
                    tempList.append(prefix)
            candidatesDictionary[verticalSymbol][suffixH] = tempList
             
        
    return candidatesDictionary

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
    
    joinConstraintList = {'gcMin': 0.2, 'gcMax': 0.80, 'runLength': 8, 
                      'regex1': 'GAGTC',
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
                      'regex19': sp12}#,
    matrix, seqBinary, seqBases, V, connectivityDictionary = connectivityMatrix(4, joinConstraintList)

    plotConnectivityMatrix(matrix, seqBases, seqBases, xSize = 28, ySize = 28)
    
    candidates = generateCandidates(matrix, connectivityDictionary, seqBinary, seqBinary, ['0', '1'])
    
    return candidates, connectivityDictionary

def getStats(candidates):
    count = []
    for key in candidates:
        #There is only one subkey
        for subKey in candidates[key]:
            count.append(len(candidates[key][subKey]))
    return count

def makeFSM(candidates):
    #An FSM is made of:
    #states, triggers, outputTable, transitionTable, initialState
    count = getStats(candidates)
    states = []
    for key in candidates:
        states.append(key)