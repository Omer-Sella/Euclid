# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 10:10:30 2022

@author: omers
      
"""
from binaryToDNA import reduceMatrix
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

def generateCandidates(inMatrix, inList, verticalSymbols, horizontalSymbols, alphabet):
    numberOfReservedSymbols = 0
    fullyConnected = False
    reducedMatrix = copy.copy(inMatrix)
    while not fullyConnected and (np.mod(reducedMatrix.shape[1], 2) == 0):
        reducedMatrix, fullyConnected = reduceMatrix(reducedMatrix)
        numberOfReservedSymbols = numberOfReservedSymbols + 1
    
    #Generate pool of suffixes:
    combinationsIter = combinations_with_replacement(alphabet, numberOfReservedSymbols)
    prefixes = []
    combinationsList = list(combinationsIter)
    for c in combinationsList:
        newPrefix = "".join(c)
        prefixes.append(newPrefix)
    
    candidates = {}
    newHorizontalSymbols = []
    for h in horizontalSymbols:
        suffixH = h[numberOfReservedSymbols :]
        temp = []
        for vSymbol in verticalSymbols:
            for prefix in prefixes:
                hSymbol = prefix + suffixH
                if inList[hSymbol][vSymbol] == 1:
                    temp.append(prefix)
        candidates[suffixH] = temp
        
    return candidates

def randomEuclidEncoding(symbolSize, constraintList):
    # First we generate a connectivity matrix based on the symbol size and the constraints
    connectivityMatrix, sequencesBinary, sequencesBases, violationsMatrix = connectivityMatrix(symbolSize, constraintList)
    reducedMatrix = copy.deepcopy(connectivityMatrix)
    isFullyConnected = np.all(reducedMatrix)
    numberOfReservedBits = 0
    while not isFullyConnected and (np.mod(reducedMatrix.shap[1], 2) == 0):
        reducedMatrix, isFullyConnected = reduceMatrix(reducedMatrix)
        numberOfReservedBits = numberOfReservedBits + 1
    
    