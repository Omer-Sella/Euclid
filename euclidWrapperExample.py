# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 09:08:52 2022

@author: omers
"""
import os
import sys
projectDir = os.environ.get('EUCLID')
assert(projectDir != None)
sys.path.insert(1, projectDir)
import numpy as np
from produceEuclidEncoding import makeFSM, trackGClevel, completelyRandom
import mapping
import convolutional
def encodeUsingEuclid(euclidFilePath, assignmentMechanism, data):
    """

    Parameters
    ----------
    euclidFilePath : TYPE
        DESCRIPTION.
    mechanism : is a sting that matches one of the assignment mechanisms in produceEuclidEncoding.py. There is currently no safety.
        DESCRIPTION.
    data : Assumed to be a string of nucleotides (A,C,T,G) capital letters. There is an option to accept binary data,
    and in fact, later on we convert the nucleotide stream to binary, since the euclidFSM we use assumes the data is a binary string.
    Returns 
    codedData : a string of nucleotides (and could return binary if Ryan decides to do so)
    -------

    """
  
      
    if assignmentMechanism == 'trackGC':
        mechanism = trackGClevel
    elif assignmentMechanism == 'random':
        mechanism = completelyRandom
    else:
        #Omer Sella: need to make sure that the mechanism is a valid mechanism from produceEuclidEncoding.
        #This is a place holder.
        mechanism = completelyRandom
        
    workSpaceDictionary = np.load(euclidFilePath ,allow_pickle='TRUE').item()
    c = workSpaceDictionary['c']
    vsym = workSpaceDictionary['vsym']
    hsym = workSpaceDictionary['hsym']
    numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols = makeFSM(c, vsym, hsym, mechanism)
    triggerLength = len(horizontalSymbols[0])
    symbolSize = len(vsym[0])
    print(symbolSize)
    initialState = '0'*symbolSize
    print(initialState)
    #####################################################################################################################################################
    euclidFSM = convolutional.makeEuclidFSM(verticalSymbols, horizontalSymbols, outputFSM, initialState)
    ############### What we need is an object called euclidFSM - if it's given to us, there is no need to do everything above this line #################
    

    ##################### This is where data is being encoded ######################################
    # Assuming data is a long string of bases (A,C,T,G all capital letters) we convert it to a binary string
    status, binaryTextLine, binaryLineNumerical = mapping.dnaToBinary(data)
    # Prepreparation: in case the binary stream is not an integer number of triggers, we pad it with 0s
    if (len(binaryTextLine) % triggerLength) != 0:
        padding = '0' * (triggerLength - (len(binaryTextLine) % triggerLength))
        binaryTextLine = binaryTextLine + padding
    
    # Now comes the encoding
    encodedStream = convolutional.FSMdictionaryEncoder(binaryTextLine, euclidFSM)
    # encodedStream is actually a LIST of output symbols, meaning a list of small binary strings.
    # So the next lines are to make it into a single line.
    flatStream = ''
    for sublist in encodedStream:
        flatStream = flatStream + sublist
    # Now we map the binary flat stream back into bases (A,C,T,G)
    codedData = mapping.binaryStreamToBases(flatStream)
    return codedData