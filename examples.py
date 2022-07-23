# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:26:38 2022

@author: omers
"""

from produceEuclidEncoding import *
from analysisFunctions import windowedGCcontent
import matplotlib.pyplot as plt
import mapping
import produceEuclidEncoding
import convolutional

EXAMPLE_SEQUENCE = 'TTTTTCGATTTTTGACCACGAAACCACCCTGACGGCCGCGAAAAAAAACATCCATCAAAAAGGGAAAAAAAAAAAAAAGAAAAAAAGGAAACAATTAAAAAAAGAAAAAAAAAAAAAACTAAAAAAAAAAAAGACGAAACACAAAAAAAAAGAAAAAAAAAAAAAAGTAAAAAAAAAAAAGATGAAACACGGAAAAAACCAAAAAAAAAAAAAAACAAAAAAAAAAAAGCGGAAACACGTAAAAAACCAAAAAAAAAAAAAAACAA'

EXAMPLE_HOMOPOLYMER = 'A' * 266

EXAMPLE_STR1 = 'AG' * 133

def dissertationExample1(windowSize = 20, symbolSize = 8):
    example1ConstraintList = {'gcMin': 0.25, 'gcMax': 0.75, 'runLength': 10}#,
    c, cd, vsym, hsym = euclidCandidates(constraintList = example1ConstraintList, symbolSize = symbolSize) 
    #numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols = makeFSM(c, vsym, hsym, minimiseReservedValue)
    numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols = makeFSM(c, vsym, hsym, trackGClevel)
    triggerLength = len(horizontalSymbols[0])
    status, binaryTextLine, binaryLineNumerical = mapping.dnaToBinary(EXAMPLE_SEQUENCE)
    nucStream = mapping.binaryStreamToBases(binaryTextLine)
    #print(nucStream == mapping.binaryStreamToBases(binaryTextLine))
    initialState = '0'*2*symbolSize
    euclidFSM = convolutional.makeEuclidFSM(verticalSymbols, horizontalSymbols, outputFSM, initialState)
    if (len(binaryTextLine) % triggerLength) != 0:
        padding = '0' * (triggerLength - (len(binaryTextLine) % triggerLength))
        binaryTextLine = binaryTextLine + padding
    encodedStream = convolutional.FSMdictionaryEncoder(binaryTextLine, euclidFSM)
    flatStream = ''
    for sublist in encodedStream:
        flatStream = flatStream + sublist
    encodedNucStream = mapping.binaryStreamToBases(flatStream)
    
    fig, ax = plt.subplots()
    status, slidingPointsSource, gcContentSource, GC, AT, other = windowedGCcontent(EXAMPLE_SEQUENCE, 20)
    txt = 'Source sequence, window size ' +str(windowSize)
    ax.plot(slidingPointsSource, gcContentSource, linewidth=2.0, label = txt)
    
    status, slidingPointsEncoded1, gcContentEncoded1, GC, AT, other = windowedGCcontent(encodedNucStream, 20)
    txt = 'Encoded, window size ' +str(windowSize)
    ax.plot(slidingPointsEncoded1, gcContentEncoded1, linewidth=2.0, label = txt)
    
    
    ax.set(ylim=(0, 1))
    ax.grid(True)
    plt.legend()
    plt.show()
    return numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols, encodedNucStream






