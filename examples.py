# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:26:38 2022

@author: omers
"""
import os
import sys
projectDir = os.environ.get('EUCLID')
if projectDir == None:
     projectDir = "D:/Euclid/"
sys.path.insert(1, projectDir)
from produceEuclidEncoding import *
from analysisFunctions import windowedGCcontent, gcVariance
import matplotlib.pyplot as plt
import mapping
import produceEuclidEncoding
import convolutional
import scipy.io
import numpy as np
import binaryToDna

EXAMPLE_SEQUENCE = 'TTTTTCGATTTTTGACCACGAAACCACCCTGACGGCCGCGAAAAAAAACATCCATCAAAAAGGGAAAAAAAAAAAAAAGAAAAAAAGGAAACAATTAAAAAAAGAAAAAAAAAAAAAACTAAAAAAAAAAAAGACGAAACACAAAAAAAAAGAAAAAAAAAAAAAAGTAAAAAAAAAAAAGATGAAACACGGAAAAAACCAAAAAAAAAAAAAAACAAAAAAAAAAAAGCGGAAACACGTAAAAAACCAAAAAAAAAAAAAAACAA'

EXAMPLE_HOMOPOLYMER = 'A' * 266

EXAMPLE_HOMOPOLYMERS = 'A' * 66 + 'C' * 66 + 'T' * 66 + 'G' * 66 + 'ACTG'

EXAMPLE_STR1 = 'ACTG' * 133

def saveFSM(states, triggers, transitionDictionary, path = projectDir, fileName = 'example'):
    workspaceDict = {}
    workspaceDict['states'] = states
    workspaceDict['triggers'] = triggers
    workspaceDict['transitionDictionary'] = transitionDictionary
    fileNameWithPath = path + "/" + fileName + '.mat'
    scipy.io.savemat(fileNameWithPath, workspaceDict)
    return fileNameWithPath
    

def dissertationExample1(windowSize = 20, symbolSize = 5, mechanism = trackGClevel, uncodedSequence = None, fileName = "example", gcMin = 0.25, gcMax = 0.75, runLength = 10):
    workspaceDict = {}
    
    example1ConstraintList = {'gcMin': gcMin, 'gcMax': gcMax, 'runLength': runLength}#,
    workspaceDict['constraintList'] = example1ConstraintList
    c, cd, vsym, hsym = euclidCandidates(constraintList = example1ConstraintList, symbolSize = symbolSize) 
    
    #numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols = makeFSM(c, vsym, hsym, minimiseReservedValue)
    numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols = makeFSM(c, vsym, hsym, mechanism)
    
    triggerLength = len(horizontalSymbols[0])
    if uncodedSequence == None:
        uncodedSequence = EXAMPLE_SEQUENCE

    status, binaryTextLine, binaryLineNumerical = mapping.dnaToBinary(uncodedSequence)
    
    
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
    status, slidingPointsSource, gcContentSource, GC, AT, other = windowedGCcontent(EXAMPLE_SEQUENCE, windowSize)
    txt = 'Source sequence, window size ' +str(windowSize)
    ax.plot(slidingPointsSource, gcContentSource, linewidth=2.0, label = txt)
    status, slidingPointsEncoded1, gcContentEncoded1, GC, AT, other = windowedGCcontent(encodedNucStream, windowSize)
    txt = 'Encoded, window size ' +str(windowSize)
    ax.plot(slidingPointsEncoded1, gcContentEncoded1, linewidth=2.0, label = txt)
    ax.set(ylim=(0, 1))
    ax.grid(True)
    ax.set_xlabel("Index of window start position. Window size = " + str(windowSize) + "bases", size = 24)
    ax.set_ylabel("GC content normalised", size = 24)
    ax.set_title('GC content calculated on a sliding windows', size = 24)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 24)
    ax.vlines(slidingPointsSource[-1], 0, 1, colors = 'red', linewidth = 2.5, linestyle = '--')
    ax.vlines(slidingPointsEncoded1[-1], 0, 1, colors = 'red', linewidth = 2.5, linestyle = '--')
    averageGCSourceTxt = "Average sliding window GC content " + str(np.average(gcContentSource))
    ax.hlines(np.average(gcContentSource), slidingPointsSource[0], slidingPointsSource[-1], colors = 'red', linestyle = '-', linewidth = 2.0, label = averageGCSourceTxt)
    #ax.hlines(np.average(gcContentSource), slidingPointsSource[0], slidingPointsSource[-1], colors = 'red', linestyle = '-', linewidth = 2.0)
    averageGCEncodedTxt = "Average sliding window GC content " + str(np.average(gcContentEncoded1))
    ax.hlines(np.average(gcContentEncoded1), slidingPointsEncoded1[0], slidingPointsEncoded1[-1], colors = 'green', linestyle = '-', linewidth = 2.0, label = averageGCEncodedTxt)
    #ax.hlines(max(gcContentEncoded1), slidingPointsEncoded1[0], slidingPointsEncoded1[-1], colors = 'red', linestyle = '-', linewidth = 2.0)
    plt.legend(fontsize = 24)
    
    workspaceDict['mechanism'] = str(mechanism)
    workspaceDict['gcMin'] = gcMin
    workspaceDict['gcMax'] = gcMax
    workspaceDict['runLength'] = runLength
    workspaceDict['constraintList'] = example1ConstraintList
    workspaceDict['c'] = c
    workspaceDict['cd'] = cd
    workspaceDict['vsym'] = vsym
    workspaceDict['hsym'] = hsym
    workspaceDict['uncodedSequence'] = uncodedSequence
    workspaceDict['binaryTextLine'] = binaryTextLine
    workspaceDict['binaryLineNumerical'] = binaryLineNumerical
    workspaceDict['nucStream'] = nucStream
    workspaceDict['encodedStream'] = encodedStream
    workspaceDict['encodedNucStream'] = encodedNucStream
    workspaceDict['outputDictionary'] = outputDictionary
    workspaceDict['outputFSM'] = outputFSM
    fileNameWithPath = projectDir + "/" + fileName + '.mat'
    scipy.io.savemat(fileNameWithPath, workspaceDict)
    #plt.show()
    plt.tight_layout()
    figureFileNameWithPath = projectDir + "/" + fileName + '.png'
    plt.savefig(fname = figureFileNameWithPath)
    return numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols, encodedNucStream

def openVirusFile(filePath):
    # Function assumes a path to a txt file that contains DNA, possibly with newline breaks \n, possibly with spaces " " and possibly lowercase letters.
    # No safety provided (so if you have line numbers for example, remove them)
    dna = ''
    with open(filePath, "r") as virusFile:
        for newLine in virusFile:
            lineList = virusFile.readline()
            lineList = lineList.split(" ")
            print(lineList)
            #line = line.strip(" ")
            #print(line)
            line = ''
            for word in lineList:
                temp = word.strip('\n')
                line = line + temp.upper()
            dna = dna + line
    return dna



def graphicsForExamples(slidingPoints, gcContent, windowSize, fileName, title = None):
    fig, ax = plt.subplots()
    txt = 'Source sequence, window size ' +str(windowSize)
    ax.plot(slidingPoints, gcContent, linewidth=2.0, label = txt)
    ax.set(ylim=(0, 1))
    ax.grid(True)
    ax.set_xlabel("Index of window start position. Window size = " + str(windowSize) + "bases", size = 24)
    ax.set_ylabel("GC content normalised", size = 24)
    if title == None:
        title = 'GC content calculated on a sliding windows'
    ax.set_title(title, size = 24)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 24)
    #ax.vlines(slidingPointsSource[-1], 0, 1, colors = 'red', linewidth = 2.5, linestyle = '--')
    #ax.vlines(slidingPointsEncoded1[-1], 0, 1, colors = 'red', linewidth = 2.5, linestyle = '--')
    averageGCSourceTxt = "Average sliding window GC content " + str(np.average(gcContent))
    ax.hlines(np.average(gcContent), slidingPoints[0], slidingPoints[-1], colors = 'red', linestyle = '-', linewidth = 2.0, label = averageGCSourceTxt)
    #ax.hlines(np.average(gcContentSource), slidingPointsSource[0], slidingPointsSource[-1], colors = 'red', linestyle = '-', linewidth = 2.0)
    #averageGCEncodedTxt = "Average sliding window GC content " + str(np.average(gcContentEncoded1))
    #ax.hlines(np.average(gcContentEncoded1), slidingPointsEncoded1[0], slidingPointsEncoded1[-1], colors = 'red', linestyle = '-', linewidth = 2.0, label = averageGCEncodedTxt)
    #ax.hlines(max(gcContentEncoded1), slidingPointsEncoded1[0], slidingPointsEncoded1[-1], colors = 'red', linestyle = '-', linewidth = 2.0)
    plt.legend(fontsize = 24)
    plt.tight_layout()
    figureFileNameWithPath = projectDir + "/" + fileName + '.png'
    plt.savefig(fname = figureFileNameWithPath)    

def produceGraphicsForMethodology(sequence = None):
    if sequence == None:
        sequence = EXAMPLE_SEQUENCE
    subSequences, occurences = binaryToDna.produceStatsFromSequence(sequence)
    binaryToDna.subSequenceGraphics(subSequences[0:4], occurences[0:4])
    binaryToDna.subSequenceGraphics(subSequences[0:20], occurences[0:20])
    binaryToDna.subSequenceGraphics(subSequences, occurences)
    return
    
    
    
def openExample(filePath):
    workSpaceDictionary = scipy.io.loadmat(filePath)
    uncodedStream = (workSpaceDictionary['uncodedSequence'])[0]
    encodedStream = (workSpaceDictionary['encodedNucStream'])[0]
    gcVariance(encodedStream, 50)
    #varianceArray, slidingPoints, gcContent, argMaxGCCOntent, argMinGCContent = gcVariance(encodedStream, 50)
    #return workSpaceDictionary, uncodedStream, encodedStream, varianceArray, slidingPoints, gcContent, argMaxGCCOntent, argMinGCContent
    return