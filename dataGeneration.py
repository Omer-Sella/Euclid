# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 12:54:08 2022

@author: omers
"""


import numpy as np
import mapping
import produceEuclidEncoding
import mapping
import convolutional


"""
#basesVector = [0,1,2,3]
np.random.seed(7134066)
basesVector = ['A','C','T','G']
probabilityVector = np.array([0.7, 0.3, 0, 0])
# Very interesting, it seems if I pass probabilityVector without writing p =, the result does not have the desired distribution.
dataParsed = np.random.choice(basesVector, 2560, p = probabilityVector)
data = ""

for d in dataParsed:
    if d == 0:
        data = data + 'A'
    elif data == 1:
        data = data + 'C'
    elif data == 2:
        data = data + 'T'
    else:
        data = data + 'G'

for d in dataParsed:
    data = data + d


subSequences, occurences = mapping.produceStatsFromDnaStream(data)

mapping.subSequenceGraphics(subSequences, occurences)
"""

def binaryFileToDnaUsingEuclid(inputFilename, lineLength, outputFilename = None, format = 'IDT'):
    
    numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols = produceEuclidEncoding.testEuclid()
    euclidFSM = convolutional.makeEuclidFSM(verticalSymbols, horizontalSymbols, outputFSM)
    triggerLength = len(horizontalSymbols[0])
    #print("Trigger length is: " + str(triggerLength))
    #print("euclidFSM thinks that trigger len is: " + str(euclidFSM.stepSize) )
    numOfBits = 8
    
    
    if lineLength % lineLength != 0:
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
            #b"" is the empty byte
            inputLine = b""
            while byte != b"":
                inputLine = inputLine + byte
                lineCounter = lineCounter + 4
                if lineCounter == actualLineLenghth:
                    
                    bitStream = (bin(int(inputLine.hex(), 16))[2:]).zfill(numOfBits)
                    #print(bitStream)
                    #print(len(bitStream))
                    #print(lineCounter)
                    #Pad with zeros so that FSMEncoder will be happy, i.e.: streamLength % triggerLength == 0
                    padding = (triggerLength - (len(bitStream) % triggerLength)) * '0'
                    #print(padding)
                    bitStream = bitStream + padding
                    #print(bitStream)
                    
                    codedLine = convolutional.FSMdictionaryEncoder(bitStream, euclidFSM)
                    #print(codedLine)
                    # Omer Sella: flatStream gives you the un-chopped encoded stream
                    flatStream = ""
                    for sublist in codedLine:
                        for item in sublist:
                            flatStream = flatStream + str(item)
                    #print(flatStream)
                    basesLine = mapping.binaryStreamToBases(flatStream)
                    line = line + basesLine
                    outFile.write(line)
                    sequenceNumber = sequenceNumber + 1
                    if format == 'IDT':
                        line = "seq" + str(sequenceNumber) + "\t"
                        
                    else:
                        line = ""
                    inputLine = b""
                    lineCounter = 0
                # Read new byte
                byte = inFile.read(1)

