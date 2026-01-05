# -*- coding: utf-8 -*-
"""

@author: omers
"""



def encodeData(stringData, pathToNpyEncodingFile = "c:/users/omer/Euclid/Euclid4.npy", pathToLeaveCodec = "c:/users/omer/Euclid/temp/"):
    """
    encode data is a wrapper for the function convolutional.FSMdictionaryEncoder
    encode the string data given in argyument stringData using a finite state machine encoding given in the function convolutional.FSMdictionaryEncoder
    according to the parameters given in pathToNpyEncodingFile
    stringData - assumed to be a string (no safety)
    pathToNpyEncodingFile - an npy file that contains variables to generate a Euclid FSM encoder
    output - string of leters from {A,C,T,G}
    """
    import numpy as np
    import convolutional
    import mapping
    from mapping import makeFakeData
    
    
    #Load work space dictionary that was generated using our previous Euclid work
    workSpaceDictionary = np.load(pathToNpyEncodingFile ,allow_pickle='TRUE').item()
    
    #Now create a finite state machine that will be used to encode the data
    euclidFSM = convolutional.makeEuclidFSM(workSpaceDictionary['vsym'], workSpaceDictionary['hsym'], workSpaceDictionary['outputFSM'], '00000000')
    
    #Process the stringData into a binary string, with possible additional zero padding so that it could be encoded by convolutional.FSMdictionaryEncoder
    asciiEncoded = stringData.encode('ascii')
    binaryEncoded = ['{0:08b}'.format(i) for i in asciiEncoded]
    binaryData = ''.join(b for b in binaryEncoded)
    triggerLength = euclidFSM.stepSize
    padding = '0' * (triggerLength - (len(binaryData) % triggerLength))
    paddedBinaryData = binaryData + padding
    #print(paddedBinaryData)
    #This is where encoding happens - the zero-padded binary-converted stringData is passed along with the FSM instance
    encodedStream = convolutional.FSMdictionaryEncoder(paddedBinaryData, euclidFSM)
    #print(encodedStream)
    #The encoded stream is a list of symbols, we're joining it to a continuous binary string 
    flatStream = "".join(encodedStream)
    #print(flatStream)
    #Convert the binary string to bases
    dnaStream = mapping.binaryStreamToBases(flatStream) 
    messageLength = len(dnaStream)
    #Add alignment markers
    #dnaStream = alignmentMarkers.BARKER11_FLAT + dnaStream + alignmentMarkers.BARKER11_FLAT
    #dnaStream = alignmentMarkers.IMPERIAL_MARKER1_FLAT + dnaStream + alignmentMarkers.BARKER11_FLAT
    
    #Padd UP with alignment markers to get 140 bases 
    #if len(dnaStream) < numberOfNucleotides:
    #    extraPadding = numberOfNucleotides * alignmentMarkers.BARKER11_FLAT
    #    #dnaStream = dnaStream + extraPadding[0 : 140 - len(dnaStream) - 1]
    #    dnaStream = dnaStream + extraPadding[0 : numberOfNucleotides - len(dnaStream)]
    
    codecDict = {}    
    codecDict['paddingLength'] = len(padding)
    codecDict['npyFilePath'] = str(pathToNpyEncodingFile)
    #codecDict['messageStart'] = len(alignmentMarkers.BARKER11_FLAT)
    #codecDict['messageEnd'] = len(alignmentMarkers.BARKER11_FLAT) + messageLength
    #codecDict['preambleAlignmentMarker'] = alignmentMarkers.IMPERIAL_MARKER1_FLAT
    #codecDict['paddingAlignmentMarker'] = alignmentMarkers.BARKER11_FLAT
    
    fileNameWithPath = pathToLeaveCodec + "/" + "codec.npy"
    np.save(fileNameWithPath, codecDict)
        
    return dnaStream

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--stringData', default = "Imperial", type = str)
    parser.add_argument('--pathToNpyEncodingFile', type=str, default = "c:/users/omer/Euclid/Euclid.npy")
    parser.add_argument('--pathToLeaveCodec', type=str, default = "c:/users/omer/Euclid/tests/")
    parser.add_argument('--generateSyntheticSequencingFiles', type = str, default = "True")
    parser.add_argument('--numberOfNucleotides', type = int, default = 140)
    args = parser.parse_args()
    if args.generateSyntheticSequencingFiles == "False":
        args.generateSyntheticSequencingFiles = False        
    
    dnaStream = encodeData(args.stringData, args.pathToNpyEncodingFile, args.pathToLeaveCodec, args.generateSyntheticSequencingFiles, args.numberOfNucleotides)
    print(dnaStream)
    
