# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 10:09:53 2021

@author: Omer Sella
"""
import numpy as np
import copy
import matplotlib.pyplot as plt
import produceEuclidEncoding as euclid
BIG_NUMBER = 1000000 #np.Inf
LOCAL_PRNG = np.random.RandomState(42)

class FSM():
    
    """
    A finite state machine (FSM) a graph made of states (vertices) and a connectivity matrix (transitionTable).
    If the FSM is at state s0, then on input (trigger) i, it will emmit an output e, and transition from state 
    s0 (the present state) to a next state, determined by the transitionTable.
    An initialState has to be given so that we know how to start.
    
    """
    def __init__(self, states, triggers, outputTable, transitionTable, initialState):
        # Safety: all triggers must be of the same fixed length
        fixedLength = len(triggers[0])
        assert all(len(t) == fixedLength for t in triggers)
    
        self.triggers = triggers
        self.numberOfStates = len(states)
        self.stepSize = fixedLength
        self.transitionTable = transitionTable
        self.outputTable = outputTable
        # Enumerate the possible triggers according to the order they were given.
        self.triggerDictionary = self.generateDictionary(triggers)
        
        # Enumerate the possible states according to the order they were given.
        self.stateDictionary = self.generateDictionary(states)
        self.presentState = initialState
        self.presentStateCoordinate = initialState
        
    def checkStreamLength(self, length):
        result = False
        if length % self.stepSize == 0:
            result = 'OK'
        return result
    
    def generateDictionary(self, keys):
        newDictionary = dict()
        #entries = np.arange(0, self.outputTable.shape[0])
        i = 0
        for t in keys:
            #assert t not in newDictionary,  'Multiple appearences for trigger'
            newDictionary[str(t)] = i
            i = i + 1
        return newDictionary
    
    def step(self, trigger):
        #Init output
        output = ''
        #Get the number of the trigger
        triggerCoordinate = self.triggerDictionary[str(trigger)]
        #Get the output from the outputTable, matching to the number of state we are in and the number of the trigger.
        output = self.outputTable[self.presentStateCoordinate][triggerCoordinate]
        nextState = self.transitionTable[self.presentStateCoordinate, triggerCoordinate]
        nextStateCoordinate = self.stateDictionary[str(nextState)]
        self.presentState = nextState
        self.presentStateCoordinate = nextStateCoordinate
        return output
    
    def stepDictionary(self, trigger):
        #Init output
        output = ''
        #Get the output from the outputTable, matching to the number of state we are in and the number of the trigger.
        #print(self.presentState)
        #print(type(self.presentState))
        output = self.outputTable[self.presentState][trigger]
        nextState = self.transitionTable[self.presentState][trigger]
        nextStateCoordinate = self.stateDictionary[nextState]
        self.presentState = nextState
        self.presentStateCoordinate = nextStateCoordinate
        return output
    
    
    
    def getNextPossibleStates(self, state):
        nextStates = self.transitionTable[state, :]
        return nextStates
    
    def getNextPossibleOutputs(self, state):
        nextOutputs = self.outputTable[state]
        return nextOutputs
    
def trellisGraphics(numberOfStates, time, stateTraversalList, scores = None):
    states = np.arange(numberOfStates)                
    timeAxis = np.arange(0,time)
    statesTiled = np.tile(states, time)
    timeRepeat = np.repeat(timeAxis, numberOfStates)
    fig, ax = plt.subplots()
    ax.set_yticks(states)
    ax.set_ylabel('States')
    ax.set_xlabel('Received symbols')
    ax.set_title('Path traversal')
    scatter = ax.scatter(timeRepeat, statesTiled)
    for t in stateTraversalList:
        pathTime = np.arange(0, len(t))
        plot = ax.plot(np.arange(0,len(t)), t)
    if scores is not None:
        assert (len(scores) == numberOfStates)
        for idx in range(numberOfStates):
            ax.text(time, idx, scores[idx])
    plt.show()
    return
    
def FSMEncoder(streamIn, FSM, graphics = False):
    """
    streamIn is a stream of (symbols). There is no safety over their content validity,
    but we do verify that the stream of symbols could be chopped into an integer number of triggers.
    For example: if the stream is made of bits, and we need 5 bits per trigger, then the stream has to have length = 5*k for some k.
    """
    assert FSM.checkStreamLength(len(streamIn)) == 'OK' , "Input to convolutional encoder must be an integer multiple of FSM.stepSize"
    numberOfSteps = len(streamIn) // FSM.stepSize
    i = 0
    encodedStream = []
    while i < numberOfSteps:
        trigger = streamIn[i * FSM.stepSize: (i + 1) * FSM.stepSize]
        #print(trigger)
        trigger = list(trigger)
        #print(trigger)
        output = FSM.step(trigger)
        #print(output)
        encodedStream.append(output)
        i = i + 1
        #print(encodedStream)
    return encodedStream

def FSMdictionaryEncoder(streamIn, FSM, graphics = False):
    """
    streamIn is a stream of (symbols). There is no safety over their content validity,
    but we do verify that the stream of symbols could be chopped into an integer number of triggers.
    For example: if the stream is made of bits, and we need 5 bits per trigger, then the stream has to have length = 5*k for some k.
    """
    assert FSM.checkStreamLength(len(streamIn)) == 'OK' , "Input to convolutional encoder must be an integer multiple of FSM.stepSize"
    numberOfSteps = len(streamIn) // FSM.stepSize
    i = 0
    encodedStream = []
    while i < numberOfSteps:
        
        trigger = streamIn[i * FSM.stepSize: (i + 1) * FSM.stepSize]
        output = FSM.stepDictionary(trigger)
        encodedStream.append(output)
        i = i + 1
    return encodedStream


class path():
    def __init__(self, initialState):
        self.traversedStates = [initialState]
        self.scores = [0]
        self.presentScore = 0
        self.pathTriggers = []
        self.pathEmitted = []
        
    def presentState(self):
        return self.traversedStates[-1]
    
    
    def appendToPath(self, extension):
        nextState = extension[0]
        trigger = extension[1]
        nextOutput = extension[2]
        addedScore = extension[3]
        #print("*** before step:")
        #print("*** " + str(self.pathTriggers))
        #print("*** " + str(self.traversedStates))
        #print("*** ")
        #print("*** ")
        self.pathTriggers.append(trigger)
        self.pathEmitted.append(nextOutput)
        self.scores.append(addedScore)
        self.traversedStates.append(nextState)   
        self.presentScore = self.presentScore + addedScore
        return
        


def viterbiDecoder(numberOfStates, initialState, fanOutFunction, observedSequence, symbolsPerStateTransition):
    # Viterbi decoder inspired by the implementation suggested by Todd K. Moon, programming laboratory 10, page 528.
    # More explanations on the Viterbi decoder are found on page 473 of the same book.
    # A metric function (!) that accepts a set of states p, next state q and observed stream r,
    # and returns the branch metric present state, next state and returns 
    
    # There is no safety of the content of the observedSequence, but the observedSequence has to be chopped into observable state transitions.
    # This means that *this version of the Viterbi decoder does not support insertions or deletions.
    assert len(observedSequence) % symbolsPerStateTransition == 0
    newPath = path(initialState)
    paths = [newPath]
    i = 0
    while i < (len(observedSequence) // symbolsPerStateTransition):
        observedOutput = observedSequence[i * symbolsPerStateTransition : (i + 1) * symbolsPerStateTransition]
        newPaths = []
        for p in paths:
            extensions = fanOutFunction(p.presentState(), observedOutput, i)
            for extension in extensions:
                newPath = path(0)
                newPath = copy.deepcopy(p)
                newPath.appendToPath(extension)
                newPaths.append(newPath)
        paths = newPaths
        #Omer Sella: Here there is usually pruning, i.e.: getting rid of costly candidates, but not in this version.
        i = i + 1

       
    #Omer Sella: Now let's find "the" most likely path, which is the that has the LOWEST score (so score is like loss)
    lowestScore  = BIG_NUMBER
    numberOfEquallyMostLikely = 1
    for p in paths:
        if p.presentScore < lowestScore:
            # New lowest score path found !
            # Set the lowest score to the new lowest score
            lowestScore = p.presentScore
            # Set ost likely path pointer to this path
            mostLikelyPath = p
            # Omer Sella: bug fix: if a "lowest score path" was found, then reset the number of of equally likely paths to 1
            numberOfEquallyMostLikely = 1
        else:
            if p.presentScore == lowestScore:
                numberOfEquallyMostLikely = numberOfEquallyMostLikely + 1
                

    # Omer Sella: Viterbi is supposed to return the original input, it could also return paths
    # So we first return the most likely path, if there is more than one then numberOfEquallyMostLikely will be > 1
    # Then we return all paths 
    # BUG: numberOfEquallyLikelyPaths should be 1 if no errors.
    return mostLikelyPath, numberOfEquallyMostLikely, paths
    
def genericFanOutFunction(myFSM, presentState, observedOutput, timeStep, additionalInformation):
    # some useful distances in this package, we may want to try and use it.
    from scipy.spatial import distance
    nextPossibleStates = myFSM.getNextPossibleStates(presentState)
    nextPossibleOutputs = myFSM.getNextPossibleOutputs(presentState)
    triggers = myFSM.triggers
    extensions = []    
    nextPossibleScores = []
    for output in nextPossibleOutputs:
        # compute the score of output with respect to the observedOutput
        # print("*** output is:")
        # print(output)
        # print("*** observedOutput is:")
        # print(observedOutput)
        # Omer Sella: scipy.distance.hamming(x,y) return number_of_coordinates_where_different / number_of_coordinates. Sequences must have identical length.
        score = distance.hamming(output, observedOutput)
        #print(score)
        nextPossibleScores.append(score)
    extensions = []
    # Omer Sella: safety
    assert (len(nextPossibleStates) == len(nextPossibleOutputs))
    for idx in range(len(nextPossibleStates)):
        extensions.append( [nextPossibleStates[idx], triggers[idx], nextPossibleOutputs[idx], nextPossibleScores[idx]])
    #print("*** extensions are: ")
    #print(extensions)
    return extensions


def viterbiDecoderWithFlagging(numberOfStates, initialState, fanOutFunction, observedSequence, symbolsPerStateTransition, produceGraphics = False):
    # Viterbi decoder inspired by the implementation suggested by Todd K. Moon, programming laboratory 10, page 528.
    # More explanations on the Viterbi decoder are found on page 473 of the same book.
    # A metric function (!) that accepts a set of states p, next state q and observed stream r,
    # and returns the branch metric present state, next state and returns 
    
    # There is no safety of the content of the observedSequence, but the observedSequence has to be chopped into observable state transitions.
    # This means that *this version of the Viterbi decoder does not support insertions or deletions.
    assert len(observedSequence) % symbolsPerStateTransition == 0
    
    indelFlag = False
    indelEstimatedLocation = BIG_NUMBER
    
    newPath = path(initialState)
    paths = [newPath]
    i = 0
    scoreVector = np.ones(numberOfStates) * BIG_NUMBER
    numberOfEquallyLikelyPathsVector = np.zeros(numberOfStates)
    while i < (len(observedSequence) // symbolsPerStateTransition):
        observedOutput = observedSequence[i * symbolsPerStateTransition : (i + 1) * symbolsPerStateTransition]
        newPaths = []
        for p in paths:
            extensions = fanOutFunction(p.presentState(), observedOutput, i)
            for extension in extensions:
                newPath = path(0)
                newPath = copy.deepcopy(p)
                newPath.appendToPath(extension)
                newPaths.append(newPath)
        paths = newPaths
       
        # Now we discard unlikely paths
        for s in range(numberOfStates):
            scores = []
            scoreVector[s] = BIG_NUMBER
            for p in paths:
                if p.presentState() == s:
                    scores.append(p.presentScore)
            if len(scores) > 0:
                lowestScoreForThisState = min(scores)
                scoreVector[s] = lowestScoreForThisState
                numberOfEquallyLikelyPathsVector[s] = scores.count(lowestScoreForThisState)
                for p in paths:
                    if p.presentState() == s:
                        if p.presentScore > lowestScoreForThisState:
                            paths.remove(p)
        #print("***At step " + str(i))
        #print("***most likely vector is: " + str(numberOfEquallyLikelyPathsVector))
        #print("***score vector is: " + str(scoreVector))
        
        # Omer Sella: the following line is just a guess based on observations
        if np.sum(numberOfEquallyLikelyPathsVector) > 2 * numberOfStates and indelFlag == False:
            indelFlag = True
            indelEstimatedLocation = i
        if produceGraphics:
            stateTraversalList = []
            for p in paths:
                stateTraversalList.append(p.traversedStates)
            trellisGraphics(numberOfStates, (len(observedSequence) // symbolsPerStateTransition), stateTraversalList, numberOfEquallyLikelyPathsVector)#scoreVector)
        i = i + 1
           
    lowestScore = min(scoreVector)
    mostLikelyPaths = []
    for p in paths:
        if p.presentScore == lowestScore:
            mostLikelyPaths.append(p)
    # Omer Sella: Viterbi is supposed to return the original input, it could also return paths
    # So we first return the most likely path, if there is more than one then numberOfEquallyMostLikely will be > 1
    # Then we return all paths 
    # BUG: numberOfEquallyLikelyPaths should be 1 if no errors.
    return mostLikelyPaths, scoreVector, numberOfEquallyLikelyPathsVector, paths, indelFlag, indelEstimatedLocation


def exampleOneHalfConvolutional():
    states = [0,1,2,3]
    triggers = [[0], [1]]
    nextStateTable = np.array([[0,1], [2,3], [0,1], [2,3]])
    outputTable = [[[0,0], [1,1]],
                   [[0,1], [1,0]],
                   [[1,1], [0,0]],
                   [[1,0], [0,1]]]
    symbolSize = 2
    initialState = 0
    myFSM = FSM(states, triggers, outputTable, nextStateTable, initialState)
    stream = LOCAL_PRNG.randint(0,2,10)
    encodedStream = FSMEncoder(stream, myFSM)
    print(encodedStream)
    flatStream = []
    for sublist in encodedStream:
        for item in sublist:
            flatStream.append(item)
    print(flatStream)
    def myFanOutFunction(state, observation, time):
        return genericFanOutFunction(myFSM, state, observation, time, None)

    mostLikelyPaths, scoreVector, numberOfEquallyLikelyPathsVector, pathsWithFlagging, indelFlag, indelEstimatedLocation = viterbiDecoderWithFlagging(len(states), initialState, myFanOutFunction, flatStream, symbolSize)

    return stream, encodedStream, mostLikelyPaths, scoreVector, numberOfEquallyLikelyPathsVector, pathsWithFlagging, indelFlag, indelEstimatedLocation


def exampleTwoThirdsConvolutional():
    states = [0,1,2,3,4,5,6,7]
    # Omer Sella: triggers are the raw data that we want to encode, post processing i.e.: chopped into blocks that the FSM likes (2 bits in our case).
    triggers = [[0,0], [0,1], [1,0], [1,1]]
    #Bug: the triggers came out as: [[0, 0], [1, 1], [1, 1], [0, 0], [0, 1]]
    nextStateTable = np.array([[0,1,2,3], [4,5,6,7], [1,0,3,2], [5,4,7,6], [2,3,0,1], [6,7,4,5] , [3,2,1,0], [7,6,5,4] ])
    outputTable = [[[0,0,0], [1,0,0], [0,1,0], [1,1,0]], 
                   [[0,0,1], [1,0,1], [0,1,1], [1,1,1]],
                   [[1,0,0], [0,0,0], [1,1,0], [0,1,0]],
                   [[1,0,1], [0,0,1], [1,1,1], [0,1,1]],
                   [[0,1,0], [1,1,0], [0,0,0], [1,0,0]],
                   [[0,1,1], [1,1,1], [0,0,1], [1,0,1]],
                   [[1,1,0], [0,1,0], [1,0,0], [0,0,0]],
                   [[1,1,1], [0,1,1], [1,0,1], [0,0,1]]]
    initialState = 0
    myFSM = FSM(states, triggers, outputTable, nextStateTable, initialState)
    stream = np.random.randint(0,2,10)
    encodedStream = FSMEncoder(stream, myFSM)
    

    # Omer Sella: flatStream gives you the un-chopped encoded stream
    flatStream = []
    for sublist in encodedStream:
        for item in sublist:
            flatStream.append(item)
    
    def myFanOutFunction(state, observation, time):
        return genericFanOutFunction(myFSM, state, observation, time, None)

    mostLikelyPath, numberOfEquallyLikelyPaths, paths = viterbiDecoder(8, 0, myFanOutFunction, flatStream, 3)
    mostLikelyPaths, scoreVector, numberOfEquallyLikelyPathsVector, pathsWithFlagging, indelFlag, indelEstimatedLocation = viterbiDecoderWithFlagging(8, 0, myFanOutFunction, flatStream, 3)

    return stream, encodedStream, paths, mostLikelyPath, numberOfEquallyLikelyPaths, mostLikelyPaths, scoreVector, numberOfEquallyLikelyPathsVector, pathsWithFlagging, indelFlag, indelEstimatedLocation

def testConvolutional_2_3():
    status = 'OK'
    stream, encodedStream, paths, mostLikelyPath, numberOfEquallyLikelyPaths, mostLikelyPaths, scoreVector, numberOfEquallyLikelyPathsVector, pathsWithFlagging, indelFlag, indelEstimatedLocation = exampleTwoThirdsConvolutional()
    if len(encodedStream) != len(stream)//2:
        status = 'FAIL'
    if mostLikelyPath.pathEmitted != encodedStream:
        status = 'FAIL'
    if len(mostLikelyPaths) > 1:
        status = 'FAIL'
    if mostLikelyPaths[0].presentScore != mostLikelyPath.presentScore:
        status = 'FAIL'

    return status

def testViterbiBitFlip():
    states = [0,1,2,3,4,5,6,7]
    # Omer Sella: triggers are the raw data that we want to encode, post processing i.e.: chopped into blocks that the FSM likes (2 bits in our case).
    triggers = [[0,0], [0,1], [1,0], [1,1]]
    #Bug: the triggers came out as: [[0, 0], [1, 1], [1, 1], [0, 0], [0, 1]]
    nextStateTable = np.array([[0,1,2,3], [4,5,6,7], [1,0,3,2], [5,4,7,6], [2,3,0,1], [6,7,4,5] , [3,2,1,0], [7,6,5,4] ])
    outputTable = [[[0,0,0], [1,0,0], [0,1,0], [1,1,0]], 
                   [[0,0,1], [1,0,1], [0,1,1], [1,1,1]],
                   [[1,0,0], [0,0,0], [1,1,0], [0,1,0]],
                   [[1,0,1], [0,0,1], [1,1,1], [0,1,1]],
                   [[0,1,0], [1,1,0], [0,0,0], [1,0,0]],
                   [[0,1,1], [1,1,1], [0,0,1], [1,0,1]],
                   [[1,1,0], [0,1,0], [1,0,0], [0,0,0]],
                   [[1,1,1], [0,1,1], [1,0,1], [0,0,1]]]
    initialState = 0
    myFSM = FSM(states, triggers, outputTable, nextStateTable, initialState)
    stream = LOCAL_PRNG.randint(0,2,20)
    encodedStream = FSMEncoder(stream, myFSM)

    # Omer Sella: flatStream gives you the un-chopped encoded stream
    flatStream = []
    for sublist in encodedStream:
        for item in sublist:
            flatStream.append(item)

    ############# Omer Sella: Now we flip a single bit
    corruptStream = copy.deepcopy(flatStream)
    corruptStream[0] = 1 - corruptStream[0]
    corruptStream[7] = 1 - corruptStream[0]
    corruptStream[10] = 1 - corruptStream[0]

    def myFanOutFunction(state, observation, time):
        return genericFanOutFunction(myFSM, state, observation, time, None)

    #mostLikelyPath, numberOfEquallyLikelyPaths, paths = viterbiDecoder(8, 0, myFanOutFunction, corruptStream, 3)
    mostLikelyPathsF, scoreVectorF, numberOfEquallyLikelyPathsVectorF, pathsWithFlaggingF, indelFlag, indelEstimatedLocation = viterbiDecoderWithFlagging(8, 0, myFanOutFunction, corruptStream, 3)

    return stream, encodedStream, flatStream, corruptStream, mostLikelyPathsF, scoreVectorF, numberOfEquallyLikelyPathsVectorF, pathsWithFlaggingF, indelFlag, indelEstimatedLocation


def testViterbiDeletion():
    states = [0,1,2,3,4,5,6,7]
    # Omer Sella: triggers are the raw data that we want to encode, post processing i.e.: chopped into blocks that the FSM likes (2 bits in our case).
    triggers = [[0,0], [0,1], [1,0], [1,1]]
    nextStateTable = np.array([[0,1,2,3], [4,5,6,7], [1,0,3,2], [5,4,7,6], [2,3,0,1], [6,7,4,5] , [3,2,1,0], [7,6,5,4] ])
    outputTable = [[[0,0,0], [1,0,0], [0,1,0], [1,1,0]], 
                   [[0,0,1], [1,0,1], [0,1,1], [1,1,1]],
                   [[1,0,0], [0,0,0], [1,1,0], [0,1,0]],
                   [[1,0,1], [0,0,1], [1,1,1], [0,1,1]],
                   [[0,1,0], [1,1,0], [0,0,0], [1,0,0]],
                   [[0,1,1], [1,1,1], [0,0,1], [1,0,1]],
                   [[1,1,0], [0,1,0], [1,0,0], [0,0,0]],
                   [[1,1,1], [0,1,1], [1,0,1], [0,0,1]]]
    initialState = 0
    myFSM = FSM(states, triggers, outputTable, nextStateTable, initialState)
    stream = LOCAL_PRNG.randint(0,2,20)
    encodedStream = FSMEncoder(stream, myFSM)

    # Omer Sella: flatStream gives you the un-chopped encoded stream
    flatStream = []
    for sublist in encodedStream:
        for item in sublist:
            flatStream.append(item)

    ############# Omer Sella: Now we delete some bits and pad with 0s at the end
    corruptStream = np.zeros(len(flatStream), dtype = np.int32)
    
    deletionLocation = 6
    numberOfBitsDeleted = 2
    corruptStream[0 : deletionLocation] =  flatStream[0 : deletionLocation]
    corruptStream[deletionLocation : len(flatStream)] = np.roll(flatStream[deletionLocation : ],  - numberOfBitsDeleted)
    
    corruptStream[len(flatStream) - numberOfBitsDeleted : len(flatStream)] = 0
    print(len(flatStream))
    print(len(corruptStream))
    assert(len(corruptStream) == len(flatStream))

    def myFanOutFunction(state, observation, time):
        return genericFanOutFunction(myFSM, state, observation, time, None)

    #mostLikelyPathsF, scoreVectorF, numberOfEquallyLikelyPathsVectorF, pathsWithFlaggingF, indelFlag, indelEstimatedLocation = viterbiDecoderWithFlagging(8, 0, myFanOutFunction, flatStream, 3, produceGraphics = True)
    # Omer Sella: the line below is WITH a deletion sent to viterbi
    mostLikelyPathsF, scoreVectorF, numberOfEquallyLikelyPathsVectorF, pathsWithFlaggingF, indelFlag, indelEstimatedLocation = viterbiDecoderWithFlagging(8, 0, myFanOutFunction, corruptStream, 3, produceGraphics = False)

    return stream, encodedStream, flatStream, corruptStream, mostLikelyPathsF, scoreVectorF, numberOfEquallyLikelyPathsVectorF, pathsWithFlaggingF, indelFlag, indelEstimatedLocation

def exampleOneThirdConvolutional():
    #Convolutional code taken from here: https://www.researchgate.net/publication/235958269_The_Viterbi_algorithm_demystified/figures?lo=1
    states = [0,1,2,3,4,5,6,7]
    # Omer Sella: triggers are the raw data that we want to encode, post processing i.e.: chopped into blocks that the FSM likes (2 bits in our case).
    triggers = [[0], [1]]
    #Bug: the triggers came out as: [[0, 0], [1, 1], [1, 1], [0, 0], [0, 1]]
    nextStateTable = np.array([[0,4], [0,4], [1,5], [1,5], [2,6], [2,6] , [3,7], [3,7] ])
    outputTable = [[[0,0,0], [1,0,0]], 
                   [[0,1,1], [1,1,1]],
                   [[0,1,0], [1,1,0]],
                   [[0,0,1], [1,0,1]],
                   [[0,0,1], [1,0,1]],
                   [[0,1,0], [1,1,0]],
                   [[0,1,1], [1,1,1]],
                   [[0,0,0], [1,0,0]]]
    initialState = 0
    myFSM = FSM(states, triggers, outputTable, nextStateTable, initialState)
    stream = np.random.randint(0,2,10)
    encodedStream = FSMEncoder(stream, myFSM)
    

    # Omer Sella: flatStream gives you the un-chopped encoded stream
    flatStream = []
    for sublist in encodedStream:
        for item in sublist:
            flatStream.append(item)
    
    ############# Omer Sella: Now we flip a single bit
    corruptStream = copy.deepcopy(flatStream)
    #corruptStream[0] = 1 - corruptStream[0]
    #corruptStream[7] = 1 - corruptStream[0]
    #corruptStream[10] = 1 - corruptStream[0]

    deletionLocation = 11
    numberOfBitsDeleted = 2
    corruptStream[0 : deletionLocation] =  flatStream[0 : deletionLocation]
    corruptStream[deletionLocation : len(flatStream)] = np.roll(flatStream[deletionLocation : ],  - numberOfBitsDeleted)

    def myFanOutFunction(state, observation, time):
        return genericFanOutFunction(myFSM, state, observation, time, None)

    #mostLikelyPath, numberOfEquallyLikelyPaths, paths = viterbiDecoder(8, 0, myFanOutFunction, flatStream, 3)
    mostLikelyPaths, scoreVector, numberOfEquallyLikelyPathsVector, pathsWithFlagging, indelFlag, indelEstimatedLocation = viterbiDecoderWithFlagging(8, 0, myFanOutFunction, corruptStream, 3)

    return stream, encodedStream, flatStream, corruptStream, mostLikelyPaths, scoreVector, numberOfEquallyLikelyPathsVector, pathsWithFlagging, indelFlag, indelEstimatedLocation


def makeEuclidFSM(horizontalSymbols, verticalSymbols, outputFSM, initState = '00000000'):
    euclidFSM = FSM(verticalSymbols, horizontalSymbols, outputFSM, outputFSM, initState)
    return euclidFSM

if __name__ == '__main__':
    #stream, encodedStream, paths, mostLikelyPath, numberOfEquallyLikelyPaths, mostLikelyPaths, scoreVector, numberOfEquallyLikelyPathsVector, pathsWithFlagging, indelFlag, indelEstimatedLocation = exampleTwoThirdsConvolutional()
    #stream, encodedStream, flatEncodedStream, corruptStream, mostLikelyPathsF, scoreVectorF, numberOfEquallyLikelyPathsVectorF, pathsWithFlaggingF, indelFlag, indelEstimatedLocation  = testViterbiBitFlip()
    #stream, encodedStream, mostLikelyPaths, scoreVector, numberOfEquallyLikelyPathsVector, pathsWithFlagging, indelFlag, indelEstimatedLocation = exampleOneHalfConvolutional()
    #stream, encodedStream, flatEncodedStream, corruptStream, mostLikelyPathsF, scoreVectorF, numberOfEquallyLikelyPathsVectorF, pathsWithFlaggingF, indelFlag, indelEstimatedLocation = testViterbiDeletion()
    #stream, encodedStream, flatEncodedStream, corruptStream, mostLikelyPathsF, scoreVector, numberOfEquallyLikelyPathsVector, pathsWithFlagging, indelFlag, indelEstimatedLocation = exampleOneThirdConvolutional()
    streamStr = '010000111011110010001001001000000000'
    initState = '00000000'
    numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols = euclid.testEuclid()
    myFSM = FSM(verticalSymbols, horizontalSymbols, outputFSM, outputFSM, initState)
    es = FSMdictionaryEncoder(streamStr, myFSM)
    print(es)
    
    
