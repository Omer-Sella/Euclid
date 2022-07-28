# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 13:50:45 2022

@author: omers
"""

#import matplotlib.pyplot as plt
import numpy as np


def countGC(sequence):
    # No safety, sequence is assumed to be composed of A,C,T,G
    GC = 0
    AT = 0
    other = 0
    for base in sequence:
        if base == 'G':
            GC = GC + 1
        elif base =='C':
            GC = GC + 1
        elif base == 'A':
            AT = AT + 1
        elif base == 'T':
            AT = AT + 1
        else:
            other = other + 1
    if (GC + AT) == len(sequence):
        status = 'OK'
    else:
        status = 'BAD'
    return status, GC, AT, other
            

def windowedGCcontent(sequence, windowSize = 5):
    #No safety, sequence is assumed to be composed of A,C,T,G
    status, GC, AT, other = countGC(sequence[0 : windowSize])
    slidingPoints = np.arange(windowSize, (len(sequence) - windowSize) , 1)
    gcContent = np.zeros(np.shape(slidingPoints))
    for i in slidingPoints:
        oldBase = sequence[i - 1]
        newBase = sequence[i]
    
        if newBase == 'G' or newBase == 'C':
            GC = GC + 1
        elif newBase == 'T' or newBase == 'A':
            AT = AT + 1
        else:
            other = other + 1
            
        if oldBase == 'G' or oldBase == 'C':
            GC = GC - 1
        elif oldBase == 'T' or oldBase == 'A':
            AT = AT - 1
        else:
            other = other - 1
        
        gcContent[i - windowSize] = GC
    
    gcContent = gcContent / (GC + AT + other)
    #print(gcContent)
    return status, slidingPoints, gcContent, GC, AT, other
            
    
def gcVariance(sequence, windowSize = 50):
    import seaborn as sns
    import matplotlib.pylab as plt
    plt.style.use("seaborn")
    status, slidingPoints, gcContent, GC, AT, other = windowedGCcontent(sequence, windowSize = windowSize)
    varianceArray = np.zeros((len(gcContent), len(gcContent)), dtype = np.float64)
    for i in range(len(gcContent)):
        for j in range(len(gcContent)):
            varianceArray[i,j] = np.abs(gcContent[i] - gcContent[j])
    plt.figure(figsize=(len(gcContent),len(gcContent)))
    heat_map = sns.heatmap( varianceArray, linewidth = 1 , annot = True)
    plt.title( "G-C difference (variance) between stretches of " + str(windowSize) + " nucleotides" )
    plt.show()
    
def windowedGCDistribution(sequence, windowSize):
    pass

    

    
