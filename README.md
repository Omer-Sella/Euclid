# Euclid
Euclid is a way to produce a Finite State Machine (FSM) encoding for DNA data storage, from a list of constraints.
The constraint list is made of regular expressions as well as some numeric parameters. An example for a constraint list that dictates GC levels between 0.35 and 0.65 is:
 constraintList = {'gcMin': 0.35, 'gcMax': 0.65}
 If we add maximum run length (maximal nuber of identical consecutive bases) of 8 we get:
  constraintList = {'gcMin': 0.35, 'gcMax': 0.65, 'runLength': 8}
If in addition we wanted to avoid the sequence GNGTC where N is one of {A,C,T,G}, we would add a regular expression:
 constraintList = {'gcMin': 0.35, 'gcMax': 0.65, 'runLength': 8,'regex3': 'G[ACTG]GTC'}
Once this list is formed, and a symbol size (the size of each symbol in the FSM encoding) is set, you can use the function connectivityMatrix to get 
connectivityMatrix, sequencesBinary, sequencesBases, violationsMatrix, connectivityDictionaryBinary = connectivityMatrix(constraintList, symbolSize)
From 
