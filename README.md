# Euclid
Euclid is a way to produce a Finite State Machine (FSM) encoding for DNA data storage, from a list of constraints.
The constraint list is made of regular expressions as well as some numeric parameters. An example for a constraint list that dictates GC levels between 0.35 and 0.65 is:
```python
constraintList = {'gcMin': 0.35, 'gcMax': 0.65}
```
 If we add maximum run length (maximal nuber of identical consecutive bases) of 8 we get:
```python
constraintList = {'gcMin': 0.35, 'gcMax': 0.65, 'runLength': 8}
```
If in addition we wanted to avoid the sequence GNGTC where N is one of {A,C,T,G}, we would add a regular expression:
```python
constraintList = {'gcMin': 0.35, 'gcMax': 0.65, 'runLength': 8,'regex3': 'G[ACTG]GTC'}
```
Once this list is formed, and a symbol size (the size of each symbol in the FSM encoding) is set (you need to choose the symbol size, but to run on a reasonable single machine limit yourself to a maximum of 6 or 7), you can use the function connectivityMatrix:
```python
connectivityMatrix, sequencesBinary, sequencesBases, violationsMatrix, connectivityDictionaryBinary = connectivityMatrix(constraintList, symbolSize)
```
connectivityMatrix, connectivityDictionaryBinary and sequenceBinary could then be used to form two lists of symbols (vertical and horizontal, but consider them just as names) as well as candidates, that could be plugged between a vertical and a horizontal symbol, in a way that will not violate the list of constraints: 
```python
candidates, vSymbols, hSymbols = generateCandidates(matrix = connectivityMatrix, connectivityDictionary = connectivityDictionaryBinary, seqBinary = sequencesBinary, seqBinary = sequencesBinary, ['0', '1'])
```python
Finally, use the candidates, vSymbols and hSymbols, as well as a choice mechanism (trackGC or random are the only ones currently implemented), to commit to a specific choice of symbols (making the encoding completely deterministic):
```python
numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols = makeFSM(c = candidates, vsym = vSymbols, hsym = hSymbols, trackGClevel)
```
Note that trackGClevel tracks 0.5 by default and could be changed.


