# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 16:27:50 2023

@author: Matteo De Pascali
"""

import numpy as np
import scipy.io
import re
import sys

def smartFreqResp(A,B,C,D,w):
    jw = 1j*w
    I = np.identity(len(A))
    G = scipy.linalg.solve((jw*I-A), B)
    fr = np.matmul(C,G) + D
    return fr[0,0]

def matrixExtractor(M,nx):
    A = M[0:nx,0:nx]
    B = M[0:nx,nx:]
    C = M[nx:,0:nx]
    D = M[nx:,nx:]
    return A,B,C,D

# Flags
debugProcedure = False

# Linearization file name
linFile = "dslinSCO2.mat"

# Settings for the analysis
# Possible inputs: wHot and wCold
# Possible outputd: TIT and Power
selectedInputName = 'wHot'
selectedOutputName = 'Power'

# Frequency at which the analysis is performed
w = 1/300  # hot --> 1/300   cold --> 1/30

#Importing data
matFile = scipy.io.loadmat(linFile)

# Variables definition
nx = matFile["nx"][0][0]
ABCD = np.matrix(matFile["ABCD"])
statesNames = matFile["xuyName"][0:nx]
A,B,C,D = matrixExtractor(ABCD,nx)

# "Cleaning" strings and removing indexes of vector variables
for i in range(len(statesNames)):
    statesNames[i] = statesNames[i].rstrip('\x00   ')
    statesNames[i] = re.sub(r'\[.*?\]', '', statesNames[i])
    
# Creating dictionary of unique states and their indeces
uniqueStates = dict()
maxSizeNameField = 0
for i in range(len(statesNames)):
    key = statesNames[i]
    if key not in uniqueStates.keys():
        uniqueStates[key] = [i]
    else:
        uniqueStates[key].append(i)
        if len(key) > maxSizeNameField:
            maxSizeNameField = len(key) + 4

# Normalizng A matrix
# Creating normalization matrix
normVect = []
for elem in statesNames:
    m = re.search(r'[^\.]([^.]*)$',elem)
    if m.group(0) == "ptilde" or m.group(0) == "pout" or m.group(0) == "pOut":
        normVect.append(100000)
    elif m.group(0) == "T" or m.group(0) == "Ttilde" or m.group(0) == "Tvol":
        normVect.append(1000)
    else:
        normVect.append(1)

N = np.diag(normVect)

# Normalizing matrices
A = np.linalg.inv(N)*A*N
B = np.linalg.inv(N)*B
C = C*N
D = D

# String formatting, header and separator
stringTableFormat = "{:^1} {:^" + str(maxSizeNameField) + "} {:^1} {:^11} {:^1}"
header = stringTableFormat.format('||','State','||','Error [%]','||')
separator = "=" * len(header)
    
# Select input-output pair and cutting system matrices for SISO representation
listOfInputs = {"wHot": 0, "wCold": 1}
listOfOutputs = {"TIT": 0, "Power": 1}
selectedInput = listOfInputs[selectedInputName]
selectedOutput = listOfOutputs[selectedOutputName]

Asel = A
Bsel = B[:,selectedInput] # select column related to selected input
Csel = C[selectedOutput,:] # select row related to selected output
Dsel = D[selectedOutput,selectedInput] # same for D matrix

# Check if just one state defines the output selected
if len(np.nonzero(Csel)) == 1:
    print(f"Output {selectedOutput} corresponds to a single state")
elif len(np.nonzero((Csel))) == 0:
    print(f"Output {selectedOutput} can't be observed!")
    sys.exit()
    
if len(np.nonzero(Bsel)) == 1:
    print(f"Output {selectedOutput} corresponds to a single state")
elif len(np.nonzero((Bsel))) == 0:
    print(f"Output {selectedOutput} can't be controlled!")
    sys.exit()

# Complete system response
Hfull = smartFreqResp(Asel, Bsel, Csel, Dsel, w)

# %%
print(f"The selected linearized system file is {linFile}")
print()

print(f"The selected input for the analysis is {selectedInputName}")
print(f"The selected input for the analysis is {selectedOutputName}")

print()
print(f"Frequency: {w} rad/s")

# Evaluating "reduced states set" frequency response
globalIndices = []
removedStates = []
listOfErrors = []
i = 0
while len(removedStates)<(len(uniqueStates.keys())-1):
    errorDict = dict()
    nullFR = []
    for elem in uniqueStates.keys():
        if elem not in removedStates:
            indices = uniqueStates[elem] + globalIndices
            Atemp1 = np.delete(Asel, indices, axis=0)
            Atemp = np.delete(Atemp1, indices, axis=1)
            Btemp = np.delete(Bsel, indices, axis=0)
            if len(np.nonzero(Btemp)[0]) == 0:
                nullFR.append(elem)
                pass
            Ctemp = np.delete(Csel, indices, axis=1)
            if len(np.nonzero(Ctemp)[0]) == 0:
                nullFR.append(elem)
                pass
            Dtemp = Dsel
            H = smartFreqResp(Atemp, Btemp, Ctemp, Dtemp, w)
            errorDict[elem] = abs(np.divide(Hfull-H,Hfull))*100
        else:
            continue
    
    errorDictList = sorted(errorDict, key=lambda dict_key: abs(errorDict[dict_key]))
    if (errorDictList[0] in nullFR and debugProcedure):
        print("Error! Either B or C matrices will contain just zeros in the next iteration!")
    removedStates.append(str(errorDictList[0]))
    listOfErrors.append(errorDict[errorDictList[0]])
    globalIndices += uniqueStates[errorDictList[0]]
    
    # Printing for debug
    if debugProcedure:
        if i == 0:
            print("Debugging: showing each iteration of the first part of the algorithm!")
        print(separator)
        print(header)
        print(separator)
        for elem in errorDictList:
            print(stringTableFormat.format('||',f"{elem}",'||',f"{errorDict[elem]:.2f}",'||'))
        print(separator)
        
    i+=1
if debugProcedure:
    print("\nTermination!  \n\n\n")
    
# Printing whole sequence of data
if debugProcedure:
    print("Summary of the debugging procedure of the first part of the algorithm below\n")
    print(separator)
    print(header)
    print(separator)
    
    for i in range(len(removedStates)):
        print(stringTableFormat.format('||',f"{removedStates[i]}",'||',f"{listOfErrors[i]:.2f}",'||'))
        
    print(stringTableFormat.format('||',f"{errorDictList[-1]}",'||',"-",'||'))
    print(separator)

# Writing final output
fullListOfOrderedStates = removedStates + [str(errorDictList[-1])]
fullListOfOrderedStates.reverse()
listOfErrors.reverse()
listOfErrors.append(0)
 
print("\n\nOutput of the first part of the algorithm:\n")

print(separator)
print(header)
print(separator)
for i in range(len(fullListOfOrderedStates)):
    print(stringTableFormat.format('||',f"{fullListOfOrderedStates[i]}",'||',f"{listOfErrors[i]:.2f}",'||'))
print(separator)


#%%

def updateIndices(statesToKeepList, uniqueStatesDict):
    # Update indices of the remaining unique states after elimination
    sortedUniqueKeys = sorted(uniqueStatesDict, key=lambda x: uniqueStatesDict[x][0])
    newUniqueStatesDict = dict()
    shiftLen = 0
    for state in sortedUniqueKeys:
        if state in statesToKeepList:
              newUniqueStatesDict[state] = [x-shiftLen for x in uniqueStatesDict[state]]
        else:
            shiftLen += len(uniqueStatesDict[state])
    return newUniqueStatesDict

def reduceSystem(A, B, C, D, statesToKeepList, uniqueStatesDict):
    # Indices of the states to keep form the previous analysis
    redIndices = []
    for state in statesToKeepList:    
        redIndices += uniqueStatesDict[state]
    # Composing the matrices of the reduced state from previous analysis
    Ared = A[np.ix_(sorted(redIndices),sorted(redIndices))]
    Bred = B[sorted(redIndices), :]
    Cred = C[:, sorted(redIndices)]
    Dred = D
    # Update inidices of the remaining unique states after elimination
    newUniqueStatesDict = updateIndices(statesToKeepList, uniqueStatesDict)
    return Ared, Bred, Cred, Dred, newUniqueStatesDict

def instantStatesReduction(A,B,C,D,selState,uniqueStatesDict):
    # Indices of the instantaneous states
    selIndices = sorted(uniqueStatesDict[selState])
    # Indices of the states that are conserved
    remainingIndices = sorted([x for elem in list(uniqueStatesDict.keys()) if elem != selState for x in uniqueStatesDict[elem]])

    # Creating the matrices that transfer the effect of the instantaneous states on the other equations in the state space representation
    # Matrix of the terms in the equation of the conserved states referred to the instantaneous states
    addA1 = A[np.ix_(remainingIndices, selIndices)]
    # Inverse of the matrix of the terms of the instantaneous states in the equations of the instantaneous states
    addA2 = -np.linalg.inv(A[np.ix_(selIndices, selIndices)])
    # Matrix of the termes referred to the conserved states in the instantaneous states equations
    addA3 = A[np.ix_(selIndices, remainingIndices)]
    # Matrix of the terms referred to input in the equations of the instantaneous states
    addB1 = B[np.ix_(selIndices, [0])]
    # Matrix referred to the terms of the instananeous states in the output equation 
    addC1 = C[np.ix_([0], selIndices)]

    # Building matrices that need to be added to the matrices of the conserved states
    addA = np.matmul(np.matmul(addA1, addA2), addA3)
    addB = np.matmul(np.matmul(addA1, addA2), addB1)
    addC = np.matmul(np.matmul(addC1, addA2), addA3)
    addD = np.matmul(np.matmul(addC1, addA2), addB1)

    # New system with instantaneous states
    Ared = A[np.ix_(remainingIndices,remainingIndices)] + addA
    Bred = B[np.ix_(remainingIndices, [0])] + addB
    Cred = C[np.ix_([0],remainingIndices)] + addC
    Dred = D + addD
    
    # Update dictionary
    updatedListOfStatesToKeep = list(uniqueStatesDict.keys())
    updatedListOfStatesToKeep.remove(selState)
    updatedDict = updateIndices(updatedListOfStatesToKeep, uniqueStatesDict)
    
    return Ared, Bred, Cred, Dred, updatedDict


def groupStatesInMatrices(A,B,C,D,uniqueStatesDict):
    shift = 0
    permVect = []
    newDict = {}
    for elem in uniqueStatesDict.keys():
        permVect = permVect + uniqueStatesDict[elem]
        newDict[elem] = list(range(shift, shift+len(uniqueStatesDict[elem])))
        shift = shift + len(uniqueStatesDict[elem])
        
    # permVect = [x for elem in list(uniqueStatesDict.keys()) for x in uniqueStatesDict[elem]]
    Aord = A[np.ix_(permVect,permVect)]
    Bord = B[np.ix_(permVect,[0])]
    Cord = C[np.ix_([0],permVect)]
    Dord = D
    return Aord, Bord, Cord, Dord, newDict


statesToKeep = ["hX.coldside.ptilde", "hX.wall.Tvol", "hX.coldside.Ttilde", "hX.hotside.Ttilde"]  # hot simple case

Aord,Bord,Cord,Dord,uniqueStatesOrd = groupStatesInMatrices(Asel,Bsel,Csel,Dsel,uniqueStates)
Ared, Bred, Cred, Dred, uniqueStatesRed = reduceSystem(Aord,Bord,Cord,Dord,statesToKeep,uniqueStatesOrd)
Hred = smartFreqResp(Ared, Bred, Cred, Dred, w)
ErrorRed = abs(np.divide(Hfull-Hred,Hfull))*100

statesToIterate = statesToKeep.copy()
Ainst = Ared.copy()
Binst = Bred.copy()
Cinst = Cred.copy()
Dinst = Dred.copy()
dictInst = uniqueStatesRed.copy()
instStates = []
errorInstStates = []
while len(statesToIterate) > 0:
    errorsDict = dict()
    for state in statesToIterate:
        AinstTemp,BinstTemp,CinstTemp,DinstTemp,tempDict = instantStatesReduction(Ainst,Binst,Cinst,Dinst,state,dictInst)
        HinstTemp = smartFreqResp(AinstTemp,BinstTemp,CinstTemp,DinstTemp, w)
        errorInstTemp = abs(np.divide(Hfull-HinstTemp,Hfull))*100
        errorsDict[state] = errorInstTemp
        
    errorDictList = sorted(errorsDict, key=lambda dict_key: abs(errorsDict[dict_key]))
    stateToInst = errorDictList[0]
    Ainst,Binst,Cinst,Dinst,dictInst = instantStatesReduction(Ainst,Binst,Cinst,Dinst,stateToInst,dictInst)
    statesToIterate.remove(stateToInst)
    instStates.append(stateToInst)
    errorInstStates.append(errorsDict[stateToInst])
    

instStates.reverse()
errorInstStates.reverse()
errorInstStates.append(ErrorRed)


print("\n\nOutput of the second part of the algorithm after selecting\n")
for elem in statesToKeep:
    print(elem)
print("\nas states to keep after the first part of the algorithm\n")

print(separator)
print(header)
print(separator)
for i in range(len(instStates)):
    print(stringTableFormat.format('||',f"{instStates[i]}",'||',f"{errorInstStates[i+1]:.2f}",'||'))
print(separator)
