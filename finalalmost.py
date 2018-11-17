import numpy as np
import scipy as sp
import random
from scipy import constants
import time
from fractions import Fraction
from math import gcd


#GATES#
X = np.array([[0,1],[1,0]])
I = np.identity(2)
H = (2.0**(-0.5)) * np.array([[1,1],[1,-1]])
CNOT = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
TONC = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]])

#end of Gates






#ground and first excited state qubits
q0 = np.array([[1],[0]])
q1 = np.array([[0],[1]])




def performMeasurement(V):
    P = V
    numOfBits = len(P)
    realNum = int(np.log2(numOfBits))
    #print(realNum)
    for i in range(numOfBits):
        P[i] = abs(V[i])**2
    
    found = False
    R = random.uniform(0, 1)
    measured = -1
    total = 0
    for i in range(numOfBits):
        total += P[i]
        if R <= total and found == False:
            found = True
            measured = i
            
    print("|"+bin(measured)[2:].zfill(realNum)+">")
    
    
def performMeasurement2(V):
    print('Measuring Quantum State')
    time.sleep(2)
    P = V
    numOfBits = len(P)
    realNum = int(np.log2(numOfBits))
    #print(realNum)
    for i in range(numOfBits):
        P[i] = abs(V[i])**2
    
    found = False
    R = random.uniform(0, 1)
    measured = -1
    total = 0
    for i in range(numOfBits):
        total += P[i]
        if R <= total and found == False:
            found = True
            measured = i
            
    print("|"+bin(measured)[2:].zfill(realNum)+">")
    return (bin(measured)[2:].zfill(realNum))
    
 

def GroverOnTwo():
    
    state0=np.kron(q0,q0)
    state1=np.dot(np.kron(H,H),state0)
    print("Grover's Algorithm is a search algorithmn that can search an unordered database in O(sqrt(N)) time.")
    
    userResponse = int(input("Would you like to hide the solution in the 0th,1st,2nd, or 3rd position? "))
    userOracle = I
    if userResponse == 0:
        oracle0 = MatrixMult(np.kron(X,I), np.kron(H,X), TONC, np.kron(H,X), np.kron(X,I))
        userOracle = oracle0
    elif userResponse == 1:
        oracle1 = MatrixMult(np.kron(X,H), CNOT, np.kron(X,H))
        userOracle = oracle1
    elif userResponse == 2:
        oracle2 = MatrixMult(np.kron(H,X), TONC, np.kron(H,X))
        userOracle = oracle2
    elif userResponse == 3:
        oracle3 = MatrixMult(np.kron(I,H), CNOT, np.kron(I,H))
        userOracle = oracle3
    
    #relfectionOverS = MatrixMult(np.kron(H,H), np.kron(I,H), CNOT, np.kron(I,H), np.kron(X,H), CNOT, np.kron(X,H), np.kron(H,X), TONC, np.kron(H,X), np.kron(H,H)) 
    other = MatrixMult(np.kron(H,H), np.kron(X,I), np.kron(H,X), TONC, np.kron(H,X), np.kron(X,I), np.kron(H,H))
    state2 = MatrixMult(userOracle, state1)
    state3 = np.dot(other,state2)
    performMeasurement(state3)

def Joyza():
    
    print("The Deutsch-Joyza algorithm generalizes the Deutsch algorithm to any number of inputs.")
    print("The number of inputs will equal the number of qubits.")
    print("Due to the complexity of these functions, we only present three algorithms")
    
    numOfQubits = int(input('How many qubits would you like to run?'))
    state0 = np.kron(q0,q1)
    for i in range(numOfQubits-2):
        state0 = np.kron(q0,state0)
    
    print("Would you like to run a:")
    print("     1: Constant-1 function")
    print("     2: Constant-0 function")
    print("     3: A particular balanced function")
    
    
    userResponse = int(input("Please enter the number corresponding to a function: "))
    userOracle = I
    #Constant 0 function
    if userResponse == 2:
        deutschJoyzaOracleConst0 = np.kron(I,I)
        for i in range(numOfQubits-2):
            deutschJoyzaOracleConst0 = np.kron(I,deutschJoyzaOracleConst0)
        userOracle = deutschJoyzaOracleConst0
        
    #Constant 1 function
    if userResponse == 1:
        deutschJoyzaOracleConst1 = np.kron(I,X)
        for i in range(numOfQubits-2):
            deutschJoyzaOracleConst1 = np.kron(I,deutschJoyzaOracleConst1)
        userOracle = deutschJoyzaOracleConst1
        
    #Balanced function
    if userResponse == 3:
        deutschJoyzaOracleBalanced = CNOT
        for i in range(numOfQubits-2):
            deutschJoyzaOracleBalanced = np.kron(I,deutschJoyzaOracleBalanced)
        userOracle = deutschJoyzaOracleBalanced
        
    manyHadamard = np.kron(H,H)
    for i in range(numOfQubits-2):
        manyHadamard = np.kron(H,manyHadamard)
        
    state1 = np.dot(manyHadamard,state0)
    state2 = np.dot(userOracle,state1)
    state3 = np.dot(manyHadamard,state2)
    
    performMeasurement(state3)
    return

def MatrixMult(*args):
    lenOfArgs = len(args)
    #print(lenOfArgs)
    intial = args[0]
    current = intial
    for i in range(1,lenOfArgs):
        current = np.dot(current,args[i])
        
    return current






def Deutsch():
    
    deutschOracleConst0 = np.identity(4)
    deutschOracleConst1 = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
    deutschOracleId = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    deutschOracleNegId = np.array([[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]])
    
    print("The Deutsch Algorithm tells you if a function that accepts either a 0 or a 1 is constant or not.")
    print("Would you like to run")
    print("     1: The Constant-1 function")
    print("     2: The Identity function")
    print("     3: The Negation function") 
    print("     4: The Constant-0 function")
    
    
    
    badInput = True
    while(badInput):
        userResponse = int(input("Enter the number corresponding to the function you would like to run: "))
        if userResponse == 1:
            userOracle = deutschOracleConst1
            badInput = False
        elif userResponse == 2:
            userOracle = deutschOracleId
            badInput = False
        elif userResponse == 3:
            userOracle = deutschOracleNegId
            badInput = False
        elif userResponse == 4:
            userOracle = deutschOracleConst0
            badInput = False
        else:
            print("Please enter a number between 1 and 4")
    state0=np.kron(q0,q1)
    state1=np.dot(np.kron(H,H),state0)
    state2=np.dot(userOracle,state1)
    state3=np.dot(np.kron(H,H),state2)
    
    performMeasurement(state3)
    return



def QuantumShor(N,x,t,L):
    time.sleep(2)
    states = q1
    states = np.array(states)
    
    for i in range(L-1):
        states = np.kron(q1,states)
        
    for i in range(t):
        states = np.kron(np.dot(H,q0),states)
        
    #unitary matrix section
    
    print('Applying Unitary Gates')
    
    j = states[:2**t]*(2**(0.5*t))
    
    j = []
    
    for i in range(2**t):
        j.append(np.binary_repr(i,t))
        
    mod = []
    
    for i in range(2**t):
        temp = 0
        for k in range(len(j[i])):
            if (int(j[i][k]) == 1):
                temp += 2**(k)
        mod.append(np.binary_repr((x**temp)%N,L))
        
    total = [''] * 2**(t)
    
    for i in range(len(j)):
        total[i] = str(j[i]) + str(mod[i])
    
    states = np.zeros(2**(t+L))
    
    for i in range(len(total)):
        tempstate = []
        if (int(total[i][0]) == 0):
            tempstate = q0
        elif (int(total[i][0]) == 1):
            tempstate = q1
            
        for k in range(1,len(total[i])):
            if (int(total[i][k]) == 0):
                tempstate = np.kron(tempstate,q0)
            elif (int(total[i][k]) == 1):
                tempstate = np.kron(tempstate,q1)
                
        states = states + tempstate
        
        
    states = states/(2**(0.5*t))
    
    
    #Inverse Quantum Fourier Transform Section
    
    print('Applying Inverse Quantum Fourier Transform')
    
    qft = []
    for i in range(2**(t+L)):
        temp = []
        for j in range(2**(t+L)):
            temp.append(np.exp(2j*sp.constants.pi/(L+t))**(i*j))
        qft.append(temp)
        
    qft = np.array(qft)
    
    iqft = np.linalg.inv(qft)
    
    states = np.dot(iqft,states)
    
    states = states[0]
    
    return states

def continued_fraction(s,n,x,N):
    sr = Fraction(s/n)
    
    b = True
    
    if (sr == 0):
        return N
    
    while b:
        if (sr > 1):
            sr = 1 - sr
            if ((x**(sr.denominator))%N==1):
                return sr.denominator
            else:
                temp1 = sr.numerator
                temp2 = sr.denominator
                sr = Fraction(temp1,temp2)
        if(s == 1):
            b = False
    
    return N
    

def Shor(N,error=0.5):
    print('Creating Quantum State')
    
    if(N % 2 == 0):
        return 2

    L = len(np.binary_repr(N))
    t = 2*L + 1 + np.ceil(np.log2(2+1/(2*error)))
    t = int(t)
        
    x = 1
    
    i = 3
    
    #This should be random at some point    
    while i < N:
        x = gcd(i,N)
        if (x == 1):
            x = i
            s = performMeasurement2(QuantumShor(N,x,t,L))
            sr = 0
            for k in range(len(s)):
                if (int(s[k]) == 1):
                    sr += 2**(k)
            print('Analyzing Measuerement')
            r = continued_fraction(sr,2**(t+L),x,N)
    
            if (r % 2 == 0):
                a = (x**r)%N-1
                b = (x**r)%N+1
                if (N % a == 0):
                    return a
                elif (N % b == 0):
                    return b
        i = i + 1
        print('Randomness Failure: Retrying')
    
    return N
    
def getFactors():
    N = int(input('What number would you like to factor? '))
    if (N == 1): 
        print(N)
        return
    factors = []
    
    while N > 1:
        f = Shor(N)
        factors.append(f)
        N = int(N / f)
        
    print('Analyzing Factors')
    time.sleep(2)
        
    for i in range(len(factors)):
        print(factors[i])
        
    return
    
def customCircuit():
    print("This is place where you can design and run your own custom quantum circuits.")
    
    numOfQubits = int(input("How many qubits would you like to run?"))
    anotherLayer = True
    listOfGates = []
    
    #calculate state 0
    state0 = np.kron(q0,q0)
    for i in range(numOfQubits-2):
        state0 = np.kron(q0,state0)
    currentState = state0
    
    
    tempGate = I
    userGate = "I"
    while(anotherLayer):
        i = 0
        kronOfGates = I
        listOfGates = []
        while(i < numOfQubits):
            userGate = input("What gate would you like on qubit " + str(i) + "? ")
            
            if userGate == 'I':
                tempGate = I
            elif userGate == 'X':
                tempGate = X
            elif userGate == "H":
                tempGate = H
            elif userGate == 'CNOT':
                tempGate = CNOT
                i += 1   
            elif userGate == 'TONC':
                tempGate = TONC
                i += 1
            else:
                tempGate = H
            i += 1
            listOfGates.append(tempGate)
        
        
        kronOfGates = listOfGates[0]
        for i in range(1, len(listOfGates)):
            kronOfGates = np.kron(kronOfGates,listOfGates[i])

        currentState = np.dot(kronOfGates,currentState)
        
        seguir = input("Would you like to add another layer of gates y/n? ")
        if seguir == 'n':
            anotherLayer = False
            
    for i in range(10):
        stateToBeMeasured = currentState.copy()
        performMeasurement(stateToBeMeasured)


def initializeTheTerminal():
    print("Welcome to the USU Quantum Terminal!")
    print("Would you like to: ")
    print("     1: Run the Deutsch Algorithm")
    print("     2: Run the Deutsch-Joyza Algorithm")
    print("     3: Run Grover's Algorithm")
    print("     4: Run Shor's Factoring Algorithm")
    print("     5: Design your own quantum circuit")
    
    
    #while input is trash
    userAnswer = int(input("Enter the corresponding number: "))
    if userAnswer == 1:
        Deutsch()
    elif userAnswer == 2:
        Joyza()
    elif userAnswer == 3:
        GroverOnTwo()
    elif userAnswer == 4:
        getFactors()
    elif userAnswer == 5:
        customCircuit()
    else:
        print("Please enter a valid number")
        
    input()
    return



initializeTheTerminal()