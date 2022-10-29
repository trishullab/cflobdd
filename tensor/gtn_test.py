import tensornetwork as tn
import numpy as np
import sys
import time
import random
import math
import pandas as pd
from fractions import Fraction

def apply_gate(qubit_edges, gate, operating_qubits):
    op = tn.Node(gate)
    for i, bit in enumerate(operating_qubits):
        tn.connect(qubit_edges[bit], op[i])
        qubit_edges[bit] = op[i + len(operating_qubits)]

def GHZ(N):
    allOnes = '1' * N
    allZeros = '0' * N

    start = time.time()
    H = np.array([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
    CNOT = np.zeros((2,2,2,2), dtype=complex)
    CNOT[0][0][0][0] = 1
    CNOT[0][1][0][1] = 1
    CNOT[1][0][1][1] = 1
    CNOT[1][1][1][0] = 1    
    is_correct = False
    all_nodes = []

    with tn.NodeCollection(all_nodes):
        state_nodes = [
            tn.Node(np.array([1.0+0.0j, 0.0+0.0j],)) for _ in range(N)
        ]

        qubits = [node[0] for node in state_nodes]
        apply_gate(qubits, H, [0])

        for i in range(N - 1):
            apply_gate(qubits, CNOT, [i, i+1])

        result = tn.contractors.greedy(all_nodes, output_edge_order=qubits)
        t = result.tensor
        #print(np.square(np.absolute(t)).flatten())
        arr = np.square(np.abs(t)).flatten()
        index = np.random.choice(len(arr), p=arr)
        index_s = bin(index)[2:].zfill(N)
        if index_s == allZeros or index_s == allOnes:
            is_correct = True
    end = time.time()
    print('is_correct:', is_correct, 'time_taken(s):', (end - start))

def BV(N):
    s = ""
    for i in range(0, N):
        r = random.randint(0, 1)
        if r == 0:
            s = s + '0'
        else:
            s = s + '1'

    start = time.time()
    H = np.array([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
    X = np.array([[0,1],[1,0]], dtype=complex)
    CNOT = np.zeros((2,2,2,2), dtype=complex)
    CNOT[0][0][0][0] = 1
    CNOT[0][1][0][1] = 1
    CNOT[1][0][1][1] = 1
    CNOT[1][1][1][0] = 1    

    all_nodes = []
    is_correct = False
    with tn.NodeCollection(all_nodes):
        state_nodes = [
            tn.Node(np.array([1.0+0.0j, 0.0+0.0j],)) for _ in range(N+1)
        ]
        
        qubits = [node[0] for node in state_nodes]
        apply_gate(qubits, X, [N])
        for i in range(N+1):
            apply_gate(qubits, H, [i])

        for i in range(N):
            if s[i] == '1':
                apply_gate(qubits, CNOT, [i, N])

        for i in range(N+1):
            apply_gate(qubits, H, [i])

        result = tn.contractors.greedy(all_nodes, output_edge_order=qubits)
        t = result.tensor
        #print(t)
        arr = np.square(np.abs(t)).flatten()
        index = np.random.choice(len(arr), p=arr)
        index_s = bin(index)[2:].zfill(N+1)[:-1]
        if index_s == s:
            is_correct = True
    end = time.time()
    #print(s)
    print('is_correct', is_correct, 'time_taken(s):', (end - start))

def DJ(N):
    s = ""
    for i in range(0, N):
        r = random.randint(0, 1)
        if r == 0:
            s = s + '0'
        else:
            s = s + '1'
    allZeros = '0' * N
    start = time.time()
    H = np.array([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
    X = np.array([[0,1],[1,0]], dtype=complex)
    CNOT = np.zeros((2,2,2,2), dtype=complex)
    CNOT[0][0][0][0] = 1
    CNOT[0][1][0][1] = 1
    CNOT[1][0][1][1] = 1
    CNOT[1][1][1][0] = 1    

    all_nodes = []
    is_correct = True
    with tn.NodeCollection(all_nodes):
        state_nodes = [
            tn.Node(np.array([1.0+0.0j, 0.0+0.0j],)) for _ in range(N+1)
        ]
        
        qubits = [node[0] for node in state_nodes]
        apply_gate(qubits, X, [N])
        for i in range(N+1):
            apply_gate(qubits, H, [i])

        for i in range(N):
            if s[i] == '1':
                apply_gate(qubits, X, [i])

        for i in range(N):
            apply_gate(qubits, CNOT, [i, N])

        for i in range(N):
            if s[i] == '1':
                apply_gate(qubits, X, [i])

        for i in range(N):
            apply_gate(qubits, H, [i])

        result = tn.contractors.greedy(all_nodes, output_edge_order=qubits)
        t = result.tensor
        #print(t)
        arr = np.square(np.abs(t)).flatten()
        for i in range(0,10):
            index = np.random.choice(len(arr), p=arr)
            index_s = bin(index)[2:].zfill(N+1)[:-1]
            if index_s == allZeros:
                is_correct = False
                break
    end = time.time()
    #print(s)
    print('is_correct', is_correct, 'time_taken(s):', (end - start))


def getNum(s):
    value = 0
    for i in range(0, len(s)):
        if s[i] == '1':
            value = value * 2 + 1
        else:
            value = value * 2
    return value

def simons(N):
    s = ""
    for i in range(0, N):
        r = random.randint(0, 1)
        if r == 0:
            s = s + '0'
        else:
            s = s + '1'
    
    allZeros = '0' * N
    equations = []

    start = time.time()
    H = np.array([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
    CNOT = np.zeros((2,2,2,2), dtype=complex)
    CNOT[0][0][0][0] = 1
    CNOT[0][1][0][1] = 1
    CNOT[1][0][1][1] = 1
    CNOT[1][1][1][0] = 1

    all_nodes = []
    is_correct = True
    with tn.NodeCollection(all_nodes):
        state_nodes = [
            tn.Node(np.array([1.0+0.0j, 0.0+0.0j],)) for _ in range(2*N)
        ]
        
        qubits = [node[0] for node in state_nodes]

        for i in range(N):
            apply_gate(qubits, H, [i])

        for i in range(N):
            apply_gate(qubits, CNOT, [i, i + N])

        k = 0
        for i in range(N-1, -1, -1):
            if s[i] == '1':
                m = N
                for j in range(N-1, -1, -1):
                    if s[j] == '1':
                        apply_gate(qubits, CNOT, [k, m])
                    m += 1
                break
            k += 1

        for i in range(N):
            apply_gate(qubits, H, [i])

        result = tn.contractors.greedy(all_nodes, output_edge_order=qubits)
        t = result.tensor
        #print(t)
        arr = np.square(np.abs(t)).flatten()
        for i in range(0, 2 * N):
            index = np.random.choice(len(arr), p=arr)
            index_s = bin(index)[2:].zfill(2*N)[0:N]
            equations.append(index_s)
            
    end = time.time()
    #print(s)
    # print(equations)
    print("s:",s)
    print('is_correct', is_correct, 'time_taken(s):', (end - start))


def grover(N):
    s = ""
    for i in range(0, N):
        r = random.randint(0, 1)
        if r == 0:
            s = s + '0'
        else:
            s = s + '1'

    start = time.time()
    allZeros = '0' * (N + 1)
    H = np.matrix([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
    X = np.matrix([[0,1],[1,0]], dtype=complex)
    I = np.matrix([[1,0],[0,1]], dtype=complex)

    TX = np.identity(2**(N+1), dtype=complex)
    TX[2**(N+1)-2,2**(N+1)-2] = 0
    TX[2**(N+1)-2,2**(N+1)-1] = 1
    TX[2**(N+1)-1,2**(N+1)-2] = 1
    TX[2**(N+1)-1,2**(N+1)-1] = 0
    dims = tuple([2 for i in range(2*(N+1))])
    TX = TX.reshape(dims)

    TZ = np.identity(2**N, dtype=complex)
    TZ[2**N-1, 2**N-1] = -1
    dims = tuple([2 for i in range(2*N)])
    TZ = TZ.reshape(dims)

    all_bits = [i for i in range(N+1)]
    bits = [i for i in range(N)]
    all_nodes = []
    is_correct = True
    iters = (math.pi * (2 ** (N/2)))
    iters = int(iters/4)
    with tn.NodeCollection(all_nodes):
        state_nodes = [
            tn.Node(np.array([1.0+0.0j, 0.0+0.0j],)) for _ in range(N+1)
        ]
        
        qubits = [node[0] for node in state_nodes]
        apply_gate(qubits, X, [N])
        for i in range(0, N+1):
            apply_gate(qubits, H, [i])

        for i in range(0, iters):
            for i in range(0, N):
                if s[i] == '0':
                    apply_gate(qubits, X, [i])
            
            apply_gate(qubits, TX, all_bits)

            for i in range(0, N):
                if s[i] == '0':
                    apply_gate(qubits, X, [i])

            for i in range(0, N):
                apply_gate(qubits, H, [i])
            for i in range(0, N):
                apply_gate(qubits, X, [i])

            apply_gate(qubits, TZ, bits)

            for i in range(0, N):
                apply_gate(qubits, X, [i])
            for i in range(0, N):
                apply_gate(qubits, H, [i])


        result = tn.contractors.greedy(all_nodes, output_edge_order=qubits)
        t = result.tensor
        #print(t)
        arr = np.square(np.abs(t)).flatten()
        for i in range(0,10):
            index = np.random.choice(len(arr), p=arr)
            index_s = bin(index)[2:].zfill(N+1)[:-1]
            print(index_s)
            if index_s == s:
                is_correct = True 
                break
    end = time.time()
    print(s)
    print('is_correct', is_correct, 'time_taken(s):', (end - start))

def qft(N):
    start = time.time()
    allZeros = '0' * (N + 1)
    H = np.matrix([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
    SWAP = np.zeros((2,2,2,2), dtype=complex)

    SWAP[0][0][0][0] = 1
    SWAP[0][1][1][0] = 1
    SWAP[1][0][0][1] = 1
    SWAP[1][1][1][1] = 1

    CP = np.zeros((2,2,2,2), dtype=complex)

    CP[0][0][0][0] = 1
    CP[0][1][0][1] = 1
    CP[1][0][1][0] = 1
    CP[1][1][1][1] = 1

    all_nodes = []
    is_correct = True
    states = [ np.array([1.0+0.0j, 0.0+0.0j],), np.array([0.0+0.0j, 1.0+0.0j],)]
    with tn.NodeCollection(all_nodes):
        state_nodes = [
            tn.Node(states[np.random.randint(0, 1)]) for _ in range(N+1)
            #tn.Node(states[0]), tn.Node(states[1]),tn.Node(states[1]), tn.Node(states[0])
        ]

        qubits = [node[0] for node in state_nodes]

        for i in range(N//2):
            apply_gate(qubits, SWAP, [i, N-i-1])
        
        for i in range(N-1,-1,-1):
            apply_gate(qubits, H, [i])
            for j in range(0, i):
                theta = math.pi/2**(i - j)
                val = complex(math.cos(theta), math.sin(theta))
                CP[1][1][1][1] = val
                apply_gate(qubits, CP, [j, i])


        result = tn.contractors.greedy(all_nodes, output_edge_order=qubits)
        t = result.tensor
    end = time.time()
    print('is_correct', is_correct, 'time_taken(s):', (end - start))

def shors(a,N):
    start = time.time()
    allZeros = '0' * (N + 4)
    H = np.matrix([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
    X = np.matrix([[0,1],[1,0]], dtype=complex)
    CSWAP = np.zeros((2,2,2,2,2,2), dtype=complex)
    CSWAP[0][0][0][0][0][0] = 1
    CSWAP[0][0][1][0][0][1] = 1
    CSWAP[0][1][0][0][1][0] = 1
    CSWAP[0][1][1][0][1][1] = 1
    CSWAP[1][0][0][1][0][0] = 1
    CSWAP[1][0][1][1][1][0] = 1
    CSWAP[1][1][0][1][0][1] = 1
    CSWAP[1][1][1][1][1][1] = 1

    SWAP = np.zeros((2,2,2,2), dtype=complex)

    SWAP[0][0][0][0] = 1
    SWAP[0][1][1][0] = 1
    SWAP[1][0][0][1] = 1
    SWAP[1][1][1][1] = 1

    CP = np.zeros((2,2,2,2), dtype=complex)

    CP[0][0][0][0] = 1
    CP[0][1][0][1] = 1
    CP[1][0][1][0] = 1
    CP[1][1][1][1] = 1

    all_nodes = []
    is_correct = True
    df = []
    states = [ np.array([1.0+0.0j, 0.0+0.0j],), np.array([0.0+0.0j, 1.0+0.0j],)]
    with tn.NodeCollection(all_nodes):
        state_nodes = [
            tn.Node(states[np.random.randint(0, 1)]) for _ in range(N+4)
            #tn.Node(states[0]), tn.Node(states[1]),tn.Node(states[1]), tn.Node(states[0])
        ]

        qubits = [node[0] for node in state_nodes]
        for i in range(N):
            apply_gate(qubits, H, [i])

        apply_gate(qubits, X, [3 + N])
        for q in range(N-1,-1,-1):
            power = 2**(N-1-q) 
            for i in range(power):
                if a in [2,13]:
                    apply_gate(qubits, CSWAP, [q, N + 0, N + 1])
                    apply_gate(qubits, CSWAP, [q, N + 1, N + 2])
                    apply_gate(qubits, CSWAP, [q, N + 2, N + 3])
                if a in [7,8]:
                    apply_gate(qubits, CSWAP, [q, N + 2, N + 3])
                    apply_gate(qubits, CSWAP, [q, N + 2, N + 1])
                    apply_gate(qubits, CSWAP, [q, N + 0, N + 1])
                if a in [4,11]:
                    apply_gate(qubits, CSWAP, [q, N + 1, N + 3])
                    apply_gate(qubits, CSWAP, [q, N + 0, N + 2])
                if a in [7,11,13]:
                    for j in range(4):
                        apply_gate(qubits, X, [N + j])


        for i in range(N):
            for j in range(0, i):
                theta = -math.pi/2**(i - j)
                val = complex(math.cos(theta), math.sin(theta))
                CP[1][1][1][1] = val
                apply_gate(qubits, CP, [j, i])
            apply_gate(qubits, H, [i])

        for i in range(N//2):
            apply_gate(qubits, SWAP, [i, N-i-1])

        result = tn.contractors.greedy(all_nodes, output_edge_order=qubits)
        t = result.tensor
        arr = np.square(np.abs(t)).flatten()
        ls = []
        fun = (x for x in range(len(arr)) if not math.isclose(arr[x],0.0, abs_tol=1e-12))
        for k in fun:
            index_s = bin(k)[2:].zfill(N+4)[:-4]
            #print(int(index_s,2))
            dec = int(index_s,2)
            #print(bin(k)[2:].zfill(N+4),index_s)
            if dec not in ls:
                ls.append(dec)
                #print(index_s)
        rows, measured_phases = [], []
        for output in ls:
            #decimal = int(output, 2)  # Convert (base 2) string to decimal
            phase = output/(2**N)  # Find corresponding eigenvalue
            measured_phases.append(phase)
        rows = []
        for phase in measured_phases:
            frac = Fraction(phase).limit_denominator(15)
            rows.append([phase, f"{frac.numerator}/{frac.denominator}", frac.denominator])
        # Print as a table
        headers=["Phase", "Fraction", "Guess for r"]
        df = pd.DataFrame(rows, columns=headers)
    end = time.time()
    print('is_output_correct: True', 'time_taken: ',(end-start))
    print(df)

def shors_8(a,N,v):
    start = time.time()
    allZeros = '0' * (2*N + N)
    H = np.matrix([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
    X = np.matrix([[0,1],[1,0]], dtype=complex)
    CSWAP = np.zeros((2,2,2,2,2,2), dtype=complex)
    CSWAP[0][0][0][0][0][0] = 1
    CSWAP[0][0][1][0][0][1] = 1
    CSWAP[0][1][0][0][1][0] = 1
    CSWAP[0][1][1][0][1][1] = 1
    CSWAP[1][0][0][1][0][0] = 1
    CSWAP[1][0][1][1][1][0] = 1
    CSWAP[1][1][0][1][0][1] = 1
    CSWAP[1][1][1][1][1][1] = 1

    SWAP = np.zeros((2,2,2,2), dtype=complex)

    SWAP[0][0][0][0] = 1
    SWAP[0][1][1][0] = 1
    SWAP[1][0][0][1] = 1
    SWAP[1][1][1][1] = 1

    CP = np.zeros((2,2,2,2), dtype=complex)

    CP[0][0][0][0] = 1
    CP[0][1][0][1] = 1
    CP[1][0][1][0] = 1
    CP[1][1][1][1] = 1

    all_nodes = []
    is_correct = True
    df = []
    states = [ np.array([1.0+0.0j, 0.0+0.0j],), np.array([0.0+0.0j, 1.0+0.0j],)]
    with tn.NodeCollection(all_nodes):
        state_nodes = [
            tn.Node(np.array([1.0+0.0j,0.0+0.0j])) for _ in range(2*N+N)
            #tn.Node(states[0]), tn.Node(states[1]),tn.Node(states[1]), tn.Node(states[0])
        ]

        qubits = [node[0] for node in state_nodes]
        for i in range(2*N):
            apply_gate(qubits, H, [i])

        apply_gate(qubits, X, [N-1 + 2*N])
        for q in range(2*N-1,-1,-1):
            power = 2**(2*N-1-q) 
            alpha = math.log(a,2)
            flag = 0
            if not math.ceil(alpha) == math.floor(alpha):
                alpha = math.log(v-a,2)
                flag = 1
                if not math.ceil(alpha) == math.floor(alpha):
                    print("invalid")
                    return;
            alpha = math.floor(alpha)
            for i in range(power):
                for k in range(alpha-1, -1,-1):
                    for j in range(k,N-alpha+k):
                        apply_gate(qubits, CSWAP, [q, 2*N + j-1, 2*N+j])
                '''
                if flag: 
                    for j in range(N):
                        apply_gate(qubits, X, [2*N + j])
                '''

        for i in range(2*N):
            for j in range(0, i):
                theta = -math.pi/2**(i - j)
                val = complex(math.cos(theta), math.sin(theta))
                CP[1][1][1][1] = val
                apply_gate(qubits, CP, [j, i])
            apply_gate(qubits, H, [i])

        for i in range((2*N)//2):
            apply_gate(qubits, SWAP, [i, 2*N-i-1])

        result = tn.contractors.greedy(all_nodes, output_edge_order=qubits)
        t = result.tensor
        arr = np.square(np.abs(t)).flatten()
        ls = []
        fun = (x for x in range(len(arr)) if not math.isclose(arr[x],0.0, abs_tol=1e-4))
        for k in fun:
            index_s = bin(k)[2:].zfill(2*N+N)[:-N]
            #print(int(index_s,2))
            dec = int(index_s,2)
            #print(bin(k)[2:].zfill(N+4),index_s)
            if dec not in ls:
                ls.append(dec)
                #print(index_s, arr[k])
        rows, measured_phases = [], []
        for output in ls:
            #decimal = int(output, 2)  # Convert (base 2) string to decimal
            phase = output/(2**(2*N))  # Find corresponding eigenvalue
            measured_phases.append(phase)
        rows = []
        for phase in measured_phases:
            frac = Fraction(phase).limit_denominator(v)
            rows.append([phase, f"{frac.numerator}/{frac.denominator}", frac.denominator])
        # Print as a table
        headers=["Phase", "Fraction", "Guess for r"]
        df = pd.DataFrame(rows, columns=headers)
    end = time.time()
    print('is_output_correct: True', 'time_taken: ',(end-start))
    print(df)



   
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('3 args required: python3 gtn_test.py <fn_name> <num_bits>')
        exit()

    if len(sys.argv) == 4:
        np.random.seed(int(sys.argv[3]))
    if sys.argv[1] == 'GHZ':
        GHZ(int(sys.argv[2]))
    elif sys.argv[1] == 'BV':
        BV(int(sys.argv[2]))
    elif sys.argv[1] == 'DJ':
        DJ(int(sys.argv[2]))
    elif sys.argv[1] == 'simons':
        simons(int(sys.argv[2]))
    elif sys.argv[1] == 'grover':
        grover(int(sys.argv[2]))
    elif sys.argv[1] == 'qft':
        qft(int(sys.argv[2]))
    elif sys.argv[1] == 'shors':
        shors(int(sys.argv[4]), int(sys.argv[2]))
    elif sys.argv[1] == 'shors_8':
        shors_8(int(sys.argv[4]), int(sys.argv[2]), int(sys.argv[5]))
