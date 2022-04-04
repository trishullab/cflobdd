import random
import quimb as qu
import quimb.tensor as qtn
import time
import sys
import random
import numpy as np

def GHZ(N):
	allOnes = '1' * N
	allZeros = '0' * N
	start = time.time()
	circ = qtn.Circuit(N)

	# randomly permute the order of qubits
	regs = list(range(N))
	random.shuffle(regs)

	# hamadard on one of the qubits
	circ.apply_gate('H', regs[0])

	# chain of cnots to generate GHZ-state
	for i in range(N - 1):
	    circ.apply_gate('CNOT', regs[i], regs[i + 1])

	mid = time.time()
	# sample it once
	check_ones = False
	check_zeros = False
	correct_out = False
	for b in circ.sample(10):
		if b == allOnes:
			check_ones = True
		elif b == allZeros:
			check_zeros = True
		else:
			correct_out = False
			break
		if check_ones and check_zeros:
			correct_out = True
			break

	end = time.time()

	print('is_output_correct: ', correct_out, ' time_taken(s): ', (end - start))

def qft(N):
	start = time.time()
	circ = qtn.Circuit(N)

	for i in range(N-1, -1, -1):
		circ.apply_gate('H', i)
		for j in range(i):
			circ.apply_gate('CU1', np.pi/2**(i - j), j, i)

	for i in range(N//2):
		circ.apply_gate('SWAP', i, N-i-1)

	for b in circ.sample(10):
		print(b)
	end = time.time()
	print('is_output_correct:', True, 'time_taken(s):', (end-start))


def BV(N):

	s = ""
	for i in range(0, N):
		r = random.randint(0, 1)
		if r == 0:
			s = s + '0'
		else:
			s = s + '1'

	start = time.time()
	circ = qtn.Circuit(N+1)

	circ.apply_gate('X', N)
	for i in range(0, N+1):
		circ.apply_gate('H', i)

	for i in range(0, N):
		if s[i] == '1':
			circ.apply_gate('CNOT', i, N)

	for i in range(0, N):
		circ.apply_gate('H', i)

	sampled_s = ""
	for b in circ.sample(1):
		sampled_s = b[:-1]
		if sampled_s == s:
			break

	end = time.time()
	print('s: ', s, 'sampled_s: ', sampled_s)
	print('is_output_correct:', (sampled_s == s), 'time_taken(s):', (end-start))


def grover(N):

	s = ""
	for i in range(0, N):
		r = random.randint(0, 1)
		if r == 0:
			s = s + '0'
		else:
			s = s + '1'

	start = time.time()
	circ = qtn.Circuit(N)

	for i in range(0, N):
		circ.apply_gate('H', i)
	
	
def DJ(N):
	s = ""
	for i in range(0, N):
		r = random.randint(0, 1)
		if r == 0:
			s = s + '0'
		else:
			s = s + '1'
	allZeros = '0' * N
	is_output_correct = True
	start = time.time()
	circ = qtn.Circuit(N+1)

	circ.apply_gate('X', N)
	for i in range(N+1):
		circ.apply_gate('H', i)

	for i in range(N):
		if s[i] == '1':
			circ.apply_gate('X', i)	

	for i in range(N):
		circ.apply_gate('CNOT', i, N)

	for i in range(N):
		if s[i] == '1':
			circ.apply_gate('X', i)		

	for i in range(N):
		circ.apply_gate('H', i)

	sampled_s = ""
	for b in circ.sample(10):
		sampled_s = b[:-1]
		if sampled_s == allZeros:
			is_output_correct = False
			break

	end = time.time()
	print('s: ', s, 'sampled_s: ', sampled_s)
	print('is_output_correct:', is_output_correct, 'time_taken(s):', (end-start))

def simons(N):
	s = ""
	is_output_correct = True
	for i in range(0, N):
		r = random.randint(0, 1)
		if r == 0:
			s = s + '0'
		else:
			s = s + '1'
	
	allZeros = '0' * N
	equations = []

	start = time.time()
	circ = qtn.Circuit(2*N)

	for i in range(0, N):
		circ.apply_gate('H', i)

	for i in range(N):
		circ.apply_gate('CNOT', i, i+N)

	k = 0
	for i in range(N-1, -1, -1):
		if s[i] == '1':
			m = N
			for j in range(N-1, -1, -1):
				if s[j] == '1':
					circ.apply_gate('CNOT', k, m)
				m += 1
			break
		k += 1

	for i in range(0, N):
		circ.apply_gate('H', i)

	for b in circ.sample(2*N):
		index_s = b[0:N]
		equations.append(index_s)
			
	end = time.time()
	print('s: ', s)
	# print(equations)
	print('is_output_correct:', is_output_correct, 'time_taken(s):', (end-start))

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('3 args required: python3 test.py <fn_name> <num_bits>')
		exit()
	
	if sys.argv[1] == 'GHZ':
		GHZ(int(sys.argv[2]))
	elif sys.argv[1] == 'BV':
		BV(int(sys.argv[2]))
	elif sys.argv[1] == 'grover':
		grover(int(sys.argv[2]))
	elif sys.argv[1] == 'DJ':
		DJ(int(sys.argv[2]))
	elif sys.argv[1] == 'simons':
		simons(int(sys.argv[2]))
	elif sys.argv[1] == 'qft':
		qft(int(sys.argv[2]))

