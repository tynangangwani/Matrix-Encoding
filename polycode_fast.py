#!/usr/bin/env python
'''
Polynomial code with fast decoding
'''

from mpi4py import MPI
import numpy as np
import random
import threading
import time
import math
from fft import polycode_ifft

# Change to True for more accurate timing, sacrificing performance
barrier = True
# Change to True to imitate straggler effects
straggling = False

def loop():
  t = time.time()
  while time.time() < t + 60:
    a = 1 + 1

##################### Parameters ########################
# Use one master and N workers
N = 17

# Matrix division
m = 4
n = 4

# Field size assumed to be prime for this implementation
F = 65537

# Input matrix size - A: s by r, B: s by t
s = 4000
r = 4000
t = 4000

# Pick a primitive root 64
rt = rt = pow(3,int(((F-1)/(m*n))), F)

# Values of x_i used by 17 workers
var = [pow(rt, i, 65537) for i in range(16)] + [3]
#########################################################

A = np.matrix(np.random.random_integers(0, 255, (r, s)))
B = np.matrix(np.random.random_integers(0, 255, (t, s)))

# Split the matrices
Ap = np.split(A, m)
Bp = np.split(B, n)

# Encode the matrices
encStart=time.time()
Aenc= [sum([Ap[j] * (pow(var[i], j, F)) for j in range(m)]) % F for i in range(N)]
Benc= [sum([Bp[j] * (pow(var[i], j * m, F)) for j in range(n)]) % F for i in range(N)]
encEnd=time.time()
print("Encoding Time: "+ str(encEnd-encStart))

# for i in range(N):
lst=list(range(1,17))#node 0 is a straggler
print(lst)
CompStart=time.time()
Crtn=[(Aenc[i] * (Benc[i].getT())) % F for i in range(N)]
CompEnd=time.time()
print("Computation time: "+ str(CompEnd-CompStart))
Crtn[0]=0 #node 0 is a straggler
DecStart=time.time()
missing = set(range(m * n)) - set(lst)

# Fast decoding hard coded for m, n = 4
sig = 4
xlist = [var[i] for i in lst]

for i in missing:
    begin = time.time()
    coeff = [1] * (m * n) #coeff size of m*n all ones
    for j in range(m * n):
  # Compute coefficient
        for k in set(lst) - set([lst[j]]): #for k in the list returned except j (current one)
          coeff[j] = (coeff[j] * (var[i] - var[k]) * pow(var[lst[j]] - var[k], F - 2, F)) % F
      #coeff at j = (x_i(missing one) -x_k)*
    Crtn[i] = sum([Crtn[lst[j]] * coeff[j] for j in range(16)]) % F

Crtn=polycode_ifft(Crtn[0:16], F, 3)
DecEnd=time.time()
print("decode time: "+ str(DecEnd-DecStart))
# Verify correctness
# Bit reverse the order to match the FFT
# To obtain outputs in an ordinary order, bit reverse the order of input matrices prior to FFT
bit_reverse = [0, 2, 1, 3]
Cver = [(Ap[bit_reverse[int(i / 4)]] * Bp[bit_reverse[i % 4]].getT()) % F for i in range(m * n)]
#Cver = [(Ap[int(i %4)] * Bp[int(i / 4)].getT()) % F for i in range(m * n)]
print(Cver)
print(Crtn)
print ([np.array_equal(Crtn[i], Cver[i]) for i in range(m * n)])
