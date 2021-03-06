#!/usr/bin/env python
'''
Polynomial code with fast decoding
'''


import numpy as np
import random
import threading
import time
import math
import csv
from fft import polycode_ifft
from fft import fft
import random

def loop():
  t = time.time()
  while time.time() < t + 60:
    a = 1 + 1
EncodingTimes=[]
ComputationTimes=[]
DecodingTimes=[]
##################### Parameters ########################
# Use one master and N workers
N =32

# Matrix division
m = 4
n = 4

length=m*n
# Field size assumed to be prime for this implementation
F =1073738753#2147483137#4294957057#2147483137

# Input matrix size - A: s by r, B: s by t
s = 20
r =20
t = 20

# Pick a primitive root 64
prim_root=10
rt = rt = pow(prim_root,int(((F-1)/(N))), F)

def polyeval_A (Ap, var, i, m):
    #Aenc= [sum([Ap[j] * (pow(var[i], j, F)) for j in range(m)]) % F for i in range(N)]
    A=np.zeros_like(Ap[0])
    for j in range(m):
        A+=(( Ap[j] * (pow(var[i], j, F) ))%F)
        A=A%F
    #Aenc = (A %F)
    return A
def polyeval_B (Bp, var, i, n, m):
    #Aenc= [sum([Ap[j] * (pow(var[i], j, F)) for j in range(m)]) % F for i in range(N)]
    B=np.zeros_like(Bp[0])
    for j in range(n):
        B+= (Bp[j] * (pow(var[i], j*m, F) )) %F
        B=B%F
    #Benc=(B%F)
    return B


# Values of x_i used by 17 workers
var = [pow(rt, i, F) for i in range(N)]
#+ [pow(3, i, F) for i in range(1, (int(N-length))+1)]
#########################################################

A = np.matrix(np.random.randint(0, 2**15-1, (r, s)))
B = np.matrix(np.random.randint(0, 2**15-1, (t, s)))
#print(A)
#print(B)

# Split the matrices
Ap = np.split(A, m)
Bp = np.split(B, n)
print(Ap[0].dtype)
#print(A)

# Encode the matrices
encStart=time.time()
#Aenc= [sum([Ap[j] * (pow(var[i], j, F)) for j in range(m)]) % F for i in range(N)]#
Aenc=[polyeval_A(Ap, var, i, m) % F for i in range(N)]
#for i in range(N-m):
#    Ap.append(np.zeros_like(Ap[0]))
#print(Ap)
AP=Ap.copy()
#Aenc=fft(Ap,F,prim_root )

print("length of Aenc: "+str(len(Aenc)))

print(Aenc[0].dtype)

print("Aenc:")
#print(Aenc)
#a=polycode_ifft(Aenc, F, prim_root)
#print(a)
Benc= [polyeval_B(Bp, var, i, n, m) % F for i in range(N)]
#Benc= [sum([Bp[j] * (pow(var[i], j * m, F)) for j in range(n)]) % F for i in range(N)]
print("Benc:")
#print(Benc)
encEnd=time.time()
print("Encoding Time: "+ str(encEnd-encStart))

# for i in range(N):
lst=list(range(N))#node N-length nodes are stragglers

CompStart=time.time()
Crtn=[(Aenc[i] * (Benc[i].getT())) % F for i in range(N)]
#print (Crtn)


CompEnd=time.time()
Copy=Crtn.copy()
#no stragglers

'''
print("Computation time: "+ str(CompEnd-CompStart)) #node 0 is a straggler
stragglers=[]
for i in range(N-length):
    straggler=random.randint(0,len(lst)-1)
    stragglers.append(lst[straggler])
    Crtn[lst[straggler]]=[None]
    del lst[straggler]
print("Stragglers: "+ str(stragglers))
print(lst)

missing = set(range(m * n)) - set(lst)

# Fast decoding hard coded for m, n = 4
sig = 4
xlist = [var[i] for i in lst]

for i in missing:
  begin = time.time()
  coeff = [1] * (m * n)
  for j in range(m * n):
    # Compute coefficient
    for k in set(lst) - set([lst[j]]):

        coeff[j] = (coeff[j] * (var[i] - var[k]) * pow(var[lst[j]] - var[k], F - 2, F)) % F

  Crtn[i] = sum([ (Crtn[lst[j]] * coeff[j]) for j in range(length)]) % F
'''
DecStart=time.time()

Crtn=polycode_ifft(Crtn, F,prim_root)
#print(Crtn)

DecEnd=time.time()
print("decode time: "+ str(DecEnd-DecStart))
# Verify correctness
# Bit reverse the order to match the FFT
# To obtain outputs in an ordinary order, bit reverse the order of input matrices prior to FFT
Copy=Crtn.copy()
for i in range(N): #bit reverse
    format_str='{:0'+str(int((math.log(N,2))))+'b}'
    bit_rev=int(format_str.format(i)[::-1], 2)
    Crtn[i]=Copy[bit_rev]

#bit_reverse = [0, 2, 1, 3]
#Cver = [(Ap[bit_reverse[int(i / 4)]] * Bp[bit_reverse[i % 4]].getT()) % F for i in range(m * n)]
Cver = [(Ap[int(i %m)] * Bp[int(i / m)].getT()) % F for i in range(m * n)]
#Cver=[(AP[int(i)] * B.getT()) % F for i in range(m)]
#print(Cver)
#print(Crtn)
#print(Cver)
#print(Crtn)
print ([np.array_equal(Crtn[i], Cver[i]) for i in range(m*n)])
