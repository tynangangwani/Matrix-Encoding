#!/usr/bin/env python
'''
Polynomial code with fast decoding
'''

import statistics
import numpy as np
import random
import threading
import time
import math
import csv
from fft import polycode_ifft
import random

def loop():
  t = time.time()
  while time.time() < t + 60:
    a = 1 + 1
EncodingTimes={}
ComputationTimes={}
DecodingTimes={}
##################### Parameters ########################
# Use one master and N workers
#N = 66

# Matrix division
#m = 8
#n = 8
#list in the form (Divisions of A, Divisions of B, total number of workers)
#requirement m*n=
parameters=[(2,2, 5), (2,2, 6), (2,2, 7), (2,4, 9), (2,4,11), (2,4,13), (2,4,15), (4,4,17),
            (4,4,21), (4,4,25), (4,4,29), (4,8,33),(4,8,41), (4,8,49), (4,8,57), (8,8, 65), (8,8, 65),
            (8,8, 81), (8,8, 97), (8,8, 113)]

# Field size assumed to be prime for this implementation
F = 65537

# Input matrix size - A: s by r, B: s by t
s = 128
r = 128
t = 128

# Pick a primitive root 64

for (m,n,N) in parameters:
    EncSamples=[]
    CompSamples=[]
    DecSamples=[]
    r*=m
    t*=n
    for sample in range(10):

        length=m*n
        rt = rt = pow(3,int(((F-1)/(m*n))), F)
        # Values of x_i used by 17 workers
        var = [pow(rt, i, 65537) for i in range(length)] + [pow(3, i, F) for i in range(1, (int(N-length))+1)]
        #########################################################

        A = np.matrix(np.random.randint(0, 255, (r, s)))
        B = np.matrix(np.random.randint(0, 255, (t, s)))

        # Split the matrices
        Ap = np.split(A, m)
        Bp = np.split(B, n)

        # Encode the matrices
        encStart=time.time()
        Aenc= [sum([Ap[j] * (pow(var[i], j, F)) for j in range(m)]) % F for i in range(N)]
        Benc= [sum([Bp[j] * (pow(var[i], j * m, F)) for j in range(n)]) % F for i in range(N)]
        encEnd=time.time()
        EncSamples.append(encEnd-encStart)

        # for i in range(N):
        lst=list(range(N))#node N-length nodes are stragglers

        CompStart=time.time()
        Crtn=[(Aenc[i] * (Benc[i].getT())) % F for i in range(N)]
        CompEnd=time.time()
        CompSamples.append(CompEnd-CompStart)

        stragglers=[]
        for i in range(N-length):
            straggler=random.randint(0,len(lst)-1)
            stragglers.append(lst[straggler])
            Crtn[lst[straggler]]=[None]
            del lst[straggler]

        #print(lst)
        DecStart=time.time()
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


        Crtn=polycode_ifft(Crtn[0:length], F, 3)
        DecEnd=time.time()
        DecSamples.append(DecEnd-DecStart)

        # Verify correctness
        # Bit reverse the order to match the FFT
        # To obtain outputs in an ordinary order, bit reverse the order of input matrices prior to FFT
        Copy=Crtn.copy()
        for i in range(length): #bit reverse
            format_str='{:0'+str(int((math.log(length,2))))+'b}'
            bit_rev=int(format_str.format(i)[::-1], 2)
            Crtn[i]=Copy[bit_rev]
        #bit_reverse = [0, 2, 1, 3]
        #Cver = [(Ap[bit_reverse[int(i / 4)]] * Bp[bit_reverse[i % 4]].getT()) % F for i in range(m * n)]
        Cver = [(Ap[int(i %m)] * Bp[int(i / m)].getT()) % F for i in range(m * n)]
        #print(Cver)
        #print(Crtn)
        for i in range(m * n):
            if not np.array_equal(Crtn[i], Cver[i]):
                print("FALSE RESULTS FOR: ")
                print((m,n,N))
    #print(EncSamples)
    EncodingTimeMean=statistics.mean(EncSamples)
    EncodingTimeStdev=statistics.stdev(EncSamples)
    EncodingTimes[(m,n,N)]=(EncodingTimeMean, EncodingTimeStdev)
    #print(str((EncodingTimeMean, EncodingTimeStdev))

    #print(CompSamples)
    CompTimeMean=statistics.mean(CompSamples)
    CompTimeStdev=statistics.stdev(CompSamples)
    ComputationTimes[(m,n,N)]=(CompTimeMean, CompTimeStdev)
    #print((CompTimeMean, CompTimeStdev))

    #print(DecSamples)
    DecTimeMean=statistics.mean(DecSamples)
    DecTimeStdev=statistics.stdev(DecSamples)
    DecodingTimes[(m,n,N)]=(DecTimeMean, DecTimeStdev)
    #print((CompTimeMean, CompTimeStdev))
    print("done "+ str((m,n,N)) )
w = csv.writer(open("Poly_Times512.csv", "w"))
w.writerow(["parameters", "EncTimeMean", "EncTimeStdev", "", "CompTimeMean", "CompTimestdev",(""),
"DecTimeMean", "DecTimeStdev" ] )
for key, val in EncodingTimes.items():
    w.writerow([key, val[0], val[1], (""), ComputationTimes[key][0], ComputationTimes[key][1],(""),
    DecodingTimes[key][0], DecodingTimes[key][1] ] )
