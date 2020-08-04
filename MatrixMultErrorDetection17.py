#!/usr/bin/env python
'''
Polynomial code with fast decoding
'''
import csv
from mpi4py import MPI
import numpy as np
import random
import threading
import time
from fft import polycode_ifft
import math
import statistics

# Change to True for more accurate timing, sacrificing performance
barrier = True
# Change to True to imitate straggler effects
straggling = False

def loop():
  t = time.time()
  while time.time() < t + 60:
    a = 1 + 1

comm = MPI.COMM_WORLD
parameters=[(16,2,2), (16,2,4),(16,3,3),(16,3,4),(16,5,3), (16, 15,1)]
dimensions=[(2400,3200,3200)]
EncodingData={}
WaitData={}
SendData={}
DecData={}

for N,m,n in parameters:
    for s,r,t in dimensions:

        ##################### Parameters ########################
        # Use one master and N workers
        length=N
        # Field size assumed to be prime for this implementation
        F =65537

        # Input matrix size - A: s by r, B: s by t
        if(r%m!=0):
            r+=m-r%m
        if(t%n!=0):
            t+=n-t%n
        #s = 2400
        #r = 3200
        #t = 3200
        # Pick a primitive root 64
        prim_root=3
        rt= pow(prim_root,int(((F-1)/(N))), F)

        var = [pow(rt, i, F) for i in range(length)]
        #+ [pow(prim_root, i, F) for i in range(1, (int(N-length))+1) if (i%int(((F-1)/(m*n)))!=0)]

        #########################################################
        encTimes=[]
        sendTimes=[]
        waitTimes=[]
        decTimes=[]
        for repeat in range(20):
            if comm.rank == 0:
              # Master
              #print ("Running with %d processes:" % comm.Get_size())

              # Decide and broadcast chose straggler
              #straggler = random.randint(1, N)
              #for i in range(N):
                #comm.send(straggler, dest=i+1, tag=7)

              # Create random matrices of 8-bit ints

              A = np.matrix(np.random.randint(0,2**8-1 , (r, s)) )
              B = np.matrix(np.random.randint(0,2**8-1, (t, s)) )

              # Split the matrices
              Ap = np.split(A, m)
              Bp = np.split(B, n)

              # Encode the matrices
              encStart=time.time()
              Aenc= [sum([Ap[j] * (pow(var[i], j, F)) for j in range(m)]) % F for i in range(N)]
              Benc= [sum([Bp[j] * (pow(var[i], j * m, F)) for j in range(n)]) % F for i in range(N)]
              encEnd=time.time()
              encodingTime=encEnd-encStart
              #print(encEnd-encStart) #encoding time
              # Initialize return dictionary
              encTimes.append(encodingTime)
              Rdict = []
              for i in range(N):
                Rdict.append(np.zeros((int(r/m), int(t/n)), dtype=np.int_))

              # Start requests to send and receive
              reqA = [None] * N
              reqB = [None] * N
              reqC = [None] * N

              bp_start = time.time()

              for i in range(N):
                reqA[i] = comm.Isend([Aenc[i], MPI.INT], dest=i+1, tag=15)
                reqB[i] = comm.Isend([Benc[i], MPI.INT], dest=i+1, tag=29)
                reqC[i] = comm.Irecv([Rdict[i], MPI.INT], source=i+1, tag=42)

              MPI.Request.Waitall(reqA)
              MPI.Request.Waitall(reqB)

              # Optionally wait for all workers to receive their submatrices, for more accurate timing
              if barrier:
                comm.Barrier()

              bp_sent = time.time()
              sendTimes.append((bp_sent - bp_start))

              Crtn = [None] * N
              lst = []
              # Wait for the mn fastest workers
              for i in range(N):
                j = MPI.Request.Waitany(reqC)
                lst.append(j)
                Crtn[j] = Rdict[j]
              bp_received = time.time()

              waitTimes.append(bp_received - bp_sent)

              Crtn=polycode_ifft(Crtn, F, prim_root)
              bp_done = time.time()

              #print ((bp_done - bp_received))
              decTimes.append(bp_done - bp_received)
              format_str='{:0'+str(int((math.log(length,2))))+'b}'
               #  time spent decoding
              parity_check=[not Crtn[int(format_str.format(i)[::-1], 2)].any() for i in range(m*n, N)]
              if all(parity_check):
                  print("no errors")
              else:
                  print("errors")
              #Copy=Crtn.copy()

              #for i in range(length): #bit reverse

                  #bit_rev=int(format_str.format(i)[::-1], 2)
                  #Crtn[i]=Copy[bit_rev]

              #bit_reverse = [0, 2, 1, 3]
              #Cver = [(Ap[bit_reverse[int(i / 4)]] * Bp[bit_reverse[i % 4]].getT()) % F for i in range(m * n)]
              #Cver = [(Ap[int(i %m)] * Bp[int(i / m)].getT()) % F for i in range(m * n)]

              #print(Cver)
              #print(Crtn)
              #print ([np.array_equal(Crtn[int(format_str.format(i)[::-1], 2)], Cver[i]) for i in range(m * n)])
            else:
              # Worker
              # Receive straggler information from the master
              #straggler = comm.recv(source=0, tag=7)

              # Receive split input matrices from the master
              Ai = np.empty_like(np.matrix([[0]*s for i in range(int(r/m))]))
              Bi = np.empty_like(np.matrix([[0]*s for i in range(int(t/n))]))
              rA = comm.Irecv(Ai, source=0, tag=15)
              rB = comm.Irecv(Bi, source=0, tag=29)

              rA.wait()
              rB.wait()

              if barrier:
                comm.Barrier()
              wbp_received = time.time()

              # Start a separate thread to mimic background computation tasks if this is a straggler

              Ci = (Ai * (Bi.getT())) % F
              wbp_done = time.time()
              # print "Worker %d computing takes: %f\n" % (comm.Get_rank(), wbp_done - wbp_received)

              sC = comm.Isend(Ci, dest=0, tag=42)
              sC.Wait()
        if(comm.rank==0):
            avgEnc=statistics.mean(encTimes)
            stdEnc=statistics.stdev(encTimes)
            avgSend=statistics.mean(sendTimes)
            stdSend=statistics.stdev(sendTimes)
            avgWait=statistics.mean(waitTimes)
            stdWait=statistics.stdev(waitTimes)
            avgDec=statistics.mean(decTimes)
            stdDec=statistics.stdev(decTimes)
    if(comm.rank==0):
        EncodingData[(N,m,n,s,r,t)]=(avgEnc, stdEnc)
        SendData[(N,m,n,s,r,t)]=(avgSend, stdSend)
        WaitData[(N,m,n,s,r,t)]=(avgWait, stdWait)
        DecData[(N,m,n,s,r,t)]=(avgDec, stdDec)
if comm.rank==0:
    w = csv.writer(open("Poly_Times.csv", "w"))
    w.writerow(["parameters", (""), "EncTimeMean", "EncTimeStdev", "", "WaitTimeMean", "WaitTimestdev",(""),"SendTimeMean", "SendTimestdev",(""),
    "DecTimeMean", "DecTimeStdev" ] )
    for key, val in EncodingData.items():
        w.writerow([key, (""),  val[0], val[1], (""), WaitData[key][0], WaitData[key][1],(""),SendData[key][0], SendData[key][1],(""),
        DecData[key][0], DecData[key][1] ] )
