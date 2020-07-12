#!/usr/bin/env python
'''
Polynomial code with fast decoding
'''

from mpi4py import MPI
import numpy as np
import random
import threading
import time
from fft import polycode_ifft
import math

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
N = 35

# Matrix division
m = 4
n = 8

length=m*n
# Field size assumed to be prime for this implementation
F =23725313

# Input matrix size - A: s by r, B: s by t
s = 16
r = 16
t = 16

# Pick a primitive root 64
prim_root=3
rt= pow(prim_root,int(((F-1)/(m*n))), F)

var = [pow(rt, i, F) for i in range(length)] + [pow(prim_root, i, F) for i in range(1, (int(N-length))+1) if (i%int(((F-1)/(m*n)))!=0)]

#########################################################

comm = MPI.COMM_WORLD

if comm.rank == 0:
  # Master
  print ("Running with %d processes:" % comm.Get_size())

  # Decide and broadcast chose straggler
  straggler = random.randint(1, N)
  for i in range(N):
    comm.send(straggler, dest=i+1, tag=7)

  # Create random matrices of 8-bit ints

  A = np.matrix(np.random.randint(0,2**12-1 , (r, s)) )
  B = np.matrix(np.random.randint(0,2**12-1, (t, s)) )

  # Split the matrices
  Ap = np.split(A, m)
  Bp = np.split(B, n)

  # Encode the matrices
  encStart=time.time()
  Aenc= [sum([Ap[j] * (pow(var[i], j, F)) for j in range(m)]) % F for i in range(N)]
  Benc= [sum([Bp[j] * (pow(var[i], j * m, F)) for j in range(n)]) % F for i in range(N)]
  # Initialize return dictionary
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
  print ("Time spent sending all messages is: %f" % (bp_sent - bp_start))

  Crtn = [None] * N
  lst = []
  # Wait for the mn fastest workers
  for i in range(m * n):
    j = MPI.Request.Waitany(reqC)
    lst.append(j)
    Crtn[j] = Rdict[j]
  bp_received = time.time()
  print ("Time spent waiting for %d workers %s is: %f" % (m * n, ",".join(map(str, [x + 1 for x in lst])), (bp_received - bp_sent)))
  print(lst)
  missing = set(range(m * n)) - set(lst)
  print("Missing: "+str(missing))

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
      print(i)
      print(Crtn[i])


  Crtn=polycode_ifft(Crtn[0:length], F, prim_root)
  bp_done = time.time()
  print ("Time spent decoding is: %f" % (bp_done - bp_received))

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
  print ([np.array_equal(Crtn[i], Cver[i]) for i in range(m * n)])
else:
  # Worker
  # Receive straggler information from the master
  straggler = comm.recv(source=0, tag=7)

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
  if straggling:
    if straggler == comm.rank:
      t = threading.Thread(target=loop)
      t.start()

  Ci = (Ai * (Bi.getT())) % F
  wbp_done = time.time()
  # print "Worker %d computing takes: %f\n" % (comm.Get_rank(), wbp_done - wbp_received)

  sC = comm.Isend(Ci, dest=0, tag=42)
  sC.Wait()
