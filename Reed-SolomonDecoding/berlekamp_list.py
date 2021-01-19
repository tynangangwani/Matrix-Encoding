import numpy as np
import random
import threading
import time
from fft import polycode_ifft
from fft import fft
import math
from DecodingRSVectorized import decodeRS
from DecodingRSVectorized import normalize
import copy

def eval(n, poly, F):
    sum=0
    for i in range(len(poly)):
        sum+= (pow(n,i,F))*poly[i]%F
        sum=sum%F
    return sum

#(k,n) code, F is the field and prim root is the primitve root of the field
t=[2, 9, 1 ,4,0, 0,0,0, ]#0,0,0,0,0,0,0,0]
length=len(t)
msgLen=4
F=65537
prim_root=3
alpha=pow(prim_root, int((F-1)/length)) #wth root of unity




print("len of t"+ str(len(t)))
t=fft(t,65537,3)

F=65537
print("start")
print(t)
t[5]=t[5]+5%65537
t[6]=t[6]+234%65537
t[7]=t[7]+1234%65537

t[2]=t[2]+21234%65537

t[1]=t[1]+2314123%65537

#t[9]=t[9]+45645%65537
#print(polycode_ifft(t, 65537, 3))

t=fft(t, 65537, 3)
print(t)
L = 0
m = 1
b = 1
B=[1]
N= length-msgLen#number of syndromes
print(N)
C=t[1:N+1]
X=[1]
d=0
print("c: ", C)
  #k=
  #/* steps 2. and 6. */
  #for (n = 0; n < N; n++) {
for n in range(N):
      #/* step 2. calculate discrepancy */


    temp=[C[n]] + [(C[(n-i)%N]*X[i])%F  for i in range(1, L+1)]


    d = sum(temp)%F #\Sigma_{i=1}^L c_i * s_{n-i};


    if(d == 0):
          #/* step 3. discrepancy is zero; annihilation continues */
          m = m + 1


    elif ((2 * L) <= n):

          #  /* step 5. */
          #  /* temporary copy of C(x) */

          #  polynomial(field K) T(x) = C(x);
          T=copy.deepcopy(X)

          degree= max(len(X), m+len(B))-1
          temp = [0]*degree
         # print(B)
          B=[0]*m + B
          #print(B)
          if len(X)>len(B):
              B=B+[0]*(len(X)-len(B))
          else:
              X=X+[0]*(len(B)-len(X))

          X=[(X[i]- (d * pow(b, F-2, F)* B[i]) )%F for i in range(degree+1)]

          normalize(X)
          normalize(B)

         # C(x) - d b^{-1} x^m B(x);
          L = n + 1 - L
          B = T

          b = d
          m = 1
    else:

         degree= max(len(X), m+len(B))
         temp = [0]*degree
         Btemp=B
         B=[0]*m + B
         if len(X)>len(B):
             B=B+[0]*(len(X)-len(B))
         else:
             X=X+[0]*(len(B)-len(X))

         X=[(X[i]- ((d * pow(b, F-2, F) * B[i])%F) ) %F  for i in range(degree)]
         normalize(X)
         normalize(B)


         #C(x) = C(x) - d b^{-1} x^m B(x);
         B=Btemp
         m = m + 1
print("d:")
print(d)
print("error locator berlekamp")
print(X)
print("L:" + str(L))
for i in range(length):
    print("Evaluated at "+str(i)+"^-1")
    print(eval(pow(alpha, length-i, F),X, F))
#  return L
