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

def power_mod(a,p,F):
    return (pow(int(a),int(p),int(F)))
matrixPow=np.vectorize(power_mod)

def eval(n, poly, F):
    sum=0
    for i in range(len(poly)):
        sum+= (pow(n,i,F))*poly[i]%F
        sum=sum%F
    return sum
dimensions=100
#(k,n) code, F is the field and prim root is the primitve root of the field

length=16
msgLen=6
t=([np.random.randint(0,65536, (dimensions,dimensions)) for i in range(msgLen)]+[np.zeros((dimensions,dimensions), dtype=int)]*(length-msgLen))
F=65537
prim_root=3
alpha=pow(prim_root, int((F-1)/length)) #wth root of unity




print("len of t"+ str(len(t)))
t=fft(t,65537,3)

F=65537
print("start")
print(t)
t[5]=t[5]+np.random.randint(0,65536, (dimensions,dimensions))
t[6]=t[6]+np.random.randint(0,65536, (dimensions,dimensions))
t[7]=t[7]+np.random.randint(0,65536, (dimensions,dimensions))

t[2]=t[2]+np.random.randint(0,65536, (dimensions,dimensions))

t[1]=t[1]+np.random.randint(0,65536, (dimensions,dimensions))
#t[9]=t[9]+45645%65537
#print(polycode_ifft(t, 65537, 3))

#t=fft(t, 65537, 3)
print(t)

def berlekamp(length, msgLen, code, F, prim_root):

    dim=np.shape(code[0])
    print(dim)
    L = 0
    m = 1
    b = 1
    N= length-msgLen#number of syndromes
    B=[np.ones(dim, dtype=int)]
    t=fft(code, F, prim_root)
    #print(N)
    C=t[1:N+1]
    X=[np.ones(dim, dtype=int)]
    d=np.zeros(dim, dtype=int)
    #alpha=pow(prim_root, int((F-1)/length)) #wth root of unity

    #print("c: ", C)

    for n in range(N):
          #/* step 2. calculate discrepancy */


        temp=[C[n]] + [np.multiply(C[(n-i)%N],X[i])%F  for i in range(1, L+1)]


        d = sum(temp)%F #\Sigma_{i=1}^L c_i * s_{n-i};


        if(not d.any()):
              #/* step 3. discrepancy is zero; annihilation continues */
              m = m + 1


        elif ((2 * L) <= n):

              #  /* step 5. */
              #  /* temporary copy of C(x) */

              #  polynomial(field K) T(x) = C(x);
              T=copy.deepcopy(X)

              degree= max(len(X), m+len(B))-1
              #temp = [0]*degree
             # print(B)
              B=[np.zeros(dim, dtype=int)]*m + B
              #print(B)
              if len(X)>len(B):
                  B=B+[np.zeros(dim, dtype=int)]*(len(X)-len(B))
              else:
                  X=X+[np.zeros(dim, dtype=int)]*(len(B)-len(X))

              X=[(X[i]- (np.multiply(matrixPow(b, F-2, F), np.multiply(d, B[i]) ) ) )%F for i in range(len(X))]

              normalize(X)
              normalize(B)

             # C(x) - d b^{-1} x^m B(x);
              L = n + 1 - L
              B = T

              b = d
              m = 1
        else:

             degree= max(len(X), m+len(B))
             temp = [np.zeros(dim, dtype=int)]*degree
             Btemp=B
             B=[np.zeros(dim, dtype=int)]*m + B
             if len(X)>len(B):
                 B=B+[np.zeros(dim, dtype=int)]*(len(X)-len(B))
             else:
                 X=X+[np.zeros(dim, dtype=int)]*(len(B)-len(X))

             X=[(X[i]- (np.multiply(matrixPow(b, F-2, F), np.multiply(d, B[i]) ) ) )%F for i in range(len(X))]
             normalize(X)
             normalize(B)


             #C(x) = C(x) - d b^{-1} x^m B(x);
             B=Btemp
             m = m + 1
    print(L)
    print(len(X))
    if(L+1> len(X)):
        return [-1]
    error_locations=[]
    for i in range(length):
        #print("Evaluated at "+str(i)+"^-1")
        X_alpha=eval(pow(alpha, length-i, F),X, F)
        if np.any(X_alpha==0):
            error_locations.append(i)
    if len(error_locations)+1 != len(X):
        print(error_locations)
        return [-1]
    else:
        return error_locations


    #return (X,L,d)
error_locations=berlekamp(length, msgLen,t, 65537, 3 )
print(error_locations)
#print(d)
#print("error locator berlekamp")
#rint(X)
print("error locations")
print(error_locations)

#  return L
