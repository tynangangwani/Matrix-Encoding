import numpy as np
import random
import threading
import time
from fft import polycode_ifft
from fft import fft
import math

#dim is a list of coefficient dimensions
def power_mod(a,p,F):
    return (pow(int(a),int(p),int(F)))
matrixPow=np.vectorize(power_mod)
dim=[2,2]
def normalize(poly):
    while poly and (not poly[-1].any()):
        poly.pop()
    if poly == []:
        poly.append(np.matrix(np.zeros(dim, dtype=int)))
    return poly

def euclidianDivision(a, b, F=65537, prim_root=3):
    normalize(b)
    normalize(a)
    q=[np.matrix(np.zeros(dim, dtype=int))]*(len(a)+1-len(b))
    #print(len(a))
    #print(len(b))
    r=a
    d=len(b)
    c=b[len(b)-1] #first non-zero
    while (len(r)>=d and not((not r[0].any()) and len(r)==1)):
        print("c")
        print(c)
        print("r[len(r)-1]")
        print(r[len(r)-1])
        s= np.multiply(r[len(r)-1], matrixPow(c,F-2,F)) %F #multiplied by x^(deg(r)-d)
        print("s")
        print(s)
        q[len(r)-d]=(q[len(r)-d]+s) %F


        r=(r[:(len(r)-len(b))])+[(r[i+(len(r)-len(b))]-np.multiply(s,b[i])) %F for i in range(len(b))]
        normalize(r)
        #if(len(q)==0 and (not q[0].any())):
        #    return q,r
        print("q")
        print(q)
        print("r")
        print(r)
        print("b")
        print( b)
        if not s.any():
            return "error"
    return q,r

def polyMultiply(a,b, F=65537, prim_root=3): #uses fft to multiply polynomials in (nlogn)
    if(len(a)==0 or len(b)==0):
        return np.matrix(np.zeros(dim, dtype=int))
    if((not a[0].any()) and len(a)==1) or ((not b[0].any()) and len(b)==1):
        return np.matrix(np.zeros(dim, dtype=int))
    if(len(a)>len(b)):
        b=b+[np.matrix(np.zeros(dim, dtype=int))]*(len(a)-len(b))
    else:
        a=a+[np.matrix(np.zeros(dim, dtype=int))]*(len(b)-len(a))
    if(2**int(math.log2(len(a))) !=len(a)):

        a=a+(([np.matrix(np.zeros(dim, dtype=int))]*(2**(int(math.log2(len(a)))+1) -len(a))))
        b=b+(([np.matrix(np.zeros(dim, dtype=int))]*(2**(int(math.log2(len(b)))+1) -len(b))))

    a=a+[np.matrix(np.zeros(dim, dtype=int))]*len(a)
    b=b+[np.matrix(np.zeros(dim, dtype=int))]*len(b)

    fft_a=fft(a, F, prim_root)
    fft_b=fft(b, F, prim_root)
    prod=[(np.multiply(a[i],b[i]))%F for i in range(len(a))]
    c=polycode_ifft(prod, F, prim_root)
    normalize(c)
    return c

def extended_gcd(a,b, gcd_size, F=65537, prim_root=3): #for as+bt=gcd(a,b), returns gcd,s
#up to the size of gcd wanted
    s=[np.matrix(np.zeros(dim, dtype=int))]
    old_s = [np.ones(dim, dtype=int)]
    r = b
    old_r = a

    while len(old_r)>=(gcd_size+1):
        print("r:",r)
        print("old r:", old_r)
        (q,tempR)= euclidianDivision(old_r, r, F, prim_root)
        print("here")
        old_r=normalize(r)
        r=normalize(tempR)
        tempS=s
        sxq=polyMultiply(q , s,F, prim_root)
        if(len(sxq)>len(old_s)):
            old_s=old_s+[np.matrix(np.zeros(dim, dtype=int))]*(len(sxq)-len(old_s))
        else:
            sxq=sxq+[np.matrix(np.zeros(dim, dtype=int))]*(len(old_s)-len(sxq))

        s=[(old_s[i] -sxq[i])%F for i in range(len(sxq)) ]
        old_s=tempS

    #print("bezout s:")
    #print(old_s)
    #print("greatest common divisor:")
    #print(old_r)
    return old_r,old_s

#print(extended_gcd([8,4,2,2,4,5], [3,5,4, 8, 2], 65537, 3) )
def decodeRS(n,k,Crtn, F=65537, prim_root=3):

    G0=[np.full(dim, -1%65537)]+[np.matrix(np.zeros(dim, dtype=int))]*(n-1)+[np.ones(dim, dtype=int)]
    G1=polycode_ifft(coded_error, F, prim_root)
    #print(G1)
    #correct_msg=polycode_ifft(coded, 65537, 3)
    g,v=extended_gcd(G1,G0, int((n+k)/2))
    print("v: ")
    print(v)
    return "error decoding"
    corr,r=euclidianDivision(g, v, F, prim_root)
    if(len(corr)<=k and (not r[0].any())):
        return corr
    else:
        return "error decoding"
        print(corr,r)
n,k=8,4
for i in range(1):

    msg=[np.matrix(np.random.randint(0,10, dim) ) for i in range(k)]+[np.matrix(np.zeros(dim, dtype=int))]*(n-k)
    #print(msg)
    Correct_msg=msg.copy()
    coded=fft(msg, 65537, 3)

    coded_error=coded.copy()
    for i in range(int((n-k)/2)):
        error_loc1=random.randint(1,n-1)
        coded_error[error_loc1]=(coded[error_loc1]+np.random.randint(0,10, dim) )%65537

    corr=decodeRS(n,k,coded_error,65537, 3)
    Correct_msg=normalize(Correct_msg)
    print([np.array_equal(Correct_msg[i], corr[i]) for i in range (len(corr))])
