import numpy as np
import random
import threading
import time
from fft import polycode_ifft
from fft import fft
import math

def normalize(poly):
    while poly and poly[-1] == 0:
        poly.pop()
    if poly == []:
        poly.append(0)
    return poly

def euclidianDivision(a, b, F=65537, prim_root=3):
    normalize(b)
    normalize(a)
    q=[0]*(len(a)+1-len(b))
    r=a
    d=len(b)
    c=b[len(b)-1] #first non-zero
    while (len(r)>=d and not(r[0]==0 and len(r)==1)):

        s= r[len(r)-1]*pow(c,F-2,F) %F #multiplied by x^(deg(r)-d)
        q[len(r)-d]=(q[len(r)-d]+s) %F
        r=(r[:(len(r)-len(b))])+[(r[i+(len(r)-len(b))]-s*b[i]) %F for i in range(len(b))]
        normalize(r)
        #print(q,r, b)
    return q,r

def polyMultiply(a,b, F=65537, prim_root=3): #uses fft to multiply polynomials in (nlogn)
    if(len(a)==0 or len(b)==0):
        return [0]
    if(a[0]==0 and len(a)==1) or (b[0]==0 and len(b)==1):
        return [0]
    if(2**int(math.log2(len(a))) !=len(a)):

        a=a+(([0]*(2**(int(math.log2(len(a)))+1) -len(a))))


    if(2**int(math.log2(len(b))) !=len(b)):
        b=b+(([0]*(2**(int(math.log2(len(b)))+1) -len(b))))


    if(len(a)>len(b)):
        b=b+[0]*(len(a)-len(b))
    else:
        a=a+[0]*(len(b)-len(a))
    a=a+[0]*len(a)
    b=b+[0]*len(b)

    fft_a=fft(a, F, prim_root)
    fft_b=fft(b, F, prim_root)
    prod=[(a[i]*b[i])%F for i in range(len(a))]
    c=polycode_ifft(prod, F, prim_root)
    normalize(c)
    return c

def extended_gcd(a, b, gcd_size, F=65537, prim_root=3): #for as+bt=gcd(a,b), returns gcd,s
#up to the size of gcd wanted
    s=[0]
    old_s = [1]
    r = b
    old_r = a

    while len(old_r)>=(gcd_size+1):
        (q,tempR)= euclidianDivision(old_r, r, F, prim_root)
        old_r=normalize(r)
        r=normalize(tempR)
        tempS=s
        sxq=polyMultiply(q , s,F, prim_root)
        if(len(sxq)>len(old_s)):
            old_s=old_s+[0]*(len(sxq)-len(old_s))
        else:
            sxq=sxq+[0]*(len(old_s)-len(sxq))

        s=[(old_s[i] -sxq[i])%F for i in range(len(sxq)) ]
        old_s=tempS

    #print("bezout s:")
    #print(old_s)
    #print("greatest common divisor:")
    #print(old_r)
    return old_r,old_s

#print(extended_gcd([8,4,2,2,4,5], [3,5,4, 8, 2], 65537, 3) )
def decodeRS(n,k,Crtn, F=65537, prim_root=3):

    G0=[-1%65537]+[0]*(n-1)+[1]
    G1=polycode_ifft(coded_error, 65537, 3)
    #print(G1)
    #correct_msg=polycode_ifft(coded, 65537, 3)
    g,v=extended_gcd(G1,G0, int((n+k)/2))
    #print(v)
    corr,r=euclidianDivision(g, v, 65537, 3)
    if(len(corr)<=k and r==[0]):
        return corr
    else:
        print(corr,r)
        return "error decoding"
'''

n,k=16,8
for i in range(1):

    msg=[random.randint(0,65536) for i in range(k)]+[0]*(n-k)
    #print(msg)
    Correct_msg=msg.copy()
    coded=fft(msg, 65537, 3)

    coded_error=coded.copy()
    for i in range(int((n-k)/2)):
        error_loc1=random.randint(0,n-1)
        coded_error[error_loc1]=(coded[error_loc1]+random.randint(0,65536))%65537

    corr=decodeRS(n,k,coded_error,65537, 3)
    if(normalize(Correct_msg)!=corr):

        print("error.")
        print(Correct_msg)
        print (corr)
'''
'''
F=65537
prim_root=3
for i in range (20):
    a=[random.randint(0,10) for i in range(random.randint(3,15))]

    b=[random.randint(0,10) for i in range(random.randint(3,15))]
    c=polyMultiply(a,b)

    a_,r=euclidianDivision(c,b)
    normalize(a_)
    if a==a_:
        print("correct")
    else:
        print("error:" )
        print(a)
        print(a_)
'''
