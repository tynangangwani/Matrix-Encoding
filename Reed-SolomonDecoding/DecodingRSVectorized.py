import numpy as np
import random
import threading
import time
from fft import polycode_ifft
from fft import fft
import math

def normalize(poly):

    while len(poly)!=0 and (not poly[-1].any()):
        poly=np.delete(poly, -1)
        #print(poly)

    if len(poly)==0:
        return np.zeros([1], dtype=int)
    return poly

def euclidianDivision(a, b, F=65537, prim_root=3):
    b=normalize(b)
    a=normalize(a)
    q=np.zeros([len(a)+1-len(b)], dtype=int)
    r=a
    d=len(b)
    c=b[len(b)-1] #first non-zero
    while (len(r)>=d) and not(r[0]==0 and len(r)==1):

        s= r[len(r)-1]*pow(int(c),F-2,F) %F #multiplied by x^(deg(r)-d)

        q[len(r)-d]=(q[len(r)-d]+s) %F
        r=np.concatenate(((r[:(len(r)-len(b))]),np.asarray( [(r[i+(len(r)-len(b))]-s*b[i]) %F for i in range(len(b))])), axis=0)
        r=normalize(r)

    return q,r

def polyMultiply(a,b, F=65537, prim_root=3): #uses fft to multiply polynomials in (nlogn)
    if(len(a)==0 or len(b)==0):
        return [0]
    if(a[0]==0 and len(a)==1) or (b[0]==0 and len(b)==1):
        return [0]
    a=np.pad(a,(0,len(a)), 'constant', constant_values=(0, 0))
    b=np.pad(b,(0,len(b)), 'constant', constant_values=(0, 0))
    if(2**int(math.log2(len(a))) !=len(a)):

        a=np.pad(a, (0, (2**(int(math.log2(len(a)))+1) -len(a))), 'constant', constant_values=(0, 0))


    if(2**int(math.log2(len(b))) !=len(b)):
        b=np.pad(b,(0,(2**(int(math.log2(len(b)))+1) -len(b))), 'constant', constant_values=(0, 0))


    if(len(a)>len(b)):
        b=np.pad(b, (0,(len(a)-len(b))), 'constant', constant_values=(0, 0))
    else:
        a=np.pad(a,(0,(len(b)-len(a))), 'constant', constant_values=(0, 0))

    fft_a=fft(a, F, prim_root)
    fft_b=fft(b, F, prim_root)
    prod=[(fft_a[i]*fft_b[i])%F for i in range(len(a))]
    c=polycode_ifft(prod, F, prim_root)

    c=normalize(c)
    return c

def extended_gcd(a, b, gcd_size, F=65537, prim_root=3): #for as+bt=gcd(a,b), returns gcd,s
#up to the size of gcd wanted
    s=np.zeros([1], dtype=int)
    old_s = np.ones([1], dtype=int)
    r = b
    old_r = a

    while len(old_r)>=(gcd_size+1) or (not r.any()):

        (q,tempR)= euclidianDivision(old_r, r, F, prim_root)
        old_r=normalize(r)
        r=normalize(tempR)
        tempS=s
        sxq=polyMultiply(q , s,F, prim_root)


        if(len(sxq)>len(old_s)):
            old_s=np.pad(old_s,(0,(len(sxq)-len(old_s))), 'constant', constant_values=(0, 0))
        else:
            sxq=np.pad(sxq,(0,(len(old_s)-len(sxq))), 'constant', constant_values=(0, 0))

        s=np.asarray([(old_s[i] -sxq[i])%F for i in range(len(sxq)) ])

        old_s=tempS


    return old_r,old_s

#print(extended_gcd([8,4,2,2,4,5], [3,5,4, 8, 2], 65537, 3) )
def decodeRS(Crtn):
    n,k=8,4
    F=65537
    prim_root=3
    G0=np.asarray([-1%65537]+[0]*(n-1)+[1])
    G1=polycode_ifft(Crtn, 65537, 3)
    if not G1[k:].any():
        return normalize(G1)
    #print(G1)
    #correct_msg=polycode_ifft(coded, 65537, 3)
    g,v=extended_gcd(G1,G0, int((n+k)/2))

    corr,r=euclidianDivision(g, v, 65537, 3)
    if(len(corr)<=k and not r.any()):
        return corr
    else:
        return "error decoding"


n,k=8,4
'''for i in range(1):
    dim=100
    msg=np.asarray([np.random.randint(0,65536, (dim,dim)) for i in range(k)]+[np.zeros((dim,dim), dtype=int)]*(n-k))

    #print(msg)

    Correct_msg=msg.copy()
    coded=fft(msg, 65537, 3)
    #print(coded)
    coded_error=coded.copy()
    for i in range(int((n-k)/2)):
        error_loc1=random.randint(0,n-1)
        coded_error[error_loc1]=(coded[error_loc1]+np.random.randint(0,65536,(dim,dim)))%65537
    corr=np.empty((4,2,2))
    coded_error=np.stack(coded_error, axis=-1)
    print("start decoding")
    decode=np.vectorize(decodeRS, excluded=['n', 'k', 'F', 'prim_root'], signature='(n)->(k)')
    corr=decode(coded_error)
    corr=np.swapaxes(corr,0,-1)
    corr=np.swapaxes(corr,2,1)

    #print(corr)
    corr[:,0,1]=decodeRS(n,k,coded_error[:,0,1],65537, 3)
    corr[:,1,0]=decodeRS(n,k,coded_error[:,1,0],65537, 3)
    corr[:,1,1]=decodeRS(n,k,coded_error[:,1,1],65537, 3)
    corr[:,0,0]=decodeRS(n,k,coded_error[:,0,0],65537, 3)
    #Correct_msg=normalize(Correct_msg)
    #print(Correct_msg)
    print([np.array_equal(Correct_msg[i], corr[i]) for i in range (len(corr))])
'''
