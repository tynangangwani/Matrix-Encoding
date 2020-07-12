from nprime.pyprime import miller_rabin
import math
import sympy

# With n the number you want to t

def is_prime(n):
    if n <= 3:
        return n > 1
    elif (n%2 == 0) or (n%3 == 0):
        return False
    i=5
    while i*i <= n:
        if n%i==0 or n%(i + 2) == 0:
            return False
        i+= 6
    return True

def primeFactors(n):

    # Print the number of two's that divide n
    while n % 2 == 0:
        print (2)
        n = n / 2

    # n must be odd at this point
    # so a skip of 2 ( i = i + 2) can be used
    for i in range(3,int(math.sqrt(n))+1,2):

        # while i divides n , print i ad divide n
        while n % i== 0:
            print (i)
            n = n / i

    # Condition if n is a prime
    # number greater than 2
    if n > 2:
        print (n)

primes=[]
a=0
'''
for n in range(int((2**(24.5)))-(2**16)+1, int(2**24.5)+2, 2):
    if miller_rabin(n):
        if ((n-1)%(2**8)==0):
            primes.append(n)
            print(n)
print (primes)
'''
F=23725313# 2147483137
for i in range(F):
    if(math.gcd(i,F==1)):
        if(sympy.ntheory.residue_ntheory.is_primitive_root(i, F)):
            print(i)
            if i>10:
                break
