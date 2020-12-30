import numpy as np
import random
import threading
import time
from fft import polycode_ifft
from fft import fft
import math
from DecodingRSVectorized import decodeRS

t=[1, 2, 3 ,4, 0 ,0,0,0 ]
t=fft(t,65537,3)
#t=fft([3453249, 4106606894,439307577, 2356111468,4142483687, 22222222, 0, 0], 4294957057, 10)

print("start")
print(t)
[5]=t[5]+5%65537
t[6]=t[6]+234%65537
#print(polycode_ifft(t, 65537, 3))

t=fft(t, 65537, 3)
print(t)
L = 0
m = 1
b = 1
N=4 #number of syndromes
C=t[1:5]
print("c: ", C)
  #k=
  #/* steps 2. and 6. */
  #for (n = 0; n < N; n++) {
'''for n in range(n):
      #/* step 2. calculate discrepancy */
      d = s_n + \Sigma_{i=1}^L c_i * s_{n-i};

      if (d == 0) {
          /* step 3. discrepancy is zero; annihilation continues */
          m = m + 1;
      } else if (2 * L <= n) {
          /* step 5. */
          /* temporary copy of C(x) */
          polynomial(field K) T(x) = C(x);

          C(x) = C(x) - d b^{-1} x^m B(x);
          L = n + 1 - L;
          B(x) = T(x);
          b = d;
          m = 1;
      } else {
          /* step 4. */
          C(x) = C(x) - d b^{-1} x^m B(x);
          m = m + 1;
      }
  }
  return L;
'''
