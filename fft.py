import math


def fft(Crtn, F, prim_root): #inverse using n^-1 * sum(y*Wn^(-jk))
    length=len(Crtn)
    rt=pow(prim_root, int((F-1)/length), F)
    var=([pow(rt, i, F) for i in range(length)])
    Copy=Crtn.copy()
    for i in range(length): #bit reverse
        format_str='{:0'+str(int((math.log(length,2))))+'b}'
        bit_rev=int(format_str.format(i)[::-1], 2)

        Crtn[i]=Copy[bit_rev]

    for s in range(1,int((math.log(length,2)))+1): #generic fft algorithm
        m = pow(2, s)
        #Wm= pow(3,int((F-1)*(m-1)/m), F)
        #W= 1
        #either W, Wm can be used or var[(-j * int(length/m))%length]
        for j in range(int(m/2)):

            for k in range(j, length, m) :
                t= (var[(j * int(length/m))%length]*Crtn[(k+int(m/2))])%F #
                u=Crtn[k]%F
                Crtn[k]=(u+t)%F
                Crtn[(k+int(m/2))]=(u-t)%F

            #W=(Wm*W) %F
#    for i in range(length): #multiplying by the inverse of n
#        Crtn[i]=(Crtn[i]*pow(length,F - 2, F))%F
    return Crtn
def polycode_ifft(Crtn, F, prim_root): #inverse using polycode method
    length=len(Crtn)
    rt=pow(prim_root, int((F-1)/length), F)
    var = [pow(rt, i, F) for i in range(length)]
    two_inv= pow(2,F-2, F)


    for k in range(int((math.log(length,2)))):
        jump = int(length/pow(2,k))
        #Wm= pow(rt,int(length*(jump-1)/jump), F)
        #W= 1
        for i in range(int(jump/2)):
            block_num = 2**k #or 8/jump
            #for j in range(block_num):
            for base in range(i, length, jump): #or (i, 2*jump*block_num)
                #base = i + j * jump * 2
                Crtn[base] = ((Crtn[base] + Crtn[base + int(jump/2)]) * two_inv) % F
                Crtn[base + int(jump/2)] = ((Crtn[base] - Crtn[base + int(jump/2)]) * var[(-i * int(block_num))%length]) % F

            #W=(Wm*W) %F
    #Copy=Crtn.copy()
    '''
    for i in range(length): #bit reverse
        format_str='{:0'+str(int((math.log(length,2))))+'b}'
        bit_rev=int(format_str.format(i)[::-1], 2)

        Crtn[i]=Copy[bit_rev]
        '''

    return Crtn
'''
#print(fft([10,510,65535,65023], 4294957057, 10)) #evaluations of a degree 3 polynomial with coefficients 1,2,3,4 at 4th roots of unity
t=fft([3453249, 4106606894,439307577, 2356111468,4142483687, 22222222, 0, 0], 4294957057, 10)
print(t)
#t[5]=t[5]+1
#t[6]=t[6]+1
#t[7]=t[7]+1
print(polycode_ifft(t ,4294957057, 10))
#print(polycode_ifft([10,510,65535,65023], 65537, 3))
#10,510,65535,65023
#1,2,3,4 #expected resultss
'''
