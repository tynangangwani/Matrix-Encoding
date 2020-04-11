import torch

def polynomial_encode_py(matA, matB, totalWorkers, divA, divB, alpha, beta):
    #takes in matrix A and B, for the computation C=A^t*B, divA and divB
    # correspond to the number of parts A and B will be divided into and
    #alpha and beta are the encoding parameters, usually alpha=1 and beta=divA

    m=divA
    n=divB

    splitsA = torch.split(matA, int(matA.size()[0]/divA), dim=0)
    splitsB= torch.split(matB, int(matB.size()[0]/divB), dim=0)

    submatricesA=[]
    submatricesB=[]

    for i in range(totalWorkers):
        if i==0:
            submatricesA.append(splitsA[i])
        else:
            for j in range(m):
                if j==0:
                    Ai=torch.mul(splitsA[j],(i**(j*alpha)) )
                else:
                    torch.add(Ai, torch.mul(splitsA[j], (i**(j*alpha)) ), out=Ai )
            submatricesA.append(Ai)
    for i in range(totalWorkers):
        if i==0:
            submatricesB.append(splitsB[i])
        else:
            for j in range(n):

                if j==0:
                    Bi=torch.mul(splitsB[j],(i**(j*beta)) )

                else:
                    torch.add(torch.mul(splitsB[j],(i**(j*beta)) ), Bi, out=Bi)
            submatricesB.append(Bi)

    return submatricesA, submatricesB
