import numpy as np
def polynomial_decode(results, totalWorkers, workerIDs, m,n, alpha, beta, coeffMatrix ):
    #m,n are the divisions of the matrices
    #results are a list of m*n vectors with the order of worker ID's according to what each worker returned
    #so lenght of workerIDs should be m*n
    #worker ID's are a list of worker IDs that returned values

    encoding_matrix=coeffMatrix[workerIDs,:]
    encoding_matrix=np.asarray(encoding_matrix)

    inverse_encoding=np.linalg.inv(np.asarray(encoding_matrix))
    #print (inverse_encoding)
    #print(np.matmul(encoding_matrix,inverse_encoding))
    decodeds=[]
    #multiplies the inverse matrix with the matrix list

    for j in range(m*n):
        for k in range(m*n):

            if k==0:
                temp=inverse_encoding[j,k]*results[k]
            else:
                temp+=inverse_encoding[j,k]*results[k]
        decodeds.append(temp)

    return decodeds
