
def polynomial_decode(results, totalWorkers, workerIDs, m,n, alpha, beta ):
    #m,n are the divisions of the matrices
    #results are a list of m*n vectors with the order of worker ID's according to what each worker returned
    #so lenght of workerIDs should be m*n
    #worker ID's are a list of worker IDs that returned values
    encoding_matrix=np.empty((m*n, m*n) )
    print(workerIDs)
    for i in range(len(workerIDs)):
        for j in range(m):
            for k in range(n):
                encoding_matrix[i,j*m+k]=workerIDs[i]**(j*alpha+k*beta) #need to  be careful about the order of j and m


    inverse_encoding=np.linalg.inv(encoding_matrix)
    print (inverse_encoding)
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
