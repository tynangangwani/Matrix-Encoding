
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


    inverse_encoding=np.linalg.inv(encoding_matrix) #paper says we can find this more efficiently than using traditional inverse
    #algorithm

    #now need to multiply the matrix with the list of results (so each coefficient is multiplied by a matrix)
