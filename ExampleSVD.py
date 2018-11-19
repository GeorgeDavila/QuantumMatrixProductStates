import math
import numpy
from numpy.linalg import svd


#--------------EXAMPLE fdatavec----------------------------------------------------------------------------------------------------------------------------------
#Comment this stuff out when no longer using it 
#fdatavec = range(32) #Use this since it's length 2^b so basic structure of the algorithm holds 
fdatavec = [1]*32 #Makes it a list of 32 1's. similarly [1]*4 = [1,1,1,1] Used as compressibility comparison 
#And you can see that for  [1]*b  it's sing values decrease very rapidly so it's very compressible as we expect
#But number of sing values doesn't decrease very quickly, so might be that the catenation's off 
x=len(fdatavec)
b = math.ceil(math.log( x ) / float(math.log( 2 )))
M=[fdatavec] #wrap it in brackets to make it a list of lists to make the later algorithm easier 
#--------------EXAMPLE fdatavec---------------------------------------------------------------------------------------------------------------------------------- 

print "b is ", b, "so the length of the combined and filled data sample is ", (2**b)

#iterative B for B in range(int(b),0,-1) = [b, b-1, b-2, ..., 4, 3, 2, 1]
#Only go down to B=1 since we can't do anymore SVD's at that point anyway 

for B in range(int(b),0,-1):
    NumofRows = int((2**(b-B))) #Number of rows of the thing being catenated 
    W =[None]*(2*NumofRows)
    for n in range(NumofRows):
        W[(2*n)] = [M[n][i] for i in range(0, int(2**(B-1)))]
        W[((2*n)+1)] = [M[n][i] for i in range(int(2**(B-1)), int((2**B)))]
    
    # Matrix Multiplication Algorithm to isolate sqrt(Lambda).V of W = U.Lambda.V :
    U=svd(W)[0]
    Lambda=svd(W)[1]
    V=svd(W)[2]
    r=0*numpy.ndarray(shape=(len(U),len(V)))
    
    for i in range(len(Lambda)):
        r[i][i]=math.sqrt(Lambda[i])
        
    L = numpy.ndarray.tolist(numpy.matrix(r)) #This is now the sqrt(Lambda) matrix of proper size (although maybe rows and columns are switched and it's a 3x8, not an 8x3 as it should be)

    result=numpy.ndarray.tolist(0*numpy.matrix(r))

    for i in range(len(L)):
       # iterate through columns of Y
       for j in range(len(V[0])):
           # iterate through rows of Y
           for k in range(len(V)):
               result[i][j] += L[i][k] * V[k][j]

    M = result #Need to redefine M so the for loop acts on result each time 
    print "For the", int(b-B)+1, "-th SVD there are", len(Lambda), "singular values and they are:" #i.e. for W_(b-B) = U.Lambda_(b-B).V, the non-zero values of Lambda_(b-B)
    print Lambda #The list of singular values 
    print "_____________________________________________________________________" #Just to separate the results a bit


#print b
#print range(int(b),0,-1)


#--------Can mostly ignore the stuff below, just shows how and why the above algorithm is structured as it is --------------------------------------------------------------------

#print result #This is given by result = sqrt(Lambda).V, a 2x 2^b matrix, now want to catenate to a 4x2**(b-2) matrix and perform another SVD

#Having trouble writing the exact catenation for W1 defined in the notes
#For now construct a similarly structured 4 x 2^(b-2) matrix and ask Mucciolo for the catenation later in Fortran form. Should be very nearly the same as python form 

'''W1 =[None]*(4) #Note that if we later rig it to do decreasing iterations on b then the intervals are just the same as we defined them to get W, i.e. range(int(2**((b-1)-1)), int((2**(b-1)))) -> range(int(2**(b-1)), int((2**b)))
    W1[0] = [result[0][i] for i in range(0, int(2**(b-1)))]
    W1[1] = [result[0][i] for i in range(int(2**(b-1)), int((2**b)))]
    W1[2] = [result[1][i] for i in range(0, int(2**(b-1)))] #Note change to reading result[1]
    W1[3] = [result[1][i] for i in range(int(2**(b-1)), int((2**b)))]
    
    W2 =[None]*(8) #Note that if we later rig it to do decreasing iterations on b then the intervals are just the same as we defined them to get W, i.e. range(int(2**((b-1)-1)), int((2**(b-1)))) -> range(int(2**(b-1)), int((2**b)))
    W2[0] = [result[0][i] for i in range(0, int(2**(b-1)))]
    W2[1] = [result[0][i] for i in range(int(2**(b-1)), int((2**b)))]
    W2[2] = [result[1][i] for i in range(0, int(2**(b-1)))] #Note change to reading result[1]
    W2[3] = [result[1][i] for i in range(int(2**(b-1)), int((2**b)))]
    W2[4] = [result[2][i] for i in range(0, int(2**(b-1)))]
    W2[5] = [result[2][i] for i in range(int(2**(b-1)), int((2**b)))]
    W2[6] = [result[3][i] for i in range(0, int(2**(b-1)))] #Note change to reading result[1]
    W2[7] = [result[3][i] for i in range(int(2**(b-1)), int((2**b)))]'''
    
    
    
'''M = result 
for n in range(NumofRows):
    W[(2*n)] = [M[n][i] for i in range(0, int(2**(b-1)))]
    W[((2*n)+1)] = [M[n][i] for i in range(int(2**(b-1)), int((2**b)))]'''


#print range(10,0,-1)
#print range(10,-1,-1)
#print range(10,-1,-2)