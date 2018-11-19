import math
import numpy
from numpy.linalg import svd

#ALTERNATING CATENATION

#--------------EXAMPLE fdatavec----------------------------------------------------------------------------------------------------------------------------------
#Comment this stuff out when no longer using it 
#fdatavec = range(32) #Use this since it's length 2^b so basic structure of the algorithm holds 
fdatavec = [1]*1024 #Makes it a list of 32 1's. similarly [1]*4 = [1,1,1,1] Used as compressibility comparison 
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
        W[(2*n)] = M[n][::2] #Should it be an alternation every 2nd pt or should we do W1[0] = A[0][::Ngenes] ???????
        W[((2*n)+1)] = M[n][1::2]
    
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





#print range(16)
#print range(16,32)
#A=[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],[16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]]

#A = [[0, 1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14, 15]]
#b=3

'''W1 = [None]*4
W1[0] = [A[0][i] for i in range(0, int(2**(b-1)))]
W1[1] = [A[1][i] for i in range(0, int(2**(b-1)))]
W1[2] = [A[0][i] for i in range(int(2**(b-1)), int((2**b)))]
W1[3] = [A[1][i] for i in range(int(2**(b-1)), int((2**b)))]'''
'''W1[4] = [A[0][i] for i in range(0, int(2**(b-1)))]
W1[5] = [A[1][i] for i in range(0, int(2**(b-1)))]
W1[6] = [A[0][i] for i in range(int(2**(b-1)), int((2**b)))]
W1[7] = [A[1][i] for i in range(int(2**(b-1)), int((2**b)))]'''
'''W1 = [None]*4
W1[0] = [A[0][::2] for i in range(0, int(2**(b-1)))]
W1[1] = [A[0][1::2] for i in range(0, int(2**(b-1)))]
W1[2] = [A[1][::2] for i in range(0, int(2**(b-1)))]
W1[3] = [A[1][1::2] for i in range(0, int(2**(b-1)))]'''

#-----------------------------------------------------------Alt Catenation algorithm -------------------------------------------------------
'''W1 = [None]*4
W1[0] = A[0][::2] #Should it be an alternation every 2nd pt or should we do W1[0] = A[0][::Ngenes] ???????
W1[1] = A[0][1::2]
W1[2] = A[1][::2]
W1[3] = A[1][1::2]'''

'''for n in range(NumofRows):
    W[(2*n)] = M[n][::2]
    W[((2*n)+1)] = M[n][1::2]'''
#-----------------------------------------------------------Alt Catenation algorithm -------------------------------------------------------

    
'''print W1
#print W1[1::2]

print "A is given by:"
print A
print A[1][::2]
print A[1][::1]
print range(10)[0:4]
print range(10)[4:8]'''





