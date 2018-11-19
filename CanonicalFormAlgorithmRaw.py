import math
import numpy 
from numpy.linalg import svd
 
movieRatings = [
    [2, 5, 3],
    [1, 2, 1],
    [4, 1, 1],
    [3, 5, 2],
    [5, 3, 1],
    [4, 5, 5],
    [2, 4, 2],
    [2, 2, 5]
]
 
 
#Have an MPS    M_N  .... M_2  M_1   * Z_1      (Z_1 = 1   is a placeholder)

#Matrix M_k 

b= 4 
M = numpy.ndarray.tolist( numpy.random.randint( 4, size= b ) ) 
Z = numpy.ndarray.tolist( numpy.random.randint( 4, size= b ) ) 
#^ Just randomly filled place holder vectors
B = [None]*int(b) 


B[1] = numpy.dot( M[1] , Z[1] )

k=1  #later iterable


#=========================Orthonormalization=======================================================================================================================================
#------------------------------------start of SVDs-------------------------------------------------------------------------------------------------------------
U, singularValues, A = svd( B[1] )  #U[k+1], singularValues[k], A[k] = svd( B[k] )     #for later iteration
#NOTE: svd(movieRatings) outputs some part1,part2,part3, so the equal sign is assigned to all terms, U singularValues, and V, not just V alone 
#this defines U as output 1, singularValues as output 2, V as output 3 
#V is the main 


DELTA=0*numpy.ndarray( shape= (len(U),len(A)) )    
for i in range(len(singularValues)):
    DELTA[i][i]= singularValues[i]
    
matrixDELTA = numpy.matrix(DELTA) #This is now the DELTA matrix of proper size (although maybe rows and columns are switched and it's a 3x8, not an 8x3 as it should be)

#DELTA just of temporary use, no need to store each one. Same goes for other matrices that go into the Z. But DO need to store LAMBDA



Q, singularValues, P = svd( B[1] ) 

LAMBDA=0*numpy.ndarray( shape = (len(Q),len(P)) )    
for i in range(len(singularValues)):
    LAMBDA[i][i]= singularValues[i]
    
matrixLAMBDA = numpy.matrix(LAMBDA)

Mtilde=M #??? right shape?
Mtilde[ 1 ] = numpy.dot( P[1] , A[1] )

Z[ 2 ] = numpy.dot( U[2] , Q[1] )  #dot LAMBDA ???


#=========================Canonization=======================================================================================================================================
#------------------------------------start of SVDs-------------------------------------------------------------------------------------------------------------
U, singularValues, A = svd( B[1] )  #U[k+1], singularValues[k], A[k] = svd( B[k] )     #for later iteration
#NOTE: svd(movieRatings) outputs some part1,part2,part3, so the equal sign is assigned to all terms, U singularValues, and V, not just V alone 
#this defines U as output 1, singularValues as output 2, V as output 3 
#V is the main 


DELTA=0*numpy.ndarray( shape= (len(U),len(A)) )    
for i in range(len(singularValues)):
    DELTA[i][i]= singularValues[i]
    
matrixDELTA = numpy.matrix(DELTA) #This is now the DELTA matrix of proper size (although maybe rows and columns are switched and it's a 3x8, not an 8x3 as it should be)

#DELTA just of temporary use, no need to store each one. Same goes for other matrices that go into the Z. But DO need to store LAMBDA



Q, singularValues, P = svd( B[1] ) 

LAMBDA=0*numpy.ndarray( shape = (len(Q),len(P)) )    
for i in range(len(singularValues)):
    LAMBDA[i][i]= singularValues[i]
    
matrixLAMBDA = numpy.matrix(LAMBDA)

Mtilde=M #??? right shape?
Mtilde[ 1 ] = numpy.dot( P[1] , A[1] )

Z[ 2 ] = numpy.dot( U[2] , Q[1] )  #dot LAMBDA ???

CanonLam=Z

#====================================================Trunctation==================================
#Trunctuate 1 at a time then reconstruct
lasteigen = len( numpy.nonzero( numpy.diagonal(CanonLam) ) )

CanonLam[lasteigen] = 0.0

'''print numpy.ndarray.tolist(0*numpy.matrix(r))
#for j in range(len(Lambda)):
#    R=[None]*(len(V))
#    R[j]=r[j]

#print len(U[1])
#print list(0*numpy.arange(len(V)))

#print 0*numpy.ndarray(shape=(len(U),len(V))) #not necessarily sure which is the row and which is column, depends if rows horizontal as usual or not 

#AAA=0*numpy.ndarray(shape=(len(U)-len(Lambda),len(V)))
#print AAA
#print AAA
#print r
#print range(3,8)

##print numpy.matrix(AAA+r)'''