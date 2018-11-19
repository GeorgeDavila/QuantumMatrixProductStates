import math
import numpy 
from numpy.linalg import svd


A = [
    [1,0,0,0],
    [0,2,0,0],
    [0,0,3,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0]
    ]

print "len(A) = ", len(A)
print "len(A[0]) = ", len(A[0])

print max( len(A), len(A[0]) )

print "Matrix has ", numpy.count_nonzero(A), " non-zero values "
print numpy.diag(A)
print numpy.diag(numpy.diag(A))

print numpy.eye(8,4)
print numpy.eye( 4,8 )

print numpy.nonzero(A)

print numpy.diag(A)
print numpy.count_nonzero( numpy.diag(A) )

NumNonZero=numpy.count_nonzero( numpy.diag(A) )
DIAG= numpy.zeros( ( NumNonZero, NumNonZero ) )

for i in range(NumNonZero):
    DIAG[i][i]= (numpy.diag(A)[i]) 
    
print DIAG


 #-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------
# Given A
NumSingVals= len( numpy.diag(A) ) # = len(A)
NumNonZero=numpy.count_nonzero( numpy.diag(A) ) #= r 
DIAG= numpy.zeros( ( NumNonZero, NumNonZero ) )

for i in range(NumNonZero):
    DIAG[i][i]= (numpy.diag(A)[i]) 
    
print DIAG

diagdelta = DIAG
#Lambda[k] = numpy.diag( numpy.diag( Delta[k-1] ) )

if NumSingVals==NumNonZero:
    P = numpy.eye( len(A), len(A) )
else:
    P = numpy.eye( len(A), NumNonZero )

if NumSingVals==NumNonZero:
    Pinv = numpy.eye( len(A[0]), len(A[0]) )
else:
    Pinv = numpy.eye( NumNonZero , len(A[0]) )
# Gives terms such that Delta[k-1] = P x Lambda[k] x Pdag 
#-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------
