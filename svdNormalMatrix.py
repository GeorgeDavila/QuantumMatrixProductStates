import math
import numpy 
from numpy.linalg import svd
 
A = [
    [-1, 3, -1],
    [-3, 5, -1],
    [-3, 3, 1]
]
 
U, singularValues, V = svd(A) 
#NOTE: svd(movieRatings) outputs some part1,part2,part3, so the equal sign is assigned to all terms, U singularValues, and V, not just V alone 
#this defines U as output 1, singularValues as output 2, V as output 3 
#V is the main 
Lambda=singularValues

print "The singular values are: ", singularValues #This gives U as an output when you run the file so you can see the SVD stuff works


r= 0*numpy.ndarray(shape=(len(U),len(V))) 
for i in range(len(Lambda)):
    r[i][i]= abs(Lambda[i]) 
    
print "Here's the SVD Lambda matrix: ", numpy.matrix(r) 

#print r[0][2], " ---> Test if that's definitionally equiv to zero: ", r[0][2] == 0

#print numpy.ndarray.tolist(numpy.matrix(r))

#print "U is: ", U
#print "V is: ", V

#print numpy.dot(U,numpy.transpose(V))
print "=====================================^^^SVD stuff^^^====================================="


eigenVals, eigVecs = numpy.linalg.eig(A)

print eigenVals
print len(eigenVals)
diagonalizedForm= 0.0*numpy.ndarray(shape=(len(eigenVals),len(eigenVals))) 
for i in range(len(eigenVals)):
    diagonalizedForm[i][i]= abs(eigenVals[i]) 
    
#print diagonalizedForm
print "Here's the diagonalized form of the matrix: ", numpy.matrix(diagonalizedForm)
#numpy.ndarray.tolist( numpy.matrix(diagonalizedForm) )

print numpy.zeros( (3,3) ) #gives 3 by 3 matrix of all zeros, important for defns 