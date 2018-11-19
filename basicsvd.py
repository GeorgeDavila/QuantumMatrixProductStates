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
 
U, singularValues, V = svd(movieRatings) 
#NOTE: svd(movieRatings) outputs some part1,part2,part3, so the equal sign is assigned to all terms, U singularValues, and V, not just V alone 
#this defines U as output 1, singularValues as output 2, V as output 3 
#V is the main 
Lambda=singularValues

print singularValues #This gives U as an output when you run the file so you can see the SVD stuff works
print len(U) #This gives the number of rows (& columns if a square matrix)
print len(singularValues) #Gives number of singular values
print len(V)

print V
#print svd(movieRatings)
#print svd(movieRatings)[1]
print numpy.ndarray.tolist(svd(movieRatings)[2])
#print svd(numpy.ndarray.tolist(svd(movieRatings)[2]))
#print svd(numpy.ndarray.tolist(svd(movieRatings)[2]))[1]
# Note we don't need singularValues written as a matrix since we don't do much further with it anyway, just need the number of sing values 


print singularValues[2]
print numpy.ndarray.tolist(singularValues)[2]
print math.sqrt(singularValues[2])
print math.sqrt(16)
print len(singularValues)

r=0*numpy.ndarray(shape=(len(U),len(V)))
for i in range(len(Lambda)):
    r[i][i]=math.sqrt(Lambda[i])
    
print numpy.matrix(r) #This is now the sqrt(Lambda) matrix of proper size (although maybe rows and columns are switched and it's a 3x8, not an 8x3 as it should be)


print numpy.ndarray.tolist(0*numpy.matrix(r))
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

##print numpy.matrix(AAA+r)