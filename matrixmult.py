import math
import numpy 

# Program to multiply two matrices using nested loops

# 3x3 matrix
X = [[1,2,3],
    [4 ,5,6],
    [7 ,8,9]]
    
vec = [1,2,3]
# 3x4 matrix
Y = [[5,8,1,2],
    [6,7,3,0],
    [4,5,9,1]]
    
XTranspose=[[1, 4, 7], [2, 5, 8], [3, 6, 9]]

    
# result is 3x4
result = [[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0]]

#Y=numpy.ndarray(shape=(3,4))
#print Y

#In this example X is on left and Y is on the right 
# iterate through rows of X
for i in range(len(X)): 
   # iterate through columns of Y
   for j in range(len(Y[0])):
       # iterate through rows of Y
       for k in range(len(Y)):
           result[i][j] += X[i][k] * Y[k][j]



print result
#"By Running this for X and XTranspose we get different results. By checking these results with some from mathematica we are assured of the structure of matrices in Python, namely that they're the same as in mathematica row/column wise"
print "Structured same as mathematica: [[row1column1,row1column2,row1column3], [row2column1,row2column2,row2column3], [row3column1,row3column2,row3column3]]"

print X
print X[1]
print X[1][1]

print vec
print vec[1]


fdatavec = range(16) #Use this since it's length 2^4 so basic structure of the algorithm holds 
x=len(fdatavec)
b = math.ceil(math.log( x ) / float(math.log( 2 )))
M=[fdatavec]

print M
print M[0][3]

print "Some stuff for ExampleSVD.py program: "
print (2*(2**(0)))/2
print range((2*(2**(0)))/2)
print (2*(2**(1)))/2
print range((2*(2**(1)))/2)
print (2*(2**(2)))/2
print range((2*(2**(2)))/2)

#print 7.%3. #mod
#print 7./3. #division 

print 4.0 - 4
print (2**(4.0-3))
print [None]*2
#print [None]*(2.0) #This one produces an error since 2.0 is a float and need it to be an integer here, but [None]*int(2.0) works 

someNumber = 4 
print "this is a four:", someNumber, "isn't it?"




#print len(X)
#print list(0*numpy.arange(len(X)))
#print numpy.matrix(list(0*numpy.arange(3)))

#for r in result:
#   print(r)