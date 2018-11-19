import math
import numpy as np
import numpy.random

x = 4.0

print 0*range( int(x) )

print 2*range( int(x) )

print int(x)*[0]

#print np.sqrt(2)


#print [0.0, 1.0, 2.0, 3.0] 
#print numpy.array([0.0, 1.0, 2.0, 3.0]) + 1.0

#print numpy.ndarray.tolist(numpy.array([0.0, 1.0, 2.0, 3.0]) + 1.0)

#print range(12-1,0,-1)


#print np.ones((5, 5))


Lambda = [1,2,3]
#a= numpy.ndarray.tolist( numpy.diag([1,2,3,0,0]) )

#b= numpy.ndarray.tolist( np.ones((5, 5)) )

#print a
#print b

#c = numpy.dot( a,b )

#print c

#print numpy.resize(c, ( len(Lambda), len(c) ) )

a = [0,1,2,3,4,None]
a.remove(None)

print a

myList = [10,20,30,40,50,60,70,80,90,100]

print myList

myInt = 10
#newList = [x / myInt for x in myList]

myList[:] = [x / myInt for x in myList]

print myList
#print newList

print float(1.0)/6248

print numpy.random.randint( 1, 5, size=16 ) # makes a list of 16 random integers with values 1, 2, 3, or 4

print np.prod([1,2,3,4])


a =[1,2,3,4]

a=map(float, a)

print a

A = [[1,2,3,4], [5,6,7,8], [9,10,11,12], [13,14,15,16]]

print A[0][0]
print A[0]

def column(matrix, i):
    return [row[i] for row in matrix]

print column(A,0)

print A

print A[0:3]

print A[0][0:( len(A) / 2 )]
print A[1][0:( len(A) / 2 )]

#print numpy.concatenate(( [1, 2] , [5, 6]), axis=0)


print np.bmat([[[A[0]], [A[1]]], [[A[2]], [A[3]]]])

AAA = [ [0, 0] ]
print np.bmat( [ [ AAA ], [ [[1, 2]] ] ] )

for i in range( 1, 10 ):
    AAA = np.bmat( [ [ AAA ], [ [[i, i]] ] ] )

print a 
print a[ 0 : (len(a)/2) ]
print a[ (len(a)/2) : (len(a)) ]

print A
print numpy.ndarray.tolist( np.bmat( [  [[ A[0][0:( len(A) / 2 )] ]], [[ A[1][0:( len(A) / 2 )] ]]  ] ) )
print numpy.ndarray.tolist( np.bmat( [  [[ A[2][0:( len(A) / 2 )] ]], [[ A[3][0:( len(A) / 2 )] ]]  ] ) )
print numpy.ndarray.tolist( np.bmat( [  [[ A[0][( len(A) / 2 ) : len(A)] ]], [[ A[1][( len(A) / 2 ) : len(A)] ]]  ] ) )
print numpy.ndarray.tolist( np.bmat( [  [[ A[2][( len(A) / 2 ) : len(A)] ]], [[ A[3][( len(A) / 2 ) : len(A)] ]]  ] ) )

print '____________________________________________________________________________'

#[( 0 + (1*j) ):( len(A) / 2 )]
for i in range( 0, ( len(A) / 2 ) ):
      print numpy.ndarray.tolist( np.bmat( [  [[ A[(2*i)][0:( len(A) / 2 )] ]], [[ A[(2*i + 1)][0:( len(A) / 2 )] ]]  ] ) )
      print numpy.ndarray.tolist( np.bmat( [  [[ A[(2*i)][( len(A) / 2 ) : len(A)] ]], [[ A[(2*i + 1)][( len(A) / 2 ) : len(A)] ]]  ] ) )

print '_______________incorporating j_____________________________________________________________'

#( ( len(A) / 2 )*j ): ( ( len(A) / 2 )*( (1+j) ) )

for j in range( 0, ( len(A) / 2 ) ):
    for i in range( 0, ( len(A) / 2 ) ):
        print numpy.ndarray.tolist( np.bmat( [  [[ A[(2*i)][( ( len(A) / 2 )*j ): ( ( len(A) / 2 )*( (1+j) ) )] ]], [[ A[(2*i + 1)][( ( len(A) / 2 )*j ): ( ( len(A) / 2 )*( (1+j) ) )] ]]  ] ) )   


Ablock = [[0,1],[2,3]]

for j in range( 0, ( len(A) / 2 ) ):
    for i in range( 0, ( len(A) / 2 ) ):
        Ablock[i][j] = numpy.ndarray.tolist( np.bmat( [  [[ A[(2*i)][( ( len(A) / 2 )*j ): ( ( len(A) / 2 )*( (1+j) ) )] ]], [[ A[(2*i + 1)][( ( len(A) / 2 )*j ): ( ( len(A) / 2 )*( (1+j) ) )] ]]  ] ) )   

#print Ablock

'''
Ablock = [ [None]*( len(A) / 2 ) ]*( len(A) / 2 )
for i in range( 0, ( len(A) / 2 ) ):
    Ablock[i][0] = numpy.ndarray.tolist( np.bmat( [  [[ A[(2*i)][0:( len(A) / 2 )] ]], [[ A[(2*i + 1)][0:( len(A) / 2 )] ]]  ] ) )
    Ablock[i][1] = numpy.ndarray.tolist( np.bmat( [  [[ A[(2*i)][( len(A) / 2 ) : len(A)] ]], [[ A[(2*i + 1)][( len(A) / 2 ) : len(A)] ]]  ] ) )

print Ablock
print ' Ablock[0][0] = ', Ablock[0][0], ' Ablock[0][1] = ', Ablock[0][1], ' Ablock[1][0] = ', Ablock[1][0], ' Ablock[1][1] = ', Ablock[1][1],
'''

print numpy.outer( [[1,2],[3,4]] , [[100,0],[0,0.01]] )

#for j in range( 0, ( len(A) / 2 ) ):
#    for i in range( 0, ( len(A) / 2 ) ):
#        Ablock[i][j] = numpy.ndarray.tolist( np.bmat( [  [[ A[i][0:( len(A) / 2 )] ]], [[ A[j][0:( len(A) / 2 )] ]]  ] ) )

#for i in range( 0, len(A) ):
#    A[i][0:( len(A) / 2 )]

#print numpy.block( A )

