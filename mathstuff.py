import math
import numpy as np
import numpy.random

#x = input = len(DNAdatavector)
print float(math.log( 2 )) / math.log( 1024 ) # gives the base 2 logarithm of 1024 

print math.log( 1024 ) / float(math.log( 2 )) # gives the reciprocal of the base 2 logarithm of 1024
# we can let b = math.ceil(math.log( x ) / float(math.log( 2 )))      <--- makes be the 
#Then if 2**b = x, then x is a power of 2 already. If not then add datavector to a list of length (2**b - x) filled with whatever the filler variable is 

print math.ceil(math.log( 1024 ) / float(math.log( 2 ))) # ceiling of base 2 log of 1024
print math.ceil(math.log( 1025 ) / float(math.log( 2 ))) # ceiling of base 2 log of 1025, different b/c this one slightly > 10 
print 2**10
print math.log( 42 )
print (10 - 20) / float((100 - 10))

print 0*numpy.arange(17) # creates a list of 17 0's
print list((0*numpy.arange(17) + 4)) #creates a list of 17 4's, we will use modifications of this this as the filler vector 

print numpy.linalg.svd( [
    [2, 5, 3],
    [1, 2, 1],
    [4, 1, 1],
    [3, 5, 2],
    [5, 3, 1],
    [4, 5, 5],
    [2, 4, 2],
    [2, 2, 5],
]) #Does an SVD 

print numpy.linalg.svd( [
    [2, 5, 3],
    [1, 2, 1],
    [4, 1, 1],
    [3, 5, 2],
    [5, 3, 1],
    [4, 5, 5],
    [2, 4, 2],
    [2, 2, 5],
])[1]
#prints only the second part of the SVD, NOTE: in python A[i] gives the (i+1)-th part of the matrix or list A. so the above gives the 2nd part of the SVD, aka the singular values 

print numpy.linalg.svd( [
    [2, 5, 3, 1, 2],
    [1, 2, 1, 6, 2],
    [4, 1, 1, 1, 5],
    [3, 5, 2, 1, 2],
    [5, 3, 1, 3, 2],
    [4, 5, 5, 1, 2],
    [2, 4, 2, 2, 2],
    [2, 2, 5, 1, 2],
])

print range(0,2**4)
print range(2**4,(2**5))
print range(0,40,3)
#print range(0,7036)
print int(2.0)

print int(15 % 4) #Prints the modulo, i.e. 15 % 4 == 15 mod 4 in Python

print range(0,70) #This is a list 
print numpy.linspace(0.0, 69.0, 70.0) #This is technically a set of scalars, use linspace to structure piecewise fns. Also gives 0 - 69 (70 numbers) in 70 equally spaced intervals, which is just the sequence 0 to 69 in this example 

print numpy.linspace(0.0, 7035.0, 7036.0)

A = [range(0,10),range(10,20)]
y = numpy.linspace(0.0, 19.0, 20.0)

print (A[1])[2]
print numpy.isscalar((A[1])[2])

print numpy.piecewise(y, [(y%2) == 0.0, (y%2) == 1.0], [0,1] )

B = numpy.asarray(A)
#print B
NumofRows = 4
print math.sqrt(NumofRows**3)

#rho = 0*numpy.ndarray(shape=(math.sqrt(NumofRows**3),math.sqrt(NumofRows**3)))
#rho = numpy.ndarray.tolist(0.0*numpy.ndarray(shape=(math.sqrt(NumofRows**3), math.sqrt(NumofRows**3))))
#r=0*numpy.ndarray(shape=(8, 8))
rho = [[0]*(8)]*8

print rho

print float(2)
print int(2.0)
print int(3)



NumSampleMin = 3
NumSampleMax = 6
Ngenes = (NumSampleMax + 1) - NumSampleMin
print Ngenes
print range(NumSampleMin, (NumSampleMax + 1) )
print range(1,Ngenes+1)

print range(2,15 +1)

print (math.log(int(12)) ) / float(math.log( 2 ))
print ((1-12))**(-1)

#print list(numpy.loadtxt('PCV1_Sample_' + str(2) + '.txt'))
#print list(numpy.loadtxt('OdysseyBook1Fagle.txt'))

def reverse_numeric(x, y):
    return y - x

aaa=range(10)

print aaa
print sorted(range(10), cmp=reverse_numeric)

a = numpy.asarray(range(10))
print numpy.asarray([1,2]), "As an array"
print a, "a As an array"

a[ abs(a) < 5.0 ] = 0 #Sets a threshold so that any value with abs(a) < 0.0000000001 just goes to 0 

#print a
print numpy.random.randint( 4, size=10 )

mytestlist=[1,2,3,4,5]

print mytestlist
del mytestlist[-1] #deletes last value of mytestlist
print mytestlist

mytestlist2 = range(10)

print mytestlist2

del mytestlist2[-3:-1] #deletes 3rd to last to the 2nd to last value of mytestlist2
del mytestlist2[-1] #deletes last value of mytestlist2
print mytestlist2

mytestlist3 = range(1,24)

print mytestlist3
print 'The excised terms are ',mytestlist3[-8:-1], "and", mytestlist3[-1]

del mytestlist3[-8:-1] #deletes 3rd to last to the 2nd to last value of mytestlist2
del mytestlist3[-1] #deletes last value of mytestlist2
print mytestlist3


mytestlist4 = range(1,31)
print mytestlist4
Ngenes=3
x0= len(range(1,31))
genelength= x0/Ngenes
print genelength
print 'The excised terms are ', mytestlist4[(genelength - 1):: genelength ]
del mytestlist4[(genelength - 1):: genelength ]
print mytestlist4

#testnum = 4
#print testnum / ( math.floor(testnum) )
#if (testnum / ( math.floor(testnum) )) != 1.0:
#    print 'WARNING: genes not of same length, proceed with extreme caution or excise/add some terms. This can break/interfere with many of the other algorithms. '

'''genelength = 10
cutgenelength = 5 
Ngenes = 10
catdatavec = (range(1,101))
x0=len(catdatavec)
print catdatavec
del catdatavec[(cutgenelength*Ngenes) : x0+1 ] 
#del catdatavec[-1] #deletes the last element of each gene (Assuming genes are of same length, if they aren't a warning will be printed above)
print catdatavec
print len(catdatavec)

#food = 'bread'
#vars()[food] = 123
#print bread  # --> 123
c = 14
stringRepr = 'Lamba' + str(c)'''


matrixsample=[[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]

'''print matrixsample
print matrixsample[0]
print matrixsample[0][0]'''

LambdaMat = [[73, 0, 0, 0, 0, 0, 0, 0], [0, 15, 0, 0, 0, 0, 0, 0]]

#print "The number of singular values, i.e. the bond dimension, is ", len(Lambda)

M = [None]*len(LambdaMat)
for n in range(len(LambdaMat)):
    M[n] = matrixsample[n]
    
#print M[17]

a = [[1, 0], [0, 1]]
b = [[4, 1], [2, 2]]

#print numpy.dot(a, b)

AA=[[82,92,102,112,122,132,142,152],[92,104,116,128,140,152,164,176],[102,116,130,144,158,172,186,200],[112,128,144,160,176,192,208,224],[122,140,158,176,194,212,230,248],[132,152,172,192,212,232,252,272],[142,164,186,208,230,252,274,296],[152,176,200,224,248,272,296,320]]

BB=[[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1]]
CC=[[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1]]

#print numpy.dot(AA, BB, CC)

list1 = [1,2,3,4]


#print numpy.diag(list1)

Lambda = [73,15,4,32]
V = AA

#LL = [ ([0]*len(V) )]*len(Lambda)

LL= [None]*len(Lambda)
for i in range(len(Lambda)):
    LL[i] = [0]*len(V)

for n in range(len(Lambda)):
        LL[n][n] = Lambda[n]

print LL



'''print Lambda[0]
print LL
print LL[1]

LL[0][0] = Lambda[0]
LL[1][1] = Lambda[1]'''


print 1000e-10

AAmatrix = [[1,2,3], [4,5,6], [7,8,9] ]
#print AAmatrix[0][2]
#threshtestmatrix[ abs(threshtestmatrix) < 4 ] = 0

#print threshtestmatrix

print range(10)
print np.clip(range(10),0,7)

threshtestlist = [10000, 9999, 9998, 9997, 9996, 9995, 4, 3, 2, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, 0.0000000001]
print threshtestlist

for i in range(len(threshtestlist)):
    if threshtestlist[i] < abs(threshtestlist[0]) * (1e-09) :
        threshtestlist[i] = 0
    else:
        threshtestlist[i] = threshtestlist[i]
        
print threshtestlist
        


#print np.sqrt(2)