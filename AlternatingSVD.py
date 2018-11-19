import math
import numpy
from numpy.linalg import svd

datavec = list(numpy.loadtxt('PCV1_1genesample.txt')) + list(numpy.loadtxt('PCV1_2genesample.txt')) + list(numpy.loadtxt('PCV1_3genesample.txt')) + list(numpy.loadtxt('PCV1_4genesample.txt'))
#Now we want to concatenate the datavector a bit 
gene1 = list(numpy.loadtxt('PCV1_1genesample.txt'))
gene2 = list(numpy.loadtxt('PCV1_2genesample.txt'))
gene3 = list(numpy.loadtxt('PCV1_3genesample.txt')) 
gene4 = list(numpy.loadtxt('PCV1_4genesample.txt'))
#Name them here for ease but would later cut out this step 

Ngenes = 4 
x = len(datavec)

#Beginning of defining cdatavec <---
catdatavec = [None]*(x) #x defined as len(datavec), so defines a list of same length as the datavec (for now filled with None's)

catdatavec[::4] = gene1
catdatavec[1::4] = gene2
catdatavec[2::4] = gene3
catdatavec[3::4] = gene4
#Can catenate before filling it with 4's since 4's end up being at the 

#NOTE: eventually want to structure this in terms of number of gene samples Ngenes used <---

b = math.ceil(math.log( x ) / float(math.log( 2 ))) # this defines b so that 2^b > x, but 2^(b-1) < x, i.e. b is the nearest power of 2 above x (or equal if x is a power of 2)

print "b is ", b, "so the length of the combined and filled data sample is ", (2**b)

filler = list((0*numpy.arange(((2**b) - x)) + 4))

if (2**b) == x:
    fdatavec = catdatavec
elif (2**b) > x:
    fdatavec = catdatavec + filler
else:
    print "Something went wrong. Probably didn't define the data vector as needed. Try again"
#Now fdatavec will be a filled vector of length 2^b

M = [fdatavec]

for B in range(int(b),0,-1):
    NumofRows = int((2**(b-B))) #Number of rows of the thing being catenated 
    W =[None]*(2*NumofRows)
    for n in range(NumofRows): #Alternating Catenation 
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

