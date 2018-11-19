#Here we construct a density matrix from the filled datavector and then find its eigenvalues 

import math
import numpy
from numpy.linalg import svd

datavec = list(numpy.loadtxt('PCV1_Sample_1.txt')) + list(numpy.loadtxt('PCV1_Sample_2.txt')) #+ list(numpy.loadtxt('PCV1_3genesample.txt')) + list(numpy.loadtxt('PCV1_4genesample.txt'))
#Now we want to concatenate the datavector a bit 
gene1 = list(numpy.loadtxt('PCV1_Sample_1.txt'))
gene2 = list(numpy.loadtxt('PCV1_Sample_2.txt'))
#gene3 = list(numpy.loadtxt('PCV1_3genesample.txt')) 
#gene4 = list(numpy.loadtxt('PCV1_4genesample.txt'))
#Name them here for ease but would later cut out this step 

#Ngenes = 4 
Ngenes = 2
x = len(datavec)

#Beginning of defining cdatavec <---
catdatavec = [None]*(x) #x defined as len(datavec), so defines a list of same length as the datavec (for now filled with None's)

catdatavec[::Ngenes] = gene1
catdatavec[1::Ngenes] = gene2
#catdatavec[2::4] = gene3
#catdatavec[3::4] = gene4
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

#Need to do W0 first and separately since the way we make W0 from fdatavec is different from how we make W1,...,Wb from the sqrt(Lambda).V's 

print "Just to confirm using len(fdatavec), this should match the number above ", len(fdatavec)
print "Also we only use 2 genes here in order to get an even b "
#---------------START of construction of W hat-----------------------------------------------------------------------------------------------------------------------------------------

NumofRows = int((2**(b/2)))

Wpre=[None]*NumofRows

for n in range(NumofRows): 
    Wpre[n] = fdatavec[(NumofRows*(n)):(NumofRows*(n+1))]
    
W=numpy.int_(Wpre)
#print W
#print numpy.int_(W)
#Note that this code takes a list of length 2**b and turns it into a 2**(b/2) X 2**(b/2) square matrix 

#---------------End of construction of W hat----------------------------------------------------------------------------------------------------------------------------------------------

#print [[0]*NumofRows]*NumofRows

#Compare some of these to mathematica to make sure positions are good (Although in mathematica use the W defined here)
#print W[1-1][1-1] 
#print W[3-1][1-1]
#print W[17-1][17-1]
#print W[35-1][35-1]
#print W[64-1][64-1] 
#  W[1-1][1-1] = W[0][0], just write it with the -1 to make it clear what row/column we're looking at. e.g.  W[42-1][7-1] -> 42nd row, 7th column

'''rho = [[None]*(math.sqrt(NumofRows)**3)]*(math.sqrt(NumofRows)**3)

for k in range(NumofRows):
    rho = (W[0][k])*((numpy.ndarray.tolist(numpy.transpose(W)))[0][k])'''

'''for i in range(NumofRows):
    # iterate through rows of W
    for j in range(NumofRows):
        # iterate through rows of the transpose of W
        for k in range(NumofRows):
            rho = (W[i][k])*((numpy.ndarray.tolist(numpy.transpose(W)))[j][k])'''

'''print (W[0][0:NumofRows])
print (W[0][0:NumofRows])[2] # should be 1.0
print (W[0][2]) # should be 1.0
print "----------------------------------------------------------------------------"
print ((numpy.ndarray.tolist(numpy.transpose(W)))[0][0:NumofRows])
print ((numpy.ndarray.tolist(numpy.transpose(W)))[0][0:NumofRows])[2] # should be 2.0
print ((numpy.ndarray.tolist(numpy.transpose(W)))[0][2]) # should be 2.0
print "----------------------------------------------------------------------------"
'''

#rho = [[None]*(math.sqrt(NumofRows)**3)]*(math.sqrt(NumofRows)**3)

'''for i in range(NumofRows):
    # iterate through rows of W
    for j in range(NumofRows):
        # iterate through rows of the transpose of W
        for k in range(NumofRows):
            rho = (W[i][k])*((numpy.ndarray.tolist(numpy.transpose(W)))[j][k])'''


rho = [[0]*NumofRows]*NumofRows   

for i in range(NumofRows):
    # iterate through rows of W
    for j in range(NumofRows):
        # iterate through rows of the transpose of W
        for k in range(NumofRows):
            rho[i][j] += (W[i][k]) * ((numpy.ndarray.tolist(numpy.transpose(W)))[k][j])

#result[i][j] += X[i][k] * Y[k][j]
print rho

'''rho = [[0]*NumofRows]*NumofRows   

for i in range(NumofRows):
    # iterate through rows of W
    for j in range(NumofRows):
        # iterate through rows of the transpose of W
        for k in range(NumofRows):
            rho[i][j] += (W[i][k]) * ((numpy.ndarray.tolist(numpy.transpose(W)))[k][j])

#result[i][j] += X[i][k] * Y[k][j]
print rho'''