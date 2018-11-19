#Here we construct a density matrix from the filled datavector and then find its eigenvalues 

import math
import numpy
from numpy.linalg import svd

#datavec = list(numpy.loadtxt('PCV1_1genesample.txt')) + list(numpy.loadtxt('PCV1_2genesample.txt')) #+ list(numpy.loadtxt('PCV1_3genesample.txt')) + list(numpy.loadtxt('PCV1_4genesample.txt'))
#Now we want to concatenate the datavector a bit 
#gene1 = list(numpy.loadtxt('PCV1_1genesample.txt'))
#gene2 = list(numpy.loadtxt('PCV1_2genesample.txt'))
#gene3 = list(numpy.loadtxt('PCV1_3genesample.txt')) 
#gene4 = list(numpy.loadtxt('PCV1_4genesample.txt'))
#Name them here for ease but would later cut out this step 

datavec = list(numpy.loadtxt('H1A_genesample.txt'))  #For H1A 

#Ngenes = 4
#Ngenes = 2
Ngenes = 1

x = len(datavec)

#Beginning of defining cdatavec <---
catdatavec = [None]*(x) #x defined as len(datavec), so defines a list of same length as the datavec (for now filled with None's)

catdatavec[::Ngenes] = datavec
#catdatavec[::Ngenes] = gene1
#catdatavec[1::Ngenes] = gene2
#catdatavec[2::4] = gene3
#catdatavec[3::4] = gene4
#Can catenate before filling it with 4's since 4's end up being at the 

#Candidate function so we don't need to use gene1, gene2, etc. :
#for j in range(Ngenes):
#    catdatavec[j::Ngenes] = datavec[ (j*(x/4)) :  ((j+1)*(x/4)) ]


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

#print "Just to confirm using len(fdatavec), this should match the number above ", len(fdatavec)
#print "Also we only use 2 genes here in order to get an even b "
#---------------START of construction of W hat-----------------------------------------------------------------------------------------------------------------------------------------

NumofRows = int((2**(b/2)))

W=[None]*NumofRows

for n in range(NumofRows): 
    W[n] = fdatavec[(NumofRows*(n)):(NumofRows*(n+1))]
    
print W
#Note that this code takes a list of length 2**b and turns it into a 2**(b/2) X 2**(b/2) square matrix 

#---------------End of construction of W hat----------------------------------------------------------------------------------------------------------------------------------------------

