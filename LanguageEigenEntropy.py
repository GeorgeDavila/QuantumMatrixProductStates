#Here we construct a density matrix from the filled datavector and then find its eigenvalues 

import math
import numpy
from numpy.linalg import svd

#Try to get this in some pdf format too 

#--------------------Start of define datavec-------------------------------------------------------------------------------------------------------------------------------------

datavec = list(numpy.loadtxt('OdysseyBook1ButlerOctalTranslation.txt'))

x = len(datavec)

#print 'We are combining sample/individual number ', NumSampleMin, 'to number ', NumSampleMax
print 'The length of this datavec, len(datavec) = ', len(datavec)

#--------------------End of define datavec-------------------------------------------------------------------------------------------------------------------------------------


b = math.ceil(math.sqrt(x)) 

print "b is ", b, "so the length of the combined and filled data sample is b^2 = ", (b**2), " (its b^2 since this is for constructing a density matrix)"

filler = list((0*numpy.arange(((b**2) - x)) + 4))

if (b**2) == x:
    fdatavec = datavec   #Changed catdatavec to datavec, sinces there's no need for catenation here 
elif (b**2) > x:
    fdatavec = datavec + filler #Changed catdatavec to datavec, sinces there's no need for catenation here 
else:
    print "Something went wrong. Probably didn't define the data vector as needed. Try again"
#Now fdatavec will be a filled vector of length n^2

print "Then the length of the filled data vector is len(fdatavec) = ", len(fdatavec)

#---------------START of construction of W hat-----------------------------------------------------------------------------------------------------------------------------------------

NumofRows = int(b)

Wpre=[None]*NumofRows

for n in range(NumofRows): 
    Wpre[n] = fdatavec[(NumofRows*(n)):(NumofRows*(n+1))]
    
W=numpy.int_(Wpre)
#print W
#print numpy.int_(W)
#Note that this code takes a list of length x and turns it into a b^2 square matrix with b = math.ceil(math.sqrt(x)) 

#---------------End of construction of W hat----------------------------------------------------------------------------------------------------------------------------------------------


rho = numpy.dot(W,numpy.transpose(W))   #This is the density matrix

#print rho

a = numpy.linalg.eigvals(rho)

a[ abs(a) < 0.0000000001 ] = 0 #Sets a threshold so that any value with abs(a) < 0.0000000001 just goes to 0 

#if any term has nonzero imaginary part, numpy.imag(), the calc is invalid 
if sum(numpy.imag(a)) > 0.0:
    print "_________________________________________________________________________________________________________________________"
    print "There are nonzero imaginary parts!!!!!! Something went wrong. Probably didn't define the data vector as needed. Try again"
    print "_________________________________________________________________________________________________________________________"

#print sum(abs(a))
#print a

normalizedeigens = (abs(a) / sum(abs(a))) #normalizes the eigenvalues

print "There are b - len(normalizedeigens) = ", (b - len(normalizedeigens) ), " number of zero-valued eigenvalues"
print "Although if you notice that there seems to be an atypically sharp (>2 orders of magnitude, at least in the examples I've done so far) drop in the eigenvalues after the (desired) drop in the first few this might very well mean that those after this drop should actually be zero. This is a common issue in numpy and can be fixed by setting a threshold where the drop occurs. "
print "The normalized eigenvalues are: ", normalizedeigens  
print "The sum of the normalized eigenvalues: ", sum(normalizedeigens), " should be equal to 1.0 "

print "There are ", (numpy.nonzero(a)[-1][-1]) + 1," = ", len((numpy.nonzero(a)[0])), "(these should be equal) nonzero eigenvalues of the density matrix "

termsShannonEntropy =  0
for k in range(len((numpy.nonzero(a)[0]))): #Take logarithm, so it must sum over only nonzero eigenvalues 
    termsShannonEntropy +=  - ( normalizedeigens[k] * math.log( normalizedeigens[k] ) )

#ShannonEntropy = sum(termsShannonEntropy)

print "The 0-th Renyi Entropy (Hartley Entropy) is ", math.log(len((numpy.nonzero(a)[0])))
print "The 1-th Renyi Entropy (Shannon Entropy) is ", termsShannonEntropy

for alpha in range(2, 15 + 1):
    print "The ", alpha, "-th Renyi Entropy is ", (((1-alpha)**(-1)) * (math.log( sum( ((normalizedeigens)**alpha) ) ) / float(math.log( 2 ))))
    
