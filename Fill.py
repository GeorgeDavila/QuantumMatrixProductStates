import math
import numpy

datavec = list(numpy.loadtxt('PCV1_1genesample.txt')) + list(numpy.loadtxt('PCV1_2genesample.txt')) + list(numpy.loadtxt('PCV1_3genesample.txt')) + list(numpy.loadtxt('PCV1_4genesample.txt'))
    
x = len(datavec)

b = math.ceil(math.log( x ) / float(math.log( 2 )))

filler = list((0*numpy.arange(((2**b) - x)) + 4))

print x
print b
#print filler

if (2**b) == x:
    fdatavec = datavec
elif (2**b) > x:
    fdatavec = datavec + filler
else:
    print "Something went wrong. Probably didn't define the data vector as needed. Try again"

#print fdatavec #no need to print this, just use it for the example 
#print len(fdatavec)



#Now fdatavec will be a filled vector of length 2^b