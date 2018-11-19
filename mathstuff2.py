import math
import numpy 


x = numpy.linspace(-2.5, 2.5, 6)

numpy.piecewise(x, [x < 0, x >= 0], [-1, 1])

y = numpy.linspace(0.0, 7035.0, 7036.0)

#cdatavec = numpy.piecewise(y, [int((y%4)) == 0, int((y%4) == 1), int((y%4) == 2), int((y%4) == 3)], [int(list(numpy.loadtxt('PCV1_1genesample.txt'))[y]), int(list(numpy.loadtxt('PCV1_2genesample.txt'))[y]), int(list(numpy.loadtxt('PCV1_3genesample.txt'))[y]), int(list(numpy.loadtxt('PCV1_4genesample.txt'))[y])])
A = numpy.array(range(0,10))

#cdatavec = numpy.piecewise(y, [(y%2) == 0.0 ,(y%2) == 1.0], [0, map(A[int(y)], y)] )

#print map(int(y), y)
#print cdatavec

#print 4%4 #4%4 == 4 mod 4 = 0 in python 
#print 15%4 #15%4 == 15 mod 4 in python 

#print list(numpy.loadtxt('PCV1_1genesample.txt'))[5-1] #Remember A[4] gives the FIFTH 5th element of list A, since in python A[0] is the first element of the list 

#print math.ceil(math.log( 3 ) / float(math.log( 2 )))
#print range(1,18)

for i in range(1,20):
    print "for", i, "the nearest above power of 2 is given by 2^", math.ceil(math.log( i ) / float(math.log( 2 ))), "=", 2**(math.ceil(math.log( i ) / float(math.log( 2 ))))
    
    
print range(16)[0:(16/2)]
print range(16)[(16/2):16]  
print len(range(16))

print 32%5

#---------------START of code used in construction of density matrix---------------------------------------------------------------------------------------------------------------------------
b = 4
NumofRows = int((2**(b/2)))

fdatavec = range(2**(b))


print NumofRows

W=[None]*NumofRows

for n in range(NumofRows): 
    W[n] = fdatavec[(NumofRows*(n)):(NumofRows*(n+1))]
    
#print W
#Note that this code takes a list of length 2**b and turns it into a 2**(b/2) X 2**(b/2) square matrix 

#---------------End of code used in construction of density matrix---------------------------------------------------------------------------------------------------------------------------

print sum(fdatavec[0:7])

#print W
#print (numpy.transpose(W))
print W[0][2]
print (numpy.transpose(W))[0][2]  #Gives correct value, so this fn should work for construction of density matrix 

#rho = [None]*(NumofRows**3)

#rho = [[None]*(NumofRows**3)]*(NumofRows**3)

#for i in range(NumofRows):
#    # iterate through rows of W
#    for j in range(NumofRows):
#        # iterate through rows of the transpose of W
#        for k in range(NumofRows):
#            rho[i][j] = (W[i][k])*((numpy.transpose(W))[j][k])

'''rho = [[0]*(math.sqrt(NumofRows)**3)]*(math.sqrt(NumofRows)**3)

for i in range(NumofRows):
    # iterate through rows of W
    for j in range(NumofRows):
        # iterate through rows of the transpose of W
        for k in range(NumofRows):
            rho[i][j] = (W[i][k])*((numpy.transpose(W))[j][k])'''

'''for xB in range(NumofRows):
    rho[xA][yA] = (W[xA][xB])*((numpy.transpose(W))[yA][xB])

for i in range(NumofRows):
    # iterate through rows of W
        for k in range(NumofRows):
            rho[i] = (W[i][k])*((numpy.transpose(W))[0][k])'''

rho = [[None]*(math.sqrt(NumofRows)**3)]*(math.sqrt(NumofRows)**3)

for i in range(NumofRows):
    # iterate through rows of W
    for j in range(NumofRows):
        # iterate through rows of the transpose of W
        for k in range(NumofRows):
            rho[i][j] = int((W[i][k]))*int(((numpy.ndarray.tolist(numpy.transpose(W)))[j][k]))
            
            
print rho