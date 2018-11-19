import math
import numpy 
from numpy.linalg import svd

#Try to get this in some pdf format too 
#Want to perfrom SVD using as many samples as possible (have 38 we can use)
#2^14 = 16384 -> math.floor(16384 / 38) = 431, so use first 431 elements of each sample
#Cutting in this way doesn't affect overall compressibilty too much since we showed that its very internally random (i.e. each individual gene sample is random)

#--------------------Start of define fdatavec-------------------------------------------------------------------------------------------------------------------------------------

NumSampleMin = 1
NumSampleMax = 16

Ngenes = (NumSampleMax + 1) - NumSampleMin
ComparisonGeneLength = 2

b0 = math.floor(math.log( ComparisonGeneLength*Ngenes ) / float(math.log( 2 ))) #makes b0 a power s.t. a vector of length 2**b0 can fit however many cut samples of the gene we're using 
#print b0

x = 2**b0

#To do a totally random canonical MPS just set x to the length of whatever you're comparing it to
#This doesn't include filler since the canonical MPS we compare it to doesn't include any filler
#To use filler see RandomCanonFilled.py

b = math.floor(math.log( x ) / float(math.log( 2 ))) # this defines b so that 2^b > x, but 2^(b-1) < x, i.e. b is the nearest power of 2 above x (or equal if x is a power of 2)

print "b is ", b, "so the length of the combined and filled data sample is ", (2**b), "or x = ", x
print "This sample has ", ( (2**b) - x ), " or ", (( (2**b) - x )/( (2**b) ) )*100, "percent filler terms"

#fdatavec = int(x)*[0]  #NOT random, creates a list of length x filled with 0's 
fdatavec = numpy.ndarray.tolist( numpy.random.randint( 1, 2, size=x ) )# makes a list of 16 random integers with values 1, 2, 3, or 4  
#fdatavec = numpy.ndarray.tolist( numpy.random.randint( 4, size=x ) ) # makes a list of 16 random integers with values 0, 1, 2, or 3 
fdatavec = map(float, fdatavec) #converts entries of fdatavec to float just in case leaving them as ints causes the program to cutoff some stuff 
print fdatavec


#THIS IS THE TOTALLY RANDOM (unfilled) CASED (UNFILLED SO AS TO MATCH FACT THAT WE DON'T FILL DNA CANON DATAVEC)
#--------------------End of define fdatavec-------------------------------------------------------------------------------------------------------------------------------------

GAMMA = [None]*int(b)
LAM = [None]*int(b)
BondDimension = [None]*int(b)

#Need to do W0 first and separately since the way we make W0 from fdatavec is different from how we make W1,...,Wb from the sqrt(Lambda).V's (or V's in Canonical Form) 

W0=[None]*2
W0[0] = fdatavec[0:(len(fdatavec)/2)]
W0[1] = fdatavec[(len(fdatavec)/2):(len(fdatavec))]
print 'W0 is given by the following array: '
print W0
mpsM = [None]*int(b)
MPSmat = [None]*int(b)

svdW = svd(W0) #So it only performs the svd once, nit 3 times to define each
U= svdW[0]    #(Commented out if not needed right now)
Lambda = svdW[1] #athresh 
V= svdW[2]

reshapedLam = numpy.zeros((len(U), len(V)), float)
reshapedSQRTLam = numpy.zeros((len(U), len(V)), float)
    
thresholdvalue = 0.00001
for i in range(len(Lambda)):
    if Lambda[i] < abs(Lambda[0]) * ( thresholdvalue ) : #Sets a threshold relative to the first eigenvalue
        Lambda[i] = 0.0
        reshapedLam[i] = 0.0 #Sends corresponding 1's in reshapedLam to 0 
        reshapedSQRTLam[i] = 0.0 #Sends corresponding 1's in reshapedLam to 0 
    else:
        Lambda[i] = Lambda[i]
        reshapedLam[i] = Lambda[i]
        reshapedSQRTLam[i] = math.sqrt( Lambda[i] )

M0 = numpy.dot( U, reshapedSQRTLam )
M1 = numpy.dot( reshapedSQRTLam, V )

#-------------------------------------Start of iterative SVD of W-----------------------------------------------------------------------------------------------------------------------
#NOTE: Choosing to construct W0 as we do the other W'sgives rel good results, using for B in range(int(b),0,-1), i.e. using int(b) rather than int(b)-1 and defining M as M = [fdatavec]
for B in range(int(b)-1,0,-1):   #NOTE: Start at int(b)-1 since we did defined W0 previously, so this one should start by making W1
    #NumofRows = int((2**(b-B))) #Number of rows of the thing being catenated 
    NumofRows = len(Lambda)
    
    W =[None]*(2*NumofRows)
    
    #    for n in range(NumofRows):    #Shifted row-wise stacked Catenation
    #        W[(2*n)] = [M[n][i] for i in range(0, int(2**(B-1)))]
    #        W[((2*n)+1)] = [M[n][i] for i in range(int(2**(B-1)), int((2**B)))]

    for n in range(NumofRows): #Alternating Catenation 
        W[(2*n)] = W0[n][::2] 
        W[((2*n)+1)] = W0[n][1::2]
    
    print 'for b = ', b, 'W is given by the following array: '
    print W
    
    svdW = svd(W) #So it only performs the svd once, nit 3 times to define each
    U= svdW[0]    #(Commented out if not needed right now)
    Lambda = svdW[1] #athresh 
    V= svdW[2]
    
    
    print 'U is size ', len(U), ' by ', len( U[0] )
    
    reshapedLam = numpy.zeros((len(U), len(V)), float)
    reshapedSQRTLam = numpy.zeros((len(U), len(V)), float)
    
    for i in range(len(Lambda)):
        if Lambda[i] < abs(Lambda[0]) * ( thresholdvalue ) : #Sets a threshold relative to the first eigenvalue
            Lambda[i] = 0.0
            reshapedLam[i] = 0.0 #Sends corresponding 1's in reshapedLam to 0 
            reshapedSQRTLam[i] = 0.0 #Sends corresponding 1's in reshapedLam to 0 
        else:
            Lambda[i] = Lambda[i]
            reshapedLam[i] = Lambda[i]
            reshapedSQRTLam[i] = math.sqrt( Lambda[i] )

 
    M = numpy.dot(reshapedSQRTLam,V)   #numpy.outer( reshapedSQRTLam , V ) ???
    mpsM[int(b-B)] = M
    
    '---------------------->' #MPSmat = numpy.outer( U , reshapedSQRTLam )
    MPSmat[int(b-B)] = numpy.outer( U , reshapedSQRTLam )
    
    GAMMA[int(b-B)] = U     #Really GAMMA[int(b-B)] corresponds to the Gamma_n ^n for the n = int(b-B) + 1 -th SVD, e.g. GAMMA[1] = Gamma_2
    #LAM[int(b-B)] =  numpy.diag(Lambda)    #Really LAM[int(b-B)] corresponds to the Lambda_n ^n for the n = int(b-B) + 1 -th SVD, e.g. LAM[1] = Lambda_2
    LL= [None]*len(U)
    for i in range(len(U)):
        LL[i] = [0]*len(V)
    
    for n in range(len(Lambda)):
        LL[n][n] = Lambda[n]
        
    LAM[int(b-B)] = LL
    
    BondDimension[int(b-B)] = numpy.count_nonzero(Lambda)
    print "For the", int(b-B)+1, "-th SVD there are", len(Lambda), "singular values and they are:" #i.e. for W_(b-B) = U.Lambda_(b-B).V, the non-zero values of Lambda_(b-B)
    print Lambda #The list of singular values #Comment this line out to get it to just tell you the number of singular values and nothing more
    print "The number of nonzero singular values is:  ", numpy.count_nonzero(Lambda)
    #print "The Length of V is ", len(V), "The Length of V[0] is ", len( (V[0]) ), "The element V[0][0] is ", V[0][0]
    #print numpy.diag(reshapedLam)
    #print "reshapedLam has length ", len(numpy.diag(reshapedLam))
    print M
    #print "The Number of rows is ", NumofRows
    #print "The length of M is (number of rows of M) ", len(M), " And the length of M[0] is (number of columns of M) ", len( (M[0]) )
    
    print "_____________________________________________________________________" #Just to separate the results a bit 
#-------------------------------------End of iterative SVD of W-----------------------------------------------------------------------------------------------------------------------




#---------------------------Amplitude Calculation----------------------------------------------------------------------------

'''for B in range(int(b),0,-1):
    print ( int(b) - B )
    print 'MPSmat[', ( int(b) - B ), '] is given by '
    print (MPSmat[int(b-B)])
    print "________________________________________________________________________________________________________" #Just to separate the results a bit 
 '''   
#MPSampReconst = [None]*(int(b)/2) 
#MPSampReconst[0] = numpy.dot( MPSmat[0] , MPSmat[1] )


#print numpy.outer( ( numpy.outer( MPSmat[0] , MPSmat[1] ) ) , MPSmat[2] )

#print numpy.dot( MPSmat[2] , M )


#print numpy.outer( MPSmat[0] , numpy.outer( MPSmat[1] , numpy.outer( MPSmat[2] , M )))
#print numpy.flatnonzero( numpy.outer( MPSmat[0] , numpy.outer( MPSmat[1] , numpy.outer( MPSmat[2] , M ))) )

#print numpy.outer( MPSmat[0] , numpy.outer( MPSmat[1] , MPSmat[2]))

#AAmpcalc = numpy.outer( MPSmat[0] , numpy.outer( MPSmat[1] , numpy.outer( MPSmat[2] , M )))
#print 'AAmpcalc gives ', AAmpcalc

print 'Matrices are '
print MPSmat[0][0]
print MPSmat[1][0]
print MPSmat[2][0]
#print M

print MPSmat[0][0][0] * MPSmat[1][0][0] * MPSmat[2][0][0]

print numpy.dot( MPSmat[0][0] , numpy.dot( MPSmat[1][0] , MPSmat[2][0] ) )