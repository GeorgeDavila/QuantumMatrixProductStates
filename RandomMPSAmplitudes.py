import math
import numpy 
from numpy.linalg import svd

#Try to get this in some pdf format too 
#Want to perfrom SVD using as many samples as possible (have 38 we can use)
#2^14 = 16384 -> math.floor(16384 / 38) = 431, so use first 431 elements of each sample
#Cutting in this way doesn't affect overall compressibilty too much since we showed that its very internally random (i.e. each individual gene sample is random)

#--------------------Start of define fdatavec-------------------------------------------------------------------------------------------------------------------------------------

NumSampleMin = 1
NumSampleMax = 4

Ngenes = (NumSampleMax + 1) - NumSampleMin
ComparisonGeneLength = 4

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
fdatavec = numpy.ndarray.tolist( numpy.random.randint( 4, size=x ) ) 
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
#Note W0 is constructed differently than subsequent W's 
#print W0

mpsM = [None]*int(b)

#-------------------------------------Start of SVD of W0-----------------------------------------------------------------------------------------------------------------------
# Matrix Multiplication Algorithm to isolate sqrt(Lambda).V of W = U.Lambda.V :
U=svd(W0)[0]   #(Commented out if not needed right now)
Lambda=svd(W0)[1]
V=svd(W0)[2]

#Took out 'result', and L and r,  stuff for Canonical Form Calculation, since we perform successive SVDs of V alone, so set M=V

GAMMA[0] = U
# LAM[0] = numpy.diag(Lambda) makes LAMBDA[i] a square matrix, which doesn't always work well 

LL= [None]*len(Lambda)
for i in range(len(Lambda)):
    LL[i] = [0]*len(V)

for n in range(len(Lambda)):
        LL[n][n] = Lambda[n]
        
LAM[0] = LL

BondDimension[0] = numpy.count_nonzero(Lambda)
thresholdvalue =  1e-300#(1e-09)

aaaLam = numpy.zeros((len(U), len(V)), float)
#aaaLam = numpy.zeros((len(U), len(V)), int) #<-----This rounds to nearest int!!!!!!!!
#numpy.fill_diagonal(aaaLam, 1)
#print aaaLam

for i in range(len(Lambda)):
        if Lambda[i] < abs(Lambda[0]) * ( thresholdvalue ) : #Sets a threshold relative to the first eigenvalue
            Lambda[i] = 0.0
            aaaLam[i] = 0.0 #Sends corresponding 1's in aaaLam to 0 
        else:
            Lambda[i] = Lambda[i]
            aaaLam[i] = Lambda[i]

M = numpy.dot(aaaLam,V)
mpsM[0] = M

#-------------------------------------End of SVD of W0-----------------------------------------------------------------------------------------------------------------------


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
        W[(2*n)] = M[n][::2] 
        W[((2*n)+1)] = M[n][1::2]
    
    #print W
    
    U=svd(W)[0]    #(Commented out if not needed right now)
    Lambda = svd(W)[1] #athresh 
    V=svd(W)[2]
    
    aaaLam = numpy.zeros((len(U), len(V)), float)
    #aaaLam = numpy.zeros((len(U), len(V)), int) #<-----This rounds to nearest int!!!!!!!!
    
    for i in range(len(Lambda)):
        if Lambda[i] < abs(Lambda[0]) * ( thresholdvalue ) : #Sets a threshold relative to the first eigenvalue
            Lambda[i] = 0.0
            aaaLam[i] = 0.0 #Sends corresponding 1's in aaaLam to 0 
        else:
            Lambda[i] = Lambda[i]
            aaaLam[i] = Lambda[i]

 
    M = numpy.dot(aaaLam,V)
    mpsM[int(b-B)] = M
    
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
    print numpy.diag(aaaLam)
    print "aaaLam has length ", len(numpy.diag(aaaLam))
    #print M
    #print "The Number of rows is ", NumofRows
    print "The length of M is (number of rows of M) ", len(M), " And the length of M[0] is (number of columns of M) ", len( (M[0]) )
    
    print "_____________________________________________________________________" #Just to separate the results a bit 
#-------------------------------------End of iterative SVD of W-----------------------------------------------------------------------------------------------------------------------

    
print "The list of Bond Dimensions (number of NONzero singular values) is ", BondDimension
print "The maximum bond dimension = ", max(BondDimension), " should be strictly less or equal to 2^ ( b/2 ) = ", 2**( int(b)/2 )
if max(BondDimension) > 2**( int(b)/2 ):
    print "WARNING: bond dimensions not properly bounded, something's wrong"
    
print "The total number of bond dimensions from the 2nd SVD onward is = ", sum(BondDimension)