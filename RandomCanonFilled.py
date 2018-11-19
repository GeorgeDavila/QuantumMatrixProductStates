import math
import numpy 
from numpy.linalg import svd

#Try to get this in some pdf format too 
#Want to perfrom SVD using as many samples as possible (have 38 we can use)
#2^14 = 16384 -> math.floor(16384 / 38) = 431, so use first 431 elements of each sample
#Cutting in this way doesn't affect overall compressibilty too much since we showed that its very internally random (i.e. each individual gene sample is random)

#--------------------Start of define datavec-------------------------------------------------------------------------------------------------------------------------------------

NumSampleMin = 1
NumSampleMax = 4   #WE CURRENTLY CUT the samples in this program, don't nec need to #Using 8 gives cutlength = 2048 > the actual gene length, which is why the program fills it 

#datavec = [None]*0  #Makes it begin as an empty set, since we use a sum in a for loop, if this is not empty it will add extraneous values (e.g. using [None]*3 adds 3 None's to the final datavec 
#for k in range(NumSampleMin, (NumSampleMax + 1) ):
#    datavec += list(numpy.loadtxt('PCV1_Sample_' + str(k) + '.txt'))

#nonrandomdatavec = datavec
ComparisonGeneLength = 1759

Ngenes = (NumSampleMax + 1) - NumSampleMin

x0 = Ngenes*ComparisonGeneLength
b0 = math.floor(math.log( x0 ) / float(math.log( 2 )))
print b0
genelength= x0/Ngenes

cutgenelength= int(math.floor( (2**b0) / Ngenes ))
#cutgenelength= int(math.floor( 16384 / NumSampleMax )) #2^14 = 16384 = comp bound on size of sample with current comp power 


datavec = numpy.ndarray.tolist( numpy.random.randint( 4, size=x0 ) )  #<<<<-------------------Part that makes datavec random

if (genelength / ( math.floor(genelength) )) != 1.0:
    print 'WARNING: genes not of same length, proceed with extreme caution or excise/add some terms to make them same. This can break/interfere with many of the other algorithms. '


print 'we are combining sample/individual number ', NumSampleMin, 'to number ', NumSampleMax
#print 'The length of this datavec, len(datavec) = ', len(datavec), ', a value which should be equal to 1759*Ngenes = ', 1759*Ngenes #Since PCV1 has 1759 base pairs 

#--------------------End of define datavec-------------------------------------------------------------------------------------------------------------------------------------

#--------------------Beginning of defining cdatavec, the catenated form-------------------------------------------------------------------------------------------------------------------------------------

catdatavec = datavec

#Since this is the test of Canonical decomposition on random data, we don't apply a catenation. Only define catdatavec so we don't have to change the code too much (so less margin for error there)

x = len(catdatavec)

print 'The cut gene length here is ', cutgenelength, ' so the length of the cut and catenated sample is now ', x, ', which should be equal to cutgenelength*Ngenes = ', (cutgenelength*Ngenes) #Since cut PCV1 has 431 base pairs  
#--------------------End of defining cdatavec, the catenated form-------------------------------------------------------------------------------------------------------------------------------------


b = math.ceil(math.log( x ) / float(math.log( 2 ))) # this defines b so that 2^b > x, but 2^(b-1) < x, i.e. b is the nearest power of 2 above x (or equal if x is a power of 2)

print "b is ", b, "so the length of the combined and filled data sample is ", (2**b)
print "This sample has ", ( (2**b) - x ), " or ", (( (2**b) - x )/( (2**b) ) )*100, "percent filler terms"

filler = list((0*numpy.arange(((2**b) - x)) + 4))

if (2**b) == x:
    fdatavec = catdatavec
elif (2**b) > x:
    fdatavec = catdatavec + filler
else:
    print "Something went wrong. Probably didn't define the data vector as needed. Try again"
#Now fdatavec will be a filled vector of length 2^b

#Need to do W0 first and separately since the way we make W0 from fdatavec is different from how we make W1,...,Wb from the sqrt(Lambda).V's (or V's in Canonical Form) 

W0=[None]*2
W0[0] = fdatavec[0:(len(fdatavec)/2)]
W0[1] = fdatavec[(len(fdatavec)/2):(len(fdatavec))]
#Note W0 is constructed differently than subsequent W's 
#print W0

GAMMA = [None]*int(b)
LAM = [None]*int(b)
BondDimension = [None]*int(b)

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

aaa = numpy.zeros((len(U), len(V)), int)
numpy.fill_diagonal(aaa, 1)
#print aaa

M = numpy.dot(aaa,V)
#print M
#print aaa
#-------------------------------------End of SVD of W0-----------------------------------------------------------------------------------------------------------------------


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
    
    aaa = numpy.zeros((len(U), len(V)), int)
    numpy.fill_diagonal(aaa, 1)
    #aaa is just Lambda with the nonzero values set to be 1, we use it to effectively excise the arbitrary terms in V 
    #print aaa
    
    for i in range(len(Lambda)):
        if Lambda[i] < abs(Lambda[0]) * (1e-09) : #Sets a threshold relative to the first eigenvalue
            Lambda[i] = 0.0
            aaa[i] = 0 #Sends corresponding 1's in aaa to 0 
        else:
            Lambda[i] = Lambda[i]
    
#    for i in range(len(V)):
#        for j in range(len(V[0])):
#            if V[i][j] < abs( numpy.amax(V) ) * (1e-09) :
#                V[i][j] = 0.0
#            else:
#                V[i][j] = V[i][j]
 

   

    M = numpy.dot(aaa,V)
    
    #Effectively makes M = aaa.V BUT aaa is constructed to be (I   0I   0I   0I), in block form, i.e. all L does is reshape V so that only appropriate parts of V are retained 
    
    #print V
    
    GAMMA[int(b-B)] = U     #Really GAMMA[int(b-B)] corresponds to the Gamma_n ^n for the n = int(b-B) + 1 -th SVD, e.g. GAMMA[1] = Gamma_2
    #LAM[int(b-B)] =  numpy.diag(Lambda)    #Really LAM[int(b-B)] corresponds to the Lambda_n ^n for the n = int(b-B) + 1 -th SVD, e.g. LAM[1] = Lambda_2
    LL= [None]*len(U)
    for i in range(len(U)):
        LL[i] = [0]*len(V)
    
    for n in range(len(Lambda)):
        LL[n][n] = Lambda[n]
        
    LAM[int(b-B)] = LL
    
    #V=svd(W)[2] #<<<<<<-------------------------------PUT V HERE SO THAT LAM AND LL ARE PROPERLY STRUCTURED!!!!!<<<<<<<<<<<<
    #i.e. the decomposition and reshape gives  LAM[n].V[n] --> W[n] = GAMMA[n+1].LAM[n+1].V[n+1]
    #------------------------------------- (a X b) (b X b) --> (b X a) = (a X a)   (a X b)   (b X b)
    
    #M = V #Need to redefine M so the for loop acts on result each time 
    #M = [None]*len(Lambda)
    #for n in range(len(Lambda)):
    #    M[n] = V[n]
    
    BondDimension[int(b-B)] = numpy.count_nonzero(Lambda)
    #print M
    print "For the", int(b-B)+1, "-th SVD there are", len(Lambda), "singular values and they are:" #i.e. for W_(b-B) = U.Lambda_(b-B).V, the non-zero values of Lambda_(b-B)
    print Lambda #The list of singular values #Comment this line out to get it to just tell you the number of singular values and nothing more
    print "The number of nonzero singular values is:  ", numpy.count_nonzero(Lambda)
    #print "The Length of V is ", len(V), "The Length of V[0] is ", len( (V[0]) ), "The element V[0][0] is ", V[0][0]
    print numpy.diag(aaa)
    print "aaa has length ", len(numpy.diag(aaa))
    #print M
    #print V
    #print "The Number of rows is ", NumofRows
    print "The length of M is (number of rows of M) ", len(M), " And the length of M[0] is (number of columns of M) ", len( (M[0]) )
    
    '''#----------------------------Start of Create a .txt file of the Gammas and Lambdas-------------------------------------------------------
    list = numpy.ndarray.tolist(U)
    stringRepr = str(list)   
    GAMtxtfilename = 'Gamma' + str(NumSampleMin) + str(NumSampleMax) + 'ofSVDnum' + str( int(b-B)+1 ) + '.txt'
    fgam = open(GAMtxtfilename,"w") #opens file with name of "GAMtxtfilename"
    fgam.write(stringRepr)
    fgam.close()
    
    list2 = LAM[int(b-B)]
    stringRepr2 = str(list2)   
    LAMtxtfilename = 'Lambda' + str(NumSampleMin) + str(NumSampleMax) + 'ofSVDnum' + str( int(b-B)+1 ) + '.txt'
    flam = open(LAMtxtfilename,"w") #opens file with name of "LAMtxtfilename"
    flam.write(stringRepr2)
    flam.close()
    #----------------------------End of Create a .txt file of the Gammas and Lambdas------------------------------------------------------------
    '''
    
    
    print "_____________________________________________________________________" #Just to separate the results a bit 
    

#print nonrandomdatavec
#print datavec
#print fdatavec

print "The list of Bond Dimesions (number of NONzero singular values) is ", BondDimension