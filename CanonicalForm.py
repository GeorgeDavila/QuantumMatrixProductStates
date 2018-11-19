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

datavec = [None]*0  #Makes it begin as an empty set, since we use a sum in a for loop, if this is not empty it will add extraneous values (e.g. using [None]*3 adds 3 None's to the final datavec 
for k in range(NumSampleMin, (NumSampleMax + 1) ):
    datavec += list(numpy.loadtxt('PCV1_Sample_' + str(k) + '.txt'))

Ngenes = (NumSampleMax + 1) - NumSampleMin
x0 = len(datavec)
b0 = math.floor(math.log( x0 ) / float(math.log( 2 )))
print b0
genelength= x0/Ngenes

cutgenelength= int(math.floor( (2**b0) / Ngenes ))
#cutgenelength= int(math.floor( 16384 / NumSampleMax )) #2^14 = 16384 = comp bound on size of sample with current comp power 

if (genelength / ( math.floor(genelength) )) != 1.0:
    print 'WARNING: genes not of same length, proceed with extreme caution or excise/add some terms to make them same. This can break/interfere with many of the other algorithms. '


print 'we are combining sample/individual number ', NumSampleMin, 'to number ', NumSampleMax
#print 'The length of this datavec, len(datavec) = ', len(datavec), ', a value which should be equal to 1759*Ngenes = ', 1759*Ngenes #Since PCV1 has 1759 base pairs 

#--------------------End of define datavec-------------------------------------------------------------------------------------------------------------------------------------

#--------------------Beginning of defining cdatavec, the catenated form-------------------------------------------------------------------------------------------------------------------------------------

catdatavec = [None]*(x0) #x defined as len(datavec), so defines a list of same length as the datavec (for now filled with None's)

#catdatavec[0::Ngenes] = gene1
#catdatavec[1::Ngenes] = gene2
#catdatavec[2::Ngenes] = gene3
#catdatavec[3::Ngenes] = gene4
#Can catenate before filling it with 4's since 4's end up being at the 

for i in range(1,Ngenes+1):               # range(NumSampleMin, (NumSampleMax + 1) ) only correct when we start at NumSampleMin = 1 
    catdatavec[(i-1)::Ngenes] = list(numpy.loadtxt('PCV1_Sample_' + str(i) + '.txt'))


del catdatavec[(cutgenelength*Ngenes) : x0+1 ] #Deletes from the (cutgenelength*Ngenes)-th part to the last term

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

M = [None]*len(Lambda)
for n in range(len(Lambda)):
    M[n] = V[n]
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
    M = [None]*len(Lambda)
    for n in range(len(Lambda)):
        M[n] = V[n]
    
    #print M
    print "For the", int(b-B)+1, "-th SVD there are", len(Lambda), "singular values and they are:" #i.e. for W_(b-B) = U.Lambda_(b-B).V, the non-zero values of Lambda_(b-B)
    print Lambda #The list of singular values #Comment this line out to get it to just tell you the number of singular values and nothing more
    #print "The Length of V is ", len(V), "The Length of V[0] is ", len( (V[0]) ), "The element V[0][0] is ", V[0][0]
    #print M
    #print V
    #print "The Number of rows is ", NumofRows
    print "The length of M is (number of rows of M) ", len(M), " And the length of M[0] is (number of columns of M) ", len( (M[0]) )
    ''' #----------------------------Start of Create a .txt file of the eigenvalues-------------------------------------------------------
    list = numpy.ndarray.tolist(V)
    stringRepr = str(list) # string.translate fn should remove the [ , and ] from the string    
    SVDtxtfilename = 'Vn' + str(NumSampleMin) + str(NumSampleMax) + 'SingValsofSVDnum' + str( int(b-B)+1 ) + '.txt'
    f = open(SVDtxtfilename,"w") #opens file with name of "SVDtxtfilename"
    f.write(stringRepr)
    f.close()
    #----------------------------End of Create a .txt file of the eigenvalues-------------------------------------------------------
    '''
    print "_____________________________________________________________________" #Just to separate the results a bit 
    
print len(GAMMA)
print len(LAM)
#print GAMMA[0]
#print LAM[0]
#print numpy.dot(GAMMA[0], LAM[0])
#print GAMMA[1]
#print LAM[1]
#print GAMMA[11]

#print GAMMA[10]
#print LAM[10]

#for i in range(len(GAMMA)): #Should be all 0's unless something broken
#    print len(GAMMA[i]) - len(GAMMA[i][0])
    
#for i in range(len(LAM)): 
#    print len(LAM[i]) - len(LAM[i][0])

#for i in range(len(GAMMA)): #Should be all 0's unless something broken
#    print len(GAMMA[i]) - len(LAM[i])


#print numpy.dot(GAMMA[9], LAM[9])

print len(numpy.dot(GAMMA[8], LAM[8])[0])
print len(numpy.dot(GAMMA[9], LAM[9]))
print len(numpy.dot(GAMMA[9], LAM[9])[0])
print len(numpy.dot(GAMMA[10], LAM[10]))
#print "Technically the Lambdas should not be square matrices, they should be 2^m x 2^n matrices n not equal to m, correct for that "
#print numpy.dot(GAMMA[11], LAM[11])
#GLmatprod=[None]*len(GAMMA)
#for i in range(len(GAMMA)):
#    GLmatprod[i] = numpy.dot(GAMMA[i], LAM[i])


#numpy.dot(GAMMA[i], LAM[i])

#print GAMMA[0]
print range(12,-1,-1)
print GAMMA[0]
#print GAMMA[1]
#print GAMMA[2]
#print numpy.dot(GAMMA[0],1)
#print LAM[7]  #Sing values after 7+1 SVDs

#print LAM[7]
#print GAMMA[11]

#print LAM[1]

#print M[7]

#f.write(stringRepr)

'''MPSconstruct = [None]*14
MPSconstruct[13] = 1

for i in range(12,-1,-1):
    MPSconstruct[i] = numpy.dot( numpy.dot(GAMMA[i], LAM[i]) , MPSconstruct[i + 1] )

print MPSconstruct[0]'''