import math
import numpy
from numpy.linalg import svd
import numpy.random

#Try to get this in some pdf format too 

#--------------------Start of define datavec-------------------------------------------------------------------------------------------------------------------------------------

NumSampleMin = 1
NumSampleMax = 9
ComparisonGeneLength = 1759 #1759 for PCV1 comparisons, #781 for human H1A
#If you want no filler specify that by changing x to some even power of 2 
#Here we include filler bc want to compare our dna samples to random Marchenko-Pastur distribution with same amount of filler, i.e. random but same amount of correlation with filler
#So if the DNA samples end up being random except for the filler part, they will look like these distributions
Ngenes = (NumSampleMax + 1) - NumSampleMin

x = Ngenes*ComparisonGeneLength
datavec = numpy.ndarray.tolist( numpy.random.randint( 4, size=x ) )


print 'This is a semi Marchenko-Pastur distribution from a random sequence of length', len(datavec), ' in order to compare it to an actual sample of the same length'
print 'Here we include filler bc want to compare our dna samples to random Marchenko-Pastur distribution with same amount of filler, i.e. random but same amount of correlation with filler. So if the DNA samples end up being random except for the filler part, they will look like these distributions'

#--------------------End of define datavec-------------------------------------------------------------------------------------------------------------------------------------


b = math.ceil(math.log( x ) / float(math.log( 2 ))) # this defines b so that 2^b > x, but 2^(b-1) < x, i.e. b is the nearest power of 2 above x (or equal if x is a power of 2)

print "b is ", b, "so the length of the combined and filled data sample is ", (2**b)
print "This sample has ", ( (2**b) - 1759*Ngenes ), " or ", (( (2**b) - 1759*Ngenes )/( (2**b) ) )*100, "percent filler terms"

filler = list((0*numpy.arange(((2**b) - x)) + 4))

if (2**b) == x:
    fdatavec = datavec
elif (2**b) > x:
    fdatavec = datavec + filler
else:
    print "Something went wrong. Probably didn't define the data vector as needed. Try again"
#Now fdatavec will be a filled vector of length 2^b

#Need to do W0 first and separately since the way we make W0 from fdatavec is different from how we make W1,...,Wb from the sqrt(Lambda).V's 

W0=[None]*2
W0[0] = fdatavec[0:(len(fdatavec)/2)]
W0[1] = fdatavec[(len(fdatavec)/2):(len(fdatavec))]
#Note W0 is constructed differently than subsequent W's 
#print W0

#-------------------------------------Start of SVD of W0-----------------------------------------------------------------------------------------------------------------------
# Matrix Multiplication Algorithm to isolate sqrt(Lambda).V of W = U.Lambda.V :
U=svd(W0)[0]
Lambda=svd(W0)[1]
V=svd(W0)[2]
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

#Alt threshold attempt from # https://stackoverflow.com/questions/36719997/threshold-in-2d-numpy-array 
#threshresult = (result < 25) * result
#Alt threshold attempt

M = result
#-------------------------------------End of SVD of W0-----------------------------------------------------------------------------------------------------------------------


#NOTE: Choosing to construct W0 as we do the other W'sgives rel good results, using for B in range(int(b),0,-1), i.e. using int(b) rather than int(b)-1 and defining M as M = [fdatavec]
for B in range(int(b)-1,0,-1):   #NOTE: Start at int(b)-1 since we did defined W0 previously, so this one should start by making W1
    NumofRows = int((2**(b-B))) #Number of rows of the thing being catenated 
    W =[None]*(2*NumofRows)
    
    #    for n in range(NumofRows):    #Shifted row-wise stacked Catenation
    #        W[(2*n)] = [M[n][i] for i in range(0, int(2**(B-1)))]
    #        W[((2*n)+1)] = [M[n][i] for i in range(int(2**(B-1)), int((2**B)))]

    for n in range(NumofRows): #Alternating Catenation 
        W[(2*n)] = M[n][::2] #Should it be an alternation every 2nd pt or should we do W1[0] = A[0][::Ngenes] ???????
        W[((2*n)+1)] = M[n][1::2]
    
    # Matrix Multiplication Algorithm to isolate sqrt(Lambda).V of W = U.Lambda.V :
    U=svd(W)[0]
    #------------------------------------------------------Singular value Threshold ------------------------------------------------------------------------------------------
    #athresh = svd(W)[1] #Call it athresh since a is defined elsewhere
    #athresh[ abs(athresh) < 0.00001 ] = 0 #Sets a threshold so that any value with abs(a) < 0.0000000001 just goes to 0 
    #------------------------------------------------------Singular value Threshold ------------------------------------------------------------------------------------------
    Lambda = svd(W)[1] #athresh 
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
    print Lambda #The list of singular values #Comment this line out to get it to just tell you the number of singular values and nothing more
    print "_____________________________________________________________________" #Just to separate the results a bit 
    
    #----------------------------Start of Create a .txt file of the eigenvalues-------------------------------------------------------
    list = numpy.ndarray.tolist(Lambda)
    stringRepr = '{' + str(list).translate(None, '[]') + '}' # string.translate fn should remove the [ , and ] from the string    
    SVDtxtfilename = 'RandomMPAnalogOfC' + str(NumSampleMin) + str(NumSampleMax) + 'SingValsofSVDnum' + str( int(b-B)+1 ) + '.txt'
    f = open(SVDtxtfilename,"w") #opens file with name of "SVDtxtfilename"
    f.write(stringRepr)
    f.close()
    #----------------------------End of Create a .txt file of the eigenvalues-------------------------------------------------------

