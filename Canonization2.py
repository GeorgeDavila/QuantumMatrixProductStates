import math
import numpy 
from numpy.linalg import svd


#===================Diagonalization function================================ !! For Real matrices, note that we work with inverse in this subroutine
def diagonalization(A): 
    #Gives the diagonalized form of a diagonalizable matrix AS ARRAY
    #Also Gives unitaries P, Pdag, such that  A = P x diagonal matrix x Pdag = P x diagonalization(A) x Pdag
    eigenVals, P = numpy.linalg.eig(A)
    Pinv = numpy.linalg.inv(P) #Use just inverse since we do the diagonalization with real matrices and construct unitaries from eigenvectors
    
    diagonalizedForm=  numpy.zeros( (len(eigenVals), len(eigenVals)) )
    for i in range(len(eigenVals)):
        diagonalizedForm[i][i]= (eigenVals[i]) 
        
    return diagonalizedForm, P, Pinv  #returns arrays
#===================Diagonalization function================================

#===================Construction of Lambda matrix subroutine================================
def svdLambda(B):   #returns the tuple: Lambda, svdU, svdV   #Gives the properly-shaped (i.e. can be multiplied by U, V) singular value matrix Lambda AS ARRAY
    svdU, singularValues, svdV = svd(B) 
    
    Lambda= numpy.zeros( ( len(svdU), len(svdV) ) )
    for i in range(len(singularValues)):
        Lambda[i][i]= (singularValues[i]) 
        
    return Lambda, svdU, svdV #returns an array
#===================Construction of Lambda matrix subroutine================================

def tripledot(A,B,C): #does dot product A.B.C in that order 
    return numpy.dot( numpy.dot( A , B) , C )
    
def Structure(Num): 
    return [[None]]*int(Num)
    

#----------------End of SubRoutines---------------------------------------------------------------------------------

'''Given some MPS \ket{\Psi} = \sum_{x}  M_N^{x_N} ... M_1^{x_1} \ket{x}  '''

#N is number of qubits basically, and have 2 M's (M^0, M^1) for each N (since we split it in half into bits x_i = 0,1 for each i) i.e. for spin chain of 10 electrons have 2 states up/down for each

N=6 
#Structure= [None]*int(N)
M= Structure(N)
Mtilde= Structure(N)
C= Structure(N)

#print "This code indexes from 0 to N-1  (unlike the notes which go from 0 to N), but clearly same overall process"
#Example MPS (all 1's)
M[5] = numpy.ndarray.tolist( numpy.ones( (1,2) ) )  # = M[N-1]
M[4] = numpy.ndarray.tolist( numpy.ones( (2,4) ) )
M[3] = numpy.ndarray.tolist( numpy.ones( (4,8) ) )
M[2] = numpy.ndarray.tolist( numpy.ones( (8,4) ) )
M[1] = numpy.ndarray.tolist( numpy.ones( (4,2) ) )
M[0] = numpy.ndarray.tolist( numpy.ones( (2,1) )  ) 



#=================Orthonormalization=======================================================================================
Z = Structure(N+1)
Z[0] = 1
#-----defns---------
B = Structure(N)
A = Structure(N)
delta = Structure(N)
X = Structure(N+1)
Y = Structure(N)
Q = Structure(N)
R = Structure(N)
lam = Structure(N)
#-----defns---------

'''#-----------Rightmost----------------------------
B[0][0] = numpy.outer( M[0] , Z[0] )
B[0][1] = numpy.outer( M[0] , Z[0] )
print "B[",0,"] is ", B[0]

delta[0], X[0+1], Y[0] = svdLambda(B[0])  # B ---svd---> X delta Y
    #print svdLambda(B[0])
    #print "delta[",0,"] is ", delta[0]
    
lam[0], Q[0], R[0] = svdLambda( delta[0] )   # delta ---plain decomposition---> Q lam R
    #print "lam[",0,"] is ", lam[0]
    #print R[0]
    #print Y[0]
Mtilde[0] = numpy.dot( R[0] , Y[0] )
    #print "Mtilde[",0,"] is ", Mtilde[0]

Z[0+1] = tripledot( X[0+1], Q[0], lam[0] )'''



#-----------Rest---------------------------------
for j in range(N): # range(1,N) --> [1,...,N-1]
    B[j] = numpy.outer( M[j] , Z[j])
    
    #print "B[",j,"] is ", B[j]

    delta[j], X[j+1], Y[j] = svdLambda(B[j])  # B ---svd---> X delta Y
    #print svdLambda(B[j])
    #print "delta[",j,"] is ", delta[j]
    
    lam[j], Q[j], R[j] = svdLambda( delta[j] )   # delta ---plain decomposition---> Q lam R
    #print "lam[",j,"] is ", lam[j]
    #print R[j]
    #print Y[j]
    Mtilde[j] = numpy.dot( R[j] , Y[j] )
    #print "Mtilde[",j,"] is ", Mtilde[j]

    Z[j+1] = tripledot( X[j+1], Q[j], lam[j] )
    if j>0: 
        Z[j+1] = Z[j+1][
    else:
        Z[j+1] = Z[j+1]
    print Z[j+1]


#Need to do whole orthonormalization before beginning canonization


#=================Canonical Decomposition=======================================================================================
Lambda = Structure(N+1)
Lambda[N] = 1 #Should also have Lambda[0] =1, works as a good check

#-----defns---------
S = Structure(N)
mu = Structure(N)
Delta = Structure(N)
U = Structure(N)
Uinv = Structure(N+1)
Uinv[N] = 1 # = Uinv[-1]
P = Structure(N)
Pinv = Structure(N+1)
Pinv[N] = 1 # = Pinv[-1]
#-----defns---------

# !!! Dont really need to index intermediate terms that are just gonna get absorbed anyway !!!
for k in range(N-1,-1,-1):  #[N-1, ..., 0] technically dont get the trivial identity Lambda here
    S[k] = tripledot( Pinv[k+1], Uinv[k+1], Mtilde[k] )

    mu[k] = tripledot(numpy.transpose(S[k]), Lambda[k+1], S[k]) #Shifted Lambdas 1 over from what it is in the notes so doesnt go below Lambda[0], thats why this 'INPUT' Lambda is indexed k+1

    Delta[k-1], U[k], Uinv[k] = diagonalization( mu[k] )  # Removed some indices, originally was:  Delta[k-1], U[k], Uinv[k] = diagonalization( mu[k] )

    Lambda[k], P[k], Pinv[k] = diagonalization( Delta[k-1] ) #Shifted Lambdas 1 over from what it is in the notes so doesnt go below Lambda[0], thats why this 'OUTPUT' Lambda is indexed k

    C[k] = tripledot(S[k], U[k], P[k])
    
    print "Lambda[",k,"] is: ", Lambda[k]




