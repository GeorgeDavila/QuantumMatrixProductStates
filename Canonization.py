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
    
def MpsStructure(Num): 
    return [[None]*2]*int(Num) #2 for the binary structure of the M's
    
def dotprod(A,B): 
    if len(A[0]) != len(B):
        print "WARNING: dot product off here, "
    else:
        return numpy.dot( A , B )
#----------------End of SubRoutines---------------------------------------------------------------------------------

'''Given some MPS \ket{\Psi} = \sum_{x}  M_N^{x_N} ... M_1^{x_1} \ket{x}  '''

#N is number of qubits basically, and have 2 M's (M^0, M^1) for each N (since we split it in half into bits x_i = 0,1 for each i) i.e. for spin chain of 10 electrons have 2 states up/down for each

N=6 
#Structure= [None]*int(N)
M= MpsStructure(N)
Mtilde= MpsStructure(N)
S= MpsStructure(N)
C= MpsStructure(N)

print "This code indexes from 0 to N-1  (unlike the notes which go from 0 to N), but clearly same overall process"
#Example MPS (all 1's)
M[5][0] =  numpy.ones( (1,2) )   # = M[N-1]
M[4][0] =  numpy.ones( (2,4) ) 
M[3][0] =  numpy.ones( (4,8) ) 
M[2][0] =  numpy.ones( (8,4) ) 
M[1][0] =  numpy.ones( (4,2) ) 
M[0][0] =  numpy.ones( (2,1) )  

M[5][1] =  numpy.ones( (1,2) )  # = M[N-1]
M[4][1] =  numpy.ones( (2,4) ) 
M[3][1] =  numpy.ones( (4,8) ) 
M[2][1] =  numpy.ones( (8,4) ) 
M[1][1] =  numpy.ones( (4,2) ) 
M[0][1] =  numpy.ones( (2,1) )  

'''M[5][0] = numpy.ndarray.tolist( numpy.ones( (1,2) ) )  # = M[N-1]
M[4][0] = numpy.ndarray.tolist( numpy.ones( (2,4) ) )
M[3][0] = numpy.ndarray.tolist( numpy.ones( (4,8) ) )
M[2][0] = numpy.ndarray.tolist( numpy.ones( (8,4) ) )
M[1][0] = numpy.ndarray.tolist( numpy.ones( (4,2) ) )
M[0][0] = numpy.ndarray.tolist( numpy.ones( (2,1) )  ) 

M[5][1] = numpy.ndarray.tolist( numpy.ones( (1,2) ) )  # = M[N-1]
M[4][1] = numpy.ndarray.tolist( numpy.ones( (2,4) ) )
M[3][1] = numpy.ndarray.tolist( numpy.ones( (4,8) ) )
M[2][1] = numpy.ndarray.tolist( numpy.ones( (8,4) ) )
M[1][1] = numpy.ndarray.tolist( numpy.ones( (4,2) ) )
M[0][1] = numpy.ndarray.tolist( numpy.ones( (2,1) )  )''' 


#=================Orthonormalization=======================================================================================
Z = Structure(N+1)
Z[0] = 1
#-----defns---------
B = MpsStructure(N)
A = Structure(N)
delta = Structure(N)
X = Structure(N+1)
Y = Structure(N)
Q = Structure(N)
R = Structure(N)
lam = Structure(N)
#-----defns---------

#-----------Rightmost----------------------------

B[0][0] = M[0]
B[1][0] = M[0]
print "B[",0,"] is ", B[0]

delta[0], X[1], Y[0] = svdLambda(B[0])  # B ---svd---> X delta Y
    
print "delta[",0,"] is ", delta[0]
    
lam[0], Q[0], R[0] = diagonalization( delta[0] )   # delta ---plain decomposition---> Q lam R
print "lam[",0,"] is ", lam[0]
    
Mtilde[0] = numpy.dot( R[0] , Y[0] )
print "Mtilde[",0,"] is ", Mtilde[0]

Z[1] = tripledot( X[1], Q[0], lam[0] )
print Z[1]

#-----------Rest---------------------------------
for j in range(1,N): # range(1,N) --> [1,...,N-1]
    #M[j] = [y for x in M[j] for y in x]
    #print "M[",j,"] is ", M[j] 
    #print "Z[",j,"] is ", Z[j]
    B[j][0] = (numpy.outer( M[j][0] , Z[j] ))
    B[j][1] = (numpy.outer( M[j][1] , Z[j] ))
    
    print "B[",j,"] is ", B[j]

    delta[j], X[j+1], Y[j] = svdLambda(B[j])  # B ---svd---> X delta Y
    print svdLambda(B[j])
    print "delta[",j,"] is ", delta[j]
    
    lam[j], Q[j], R[j] = svdLambda( delta[j] )   # delta ---plain decomposition---> Q lam R
    #print "lam[",j,"] is ", lam[j]
    #print len(R[j][0])
    #print len(Y[j])
    #print R[j]
    #Y[j] = [y for x in Y[j] for y in x] #flattens Y 
    #print Y[j]
    Mtilde[j] = numpy.dot( R[j] , Y[j] )
    print "Mtilde[",j,"] is ", Mtilde[j]

    Z[j+1] = tripledot( X[j+1], Q[j], lam[j] )
    #Z[j+1] = Z[j+1][0]
    print Z[j+1]


#Need to do whole orthonormalization before beginning canonization


#=================Canonical Decomposition=======================================================================================
Lambda = Structure(N+1)
Lambda[N] = 1 #Should also have Lambda[-1 (technically)] =1, works as a good check

#-----defns---------
#S = Structure(N)
mu = Structure(N)
Delta = Structure(N)
U = Structure(N)
#Uinv = Structure(N+1)
Uinv= 1 # = Uinv[-1]
P = Structure(N)
#Pinv = Structure(N+1)
Pinv = 1 # = Pinv[-1]
#-----defns---------

# !!! Dont really need to index intermediate terms that are just gonna get absorbed anyway !!!
for k in range(N-1,-1,-1):  #[N-1, ..., 0] technically dont get the trivial identity Lambda here
    S[k] = tripledot( Pinv, Uinv, Mtilde[k] )
    print "S[",k,"] is: ", S[k]
    mu = tripledot(numpy.transpose(S[k]), Lambda[k+1], S[k]) #Shifted Lambdas 1 over from what it is in the notes so doesnt go below Lambda[0], thats why this 'INPUT' Lambda is indexed k+1
    print "mu is: ", mu
    Delta, U, Uinv = svdLambda( mu )  # Removed some indices, originally was:  Delta[k-1], U[k], Uinv[k] = diagonalization( mu[k] )
    print "Delta is: ", Delta
    Lambda[k], P, Pinv = diagonalization( Delta ) #Shifted Lambdas 1 over from what it is in the notes so doesnt go below Lambda[0], thats why this 'OUTPUT' Lambda is indexed k

    C[k] = tripledot(S[k], U, P)
    
    #print "C[",k,"] is: ", C[k]
    print "Lambda[",k,"] is: ", Lambda[k]




