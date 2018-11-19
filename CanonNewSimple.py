import math
import numpy 
from numpy.linalg import svd


#===================Diagonalization function================================ 
# !! For Real matrices, note that we work with inverse in this subroutine # i.e want unitaries (their own inverse) but deal with reals so effectively want inverses
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

#===================Construction of SQUARE Lambda matrix subroutine================================
def SQUAREsvdLambda(B):   #returns the tuple: Lambda, svdU, svdV   #Gives the properly-shaped (i.e. can be multiplied by U, V) singular value matrix Lambda AS ARRAY
    svdU, singularValues, svdV = svd(B) 
    
    Lambda= numpy.zeros( ( len(svdU), len(svdU) ) )
    for i in range(len(singularValues)):
        Lambda[i][i]= (singularValues[i]) 
        
    return Lambda, svdU, svdV #returns an array
#===================Construction of SQUARE Lambda matrix subroutine================================

def tripledot(A,B,C): #does dot product A.B.C in that order 
    return numpy.dot( numpy.dot( A , B) , C )
    
def Structure(Num): 
    return [[None]]*int(Num)
    
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
M= Structure(N)
M2= Structure(N)
Mtilde= Structure(N)
C= Structure(N)
rold=1


print "This code indexes from 0 to N-1  (unlike the notes which go from 0 to N), but clearly same overall process"
#Example MPS (all 1's)
'''M[5] = numpy.ndarray.tolist( numpy.ones( (1,2) ) )  # = M[N-1]
M[4] = numpy.ndarray.tolist( numpy.ones( (2,4) ) )
M[3] = numpy.ndarray.tolist( numpy.ones( (4,8) ) )
M[2] = numpy.ndarray.tolist( numpy.ones( (8,4) ) )
M[1] = numpy.ndarray.tolist( numpy.ones( (4,2) ) )
M[0] = numpy.ndarray.tolist( numpy.ones( (2,1) )  ) '''


#Random M matrices with values from 0 to 100:
M[5] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (1,2) ) )  # = M[N-1]
M[4] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (2,4) ) )
M[3] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (4,8) ) )
M[2] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (8,4) ) )
M[1] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (4,2) ) )
M[0] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (2,1) )  ) 

M2[5] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (1,2) ) )  # = M[N-1]
M2[4] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (2,4) ) )
M2[3] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (4,8) ) )
M2[2] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (8,4) ) )
M2[1] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (4,2) ) )
M2[0] = numpy.ndarray.tolist( 10000 * numpy.random.random_sample( (2,1) )  ) 

#=================Orthonormalization=======================================================================================
Z = Structure(N+1)
Z[0] = 1

#-----defns---------
B = Structure(N)
#-----defns---------

#-----------Rightmost----------------------------
B[0] = numpy.concatenate((M[0],M2[0]) , axis=1 )
print "B[",0,"] is ", B[0]

delta, u, A = svdLambda(B[0])  # B ---svd---> u delta A
print "delta[",0,"] is ", delta

#-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------
#For rightmost case just want to set lam[0]= delta[0]
diagdelta = delta
lam = numpy.diag(numpy.diag(diagdelta))
    
if len(diagdelta) == max( len(diagdelta), len(diagdelta[0]) ):
    Q = numpy.eye( len(diagdelta) , len(diagdelta[0]) )
else:
    Q = numpy.eye( len(diagdelta), len(diagdelta) )
    
if len(diagdelta[0]) == max( len(diagdelta), len(diagdelta[0]) ):
    p = numpy.eye( len(diagdelta), len(diagdelta[0]) )
else:
    p = numpy.eye( len(diagdelta[0]), len(diagdelta[0]) )
#-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------

print "lam[",0,"] is ", lam
print "Q[",0,"] is ", Q
print "u[",1,"] is ", u
print "A[",0,"] is ", A

Mtilde[0] = numpy.dot( p , A )    #numpy.tensordot( p , A, axes=2 )
print "Mtilde[",0,"] is ", Mtilde[0]
print "Mtilde[",0,"][0] is ", Mtilde[0][0]
Z = tripledot( u, Q, lam )
print "Z[",1,"] is ", Z

#-----------Rest---------------------------------
for j in range(1,N): # range(1,N) --> [1,...,N-1]
    B[j] = numpy.ndarray.tolist( numpy.dot( M[j] , Z ) )#numpy.concatenate(( numpy.dot( M[j] , Z ) , numpy.dot( M[j] , Z ) ) , axis=1 ) #numpy.dot( M[j] , Z )
    print "B[",j,"] is ", B[j]
    print svdLambda(B[j])[2]
    delta, u, A = svdLambda( (B[j]) )  # B ---svd---> u delta A

    print "delta[",j,"] is ", delta
    #-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------
    diagdelta = delta
    lam = numpy.diag(numpy.diag(diagdelta))
    
    if len(diagdelta) == max( len(diagdelta), len(diagdelta[0]) ):
        Q = numpy.eye( len(diagdelta) , len(diagdelta[0]) )
    else:
        Q = numpy.eye( len(diagdelta), len(diagdelta) )
    
    if len(diagdelta[0]) == max( len(diagdelta), len(diagdelta[0]) ):
        p = numpy.eye( len(diagdelta), len(diagdelta[0]) )
    else:
        p = numpy.eye( len(diagdelta[0]), len(diagdelta[0]) )
    #-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------
    print "lam[",j,"] is ", lam
    
    Mtilde[j] = numpy.dot( p , A )    #numpy.tensordot( p , A, axes=2 )
    print "Mtilde[",j,"] is ", Mtilde[j]

    Z = tripledot( u, Q, lam )
    
    print "Z[",j+1,"] = ", Z
    print "len( B[",j,"] ) = ", len(B[j])
    print "len( Mtilde[",j,"] ) = ", len(Mtilde[j])
    print "Mtilde[",j,"][0] is ", Mtilde[j][0]
    print "len( Z[",j+1,"] ) = ", len(Z)


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
    #print " pre matrices-------------------------------------------"
    #print Pinv[k+1]
    #print Uinv[k+1]
    print "Mtilde[",k,"] is ", Mtilde[k] 
    #print" start Canon iteration-------------------------------------------"
    
    S[k] = tripledot( Pinv[k+1], Uinv[k+1], Mtilde[k] )
    #print "S[",k,"] is ", S[k]
    
    mu[k] = tripledot(numpy.transpose(S[k]), Lambda[k+1], S[k]) #Shifted Lambdas 1 over from what it is in the notes so doesnt go below Lambda[0], thats why this 'INPUT' Lambda is indexed k+1
    print "mu[",k,"] is ", mu[k]
    
    Delta[k-1], U[k], Uinv[k] = diagonalization( mu[k] )  # Removed some indices, originally was:  Delta[k-1], U[k], Uinv[k] = diagonalization( mu[k] )
    print "Delta[",k-1,"] is: ",  Delta[k-1]
    print numpy.dot(U[k],Uinv[k]), "<----Should be identity matrix, if not have an error "
    #-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------
    diagdelta = Delta[k-1]
    Lambda[k] = numpy.diag( numpy.diag( Delta[k-1] ) )
    
    if len(diagdelta) == max( len(diagdelta), len(diagdelta[0]) ):
        P[k] = numpy.eye( len(diagdelta) , len(diagdelta[0]) )
    else:
        P[k] = numpy.eye( len(diagdelta), len(diagdelta) )
    
    if len(diagdelta[0]) == max( len(diagdelta), len(diagdelta[0]) ):
        Pinv[k] = numpy.eye( len(diagdelta), len(diagdelta[0]) )
    else:
        Pinv[k]  = numpy.eye( len(diagdelta[0]), len(diagdelta[0]) )
    # Gives terms such that     Delta[k-1] = P x Lambda[k] x Pdag 
    #-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------
    #Lambda[k], P[k], Pinv[k] = diagonalization( Delta[k-1] ) #Shifted Lambdas 1 over from what it is in the notes so doesnt go below Lambda[0], thats why this 'OUTPUT' Lambda is indexed k
    
    print "Lambda[",k,"] is: ", Lambda[k]
    #print "S[",k,"] is ", S[k]
    #print "U[",k,"] is ", U[k]
    print "U[",k,"] is ", U[k]
    print "Uinv[",k,"] is ", Uinv[k]
    print "P[",k,"] is ", P[k]
    print "Pinv[",k,"] is ", Pinv[k]
    
    if len(P[k]) == len((P[k])[0]) and len(P[k]) == 1:
        C[k] = numpy.dot( S[k], U[k] )
    else:
        C[k] = tripledot(S[k], U[k], P[k])
    
    #print "C[",k,"] is: ", C[k]
    print "Lambda[",k,"] is: ", Lambda[k]




