import math
import numpy 
from numpy.linalg import svd


#===================Diagonalization function================================ 
# !! For Real matrices, note that we work with inverse in this subroutine # i.e want unitaries (their own inverse) but deal with reals so effectively want inverses
'''NOTE: Not yet made to select only non-zero terms of diag'''
def diagMAT(A,rold): 
    #Gives a square matrix of the diagonal terms of the input (rectangular diag) matrix 
    #Only usefule after SVD, since SVD by construction gives a (rectangular diag) matrix
    diagtermMat = numpy.diag(numpy.diag(A))
    
    if len(A) == max( len(A), len(A[0]) ):
        Q = numpy.eye( len(A) , len(A[0]) )
    
    if len(A[0]) == max( len(A), len(A[0]) ):
        p = numpy.eye( len(A), len(A[0]) )
    
    rold=len(numpy.diag(A)) #redefines rold for the next iteration in the loop
    
    return diagtermMat, Q, p, rold
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
Mtilde= Structure(N)
C= Structure(N)
rold=1


print "This code indexes from 0 to N-1  (unlike the notes which go from 0 to N), but clearly same overall process"
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
u = Structure(N+1)

A = Structure(N)
Q = Structure(N)
p = Structure(N)
lam = Structure(N)
#-----defns---------

#-----------Rightmost----------------------------
B[0] = M[0]
print "B[",0,"] is ", B[0]

delta[0], u[1], A[0] = svdLambda(B[0])  # B ---svd---> u delta A
#u[1], delta[0], A[0] = svd(B[0])
#print svd(B[0])
print "delta[",0,"] is ", delta[0]

#-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------
#For rightmost case just want to set lam[0]= delta[0]
A = delta[0]
lam[0] = delta[0]
    
rmax = max( len(A), len(A[0]) ) #should be number of rows in rightmost case
Q[0] = numpy.eye( rmax,rmax )
    
rmin = min( len(A), len(A[0]) ) #should be number of columns in rightmost case
p[0] = numpy.eye( rmin,rmin )
#-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------

#Q[0], lam[0], p[0] = svd( delta[0] )
print "lam[",0,"] is ", lam[0]
print "Q[",0,"] is ", Q[0]
print "u[",1,"] is ", u[1]
print "A[",0,"] is ", A[0]

Mtilde[0] = numpy.dot( p[0] , A[0] )    #numpy.tensordot( p[0] , A[0], axes=2 )
print "Mtilde[",0,"] is ", Mtilde[0]

Z[1] = tripledot( u[1], Q[0], lam[0] )
    
'''print "Z[",1,"] is ", Z[1]
print "len( B[",0,"] ) = ", len(B[0])
print "len( Mtilde[",0,"] ) = ", len(Mtilde[0])
print "len( Z[",0,"] ) = ", len(Z[1])'''

#-----------Rest---------------------------------
for j in range(1,N): # range(1,N) --> [1,...,N-1]
    B[j] = numpy.dot( M[j] , Z[j] )
    print "B[",j,"] is ", B[j]

    delta[j], u[j+1], A[j] = svdLambda(B[j])  # B ---svd---> u delta A
    #u[j+1], delta[j], A[j] = svd(B[j])  # B ---svd---> u delta A
    #print svd(B[j])
    print "delta[",j,"] is ", delta[j]
    #-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------
    #For rightmost case just want to set lam[0]= delta[0]
    A = delta[j]
    lam[j] = numpy.diag(numpy.diag(A))
    
    if len(A) == max( len(A), len(A[j]) ):
        Q[j] = numpy.eye( len(A) , len(A[j]) )
    else:
        Q[j] = numpy.eye( len(A), len(A) )
    
    if len(A[j]) == max( len(A), len(A[j]) ):
        p[j] = numpy.eye( len(A), len(A[j]) )
    else:
        p[j] = numpy.eye( len(A[j]), len(A[j]) )
    #-----------------Secondary Diagonalization----------------------------------------------------------------------------------------------
    
    #lam[j], Q[j], p[j], rold = trivialdiagonalization( delta[j] ,rold)    # delta ---plain decomposition---> Q lam p
    #Q[j], lam[j], p[j] = svd( delta[j] )   # delta ---plain decomposition---> Q lam p
    print "lam[",j,"] is ", lam[j]
    
    Mtilde[j] = numpy.dot( p[j] , A[j] )    #numpy.tensordot( p[j] , A[j], axes=2 )
    print "Mtilde[",j,"] is ", Mtilde[j]

    Z[j+1] = tripledot( u[j+1], Q[j], lam[j] )
    
    print Z[j+1]
    print "len( B[",j,"] ) = ", len(B[j])
    print "len( Mtilde[",j,"] ) = ", len(Mtilde[j])
    print "len( Z[",j,"] ) = ", len(Z[j+1])


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




