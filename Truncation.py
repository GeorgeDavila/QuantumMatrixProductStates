import math
import numpy 
from numpy.linalg import svd


rho=[[4,0,0],[0,3.14159,0],[0,0,3]]

a = numpy.linalg.eigvals(rho)

print a

a[ abs(a) < 3.14 ] = 0 

print a 