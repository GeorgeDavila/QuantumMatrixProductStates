import math
import numpy 
#x = input = len(DNAdatavector)

print numpy.linalg.eigvals([[1,2,3],
    [4 ,5,6],
    [7 ,8,9]])

a = numpy.linalg.eigvals([[1,2,3],
    [4 ,5,6],
    [7 ,8,9]])

a[ abs(a) < 0.00001 ] = 0 #Sets a threshold so that any value with abs(a) < 0.00001 just goes to 0, since in this instance the last term should actually be 0 

print a


