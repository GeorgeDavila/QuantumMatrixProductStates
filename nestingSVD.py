import math
import numpy 
from numpy.linalg import svd

TestMatrix = [
    [2, 5, 3],
    [1, 2, 1],
    [4, 1, 1],
    [3, 5, 2],
    [5, 3, 1],
    [4, 5, 5],
    [2, 4, 2],
    [2, 2, 5],
]

b= 2

def repeated(f, n):
    def rfun(p):
        return reduce(lambda x, _: f(x), xrange(n), p)
    return rfun

def square(x):
    print "square(%d)" % x
    return x * x

#A = TestMatrix
def SVDfunc(A):
    print "SVDfunc(%s)" % A[3]
    return svd(A[3])

print repeated(square, 5)(3)


print repeated(SVDfunc, b)(TestMatrix)


#The _ is convention for a variable whose value you don't care about.
