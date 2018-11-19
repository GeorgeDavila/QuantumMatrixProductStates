import math
import numpy 
#Can also do some concatenations using itertools package 

datavec = list(numpy.loadtxt('PCV1_1genesample.txt')) + list(numpy.loadtxt('PCV1_2genesample.txt')) + list(numpy.loadtxt('PCV1_3genesample.txt')) + list(numpy.loadtxt('PCV1_4genesample.txt'))

''' print list(numpy.loadtxt('PCV1_1genesample.txt'))[5-1] #Remember A[4] gives the FIFTH 5th element of list A, since in python A[0] is the first element of the list 
print list(numpy.loadtxt('PCV1_1genesample.txt'))[5%4] #Remember A[4] gives the FIFTH 5th element of list A, since in python A[0] is the first element of the list

datavec = list(numpy.loadtxt('PCV1_1genesample.txt')) + list(numpy.loadtxt('PCV1_2genesample.txt')) + list(numpy.loadtxt('PCV1_3genesample.txt')) + list(numpy.loadtxt('PCV1_4genesample.txt'))
gene1 = list(numpy.loadtxt('PCV1_1genesample.txt'))
gene2 = list(numpy.loadtxt('PCV1_2genesample.txt'))
gene3 = list(numpy.loadtxt('PCV1_3genesample.txt')) 
gene4 = list(numpy.loadtxt('PCV1_4genesample.txt'))

Ngenes = 4    
x = len(datavec)

#Can do by slicing: 
list1 = ['f', 'o', 'o']
list2 = ['hello', 'world']
result = [None]*(len(list1)+len(list2))
result[::2] = list1
result[1::2] = list2
print result

print range(0,10)[1::2] # gives [1, 3, 5, 7, 9] every 2nd number starting from 1 
print range(0,10)[2::2] # gives [2, 4, 6, 8] every 2nd number starting from 2 
print range(0,10)[1::3] # gives [1, 4, 7] every 3rd number starting from 1
print range(0,10)[::3]# gives [0, 3, 6, 9] every 3rd number starting from beginning of list 


datavec[::4]
datavec[1::4]
datavec[2::4]
datavec[3::4]

catdatavec = [None]*(x) #x defined as len(datavec), so defines a list of same length as the datavec (for now filled with None's)

catdatavec[::4] = gene1
catdatavec[1::4] = gene2
catdatavec[2::4] = gene3
catdatavec[3::4] = gene4 '''

#print catdatavec

#-------------------Filling of the catenated datavector (catenated first then filled)--------------------------------------------------------
x = len(datavec)

b = math.ceil(math.log( x ) / float(math.log( 2 )))

filler = list((0*numpy.arange(((2**b) - x)) + 4))


gene1 = list(numpy.loadtxt('PCV1_1genesample.txt'))
gene2 = list(numpy.loadtxt('PCV1_2genesample.txt'))
gene3 = list(numpy.loadtxt('PCV1_3genesample.txt')) 
gene4 = list(numpy.loadtxt('PCV1_4genesample.txt'))

catdatavec = [None]*(x) #x defined as len(datavec), so defines a list of same length as the datavec (for now filled with None's)

catdatavec[::4] = gene1
catdatavec[1::4] = gene2
catdatavec[2::4] = gene3
catdatavec[3::4] = gene4

datavec=catdatavec #<<<------------------Redefines datavec as the catenated one 
#print datavec 

if (2**b) == x:
    fdatavec = datavec
elif (2**b) > x:
    fdatavec = datavec + filler
else:
    print "Something went wrong. Probably didn't define the data vector as needed. Try again"
    
print fdatavec 
print len(fdatavec)

#---------------------------------------------------------------

print "The result is the same whether or not we catenate it before or after filling it, so serves as a check that we can do these in whatever order we choose" 
print "(for data samples of the same length in which the fillings are placed in the same spots ) "

#Now this does exactly what it's meant for 

