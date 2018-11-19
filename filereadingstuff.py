import math
import numpy
#from numpy.linalg import svd

#datavec = list(numpy.loadtxt('PCV1_1genesample.txt')) + list(numpy.loadtxt('PCV1_2genesample.txt')) + list(numpy.loadtxt('PCV1_3genesample.txt')) + list(numpy.loadtxt('PCV1_4genesample.txt'))
#Now we want to concatenate the datavector a bit 
gene1 = list(numpy.loadtxt('PCV1_Sample_1.txt'))
#gene2 = list(numpy.loadtxt('PCV1_2genesample.txt'))
#gene3 = list(numpy.loadtxt('PCV1_3genesample.txt')) 
#gene4 = list(numpy.loadtxt('PCV1_4genesample.txt'))
#Name them here for ease but would later cut out this step 

#print gene1
#print range(1,50)

kthtxtfilename = 'PCV1_Sample_k.txt'

print kthtxtfilename
print kthtxtfilename[-1]
print kthtxtfilename[-5]

print 'PCV1_Sample_' + str(8*8)  + '.txt' #<------------------------------IMPORTANT!!! How to make an int input into a string file name format 
print ('PCV1_Sample_' + str(2)  + '.txt')
#print str(8*8)


# list(numpy.loadtxt('PCV1_Sample_' + str(2)  + '.txt'))

#locals() = 4
#print 'gene_' + str(2) 
#print gene_100 

'''d = {}   #dictionary, expleantion found at https://stackoverflow.com/questions/8028708/dynamically-set-local-variable 
d['xyz'] = 7*9
print(d['xyz'])

#d = {}   #dictionary, expleantion found at https://stackoverflow.com/questions/8028708/dynamically-set-local-variable 
d[('gene_' + str(100))] = 10*11
print(d['gene_100'])

d = {}   #dictionary, expleantion found at https://stackoverflow.com/questions/8028708/dynamically-set-local-variable 
d[('gene_' + str(2))] = list(numpy.loadtxt('PCV1_Sample_' + str(2)  + '.txt'))
print(d['gene_2'])'''

#print list(numpy.loadtxt('PCV1_Sample_' + str(2)  + '.txt'))

NumSampleMin = 1
NumSampleMax = 3 
Ngenes = (NumSampleMax + 1) - NumSampleMin

datavec = [None]*0  #Makes it begin as an empty set, since we use a sum in a for loop, if this is not empty it will add extraneous values (e.g. using [None]*3 adds 3 None's to the final datavec 
for k in range(NumSampleMin, (NumSampleMax + 1) ):
    datavec += list(numpy.loadtxt('PCV1_Sample_' + str(k) + '.txt'))

#print datavec
#print [None]*0

print range(NumSampleMin, (NumSampleMax + 1) )
print len(datavec)

print 'The length of this datavec, len(datavec) = ', len(datavec), ', a value which should be equal to 1759*Ngenes = ', 1759*Ngenes
