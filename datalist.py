#Want to make a file which contains the DNA data that I can call as a list in python
#Also should export these results to pdfs for better viewing
import math
import numpy

datavec = list(numpy.loadtxt('PCV1_1genesample.txt')) + list(numpy.loadtxt('PCV1_2genesample.txt')) + list(numpy.loadtxt('PCV1_3genesample.txt')) + list(numpy.loadtxt('PCV1_4genesample.txt'))
#We first call the values from text files using numpy.loadtxt('PCV1_1genesample.txt') then make them into lists, then combine them into one big list called datavec
#NOTE: can do this same proecess for .biz or .gz or other text files, just need to be unpacked a bit more 
print datavec
print len(datavec)