import math
import numpy 


# from   http://www.afterhoursprogramming.com/tutorial/Python/Writing-to-Files/
f = open("txttest.txt","w") #opens file with name of "test.txt"
f.write("I am a test file.")
f.write("Maybe someday, he will promote me to a real file.")
f.write("Man, I long to be a real file")
f.write("and hang out with all my new real file friends.") 
f.close()

list = [1, 2, 3]
stringRepr = str(list).translate(None, '[,]') # string.translate fn should remove the [ , and ] from the string    See # https://docs.python.org/2/library/stdtypes.html#str.translate

print "Look at this test line removing vowels with string.translate"
print "Look at this test line removing vowels with string.translate".translate(None, 'aeiou')

#  https://stackoverflow.com/questions/2906092/converting-a-list-to-a-string

f = open("listastxttest.txt","w") #opens file with name of "test.txt"
f.write(stringRepr)
f.close()

