from numpy import *
import numpy.linalg

x = array([[92.477, 10.202, -28.832], 
	   [1.963, 32.816, 62.414],
	   [26.821, 36.816, 57.234]])
	   #[23.2134, -86.3925, 44.693]])

print x

v, s, wt = numpy.linalg.singular_value_decomposition(x, 1)

print "v *********"
print v
print "s *********"
print s
print "wt *********"
print wt
print "*************"

cond = max(s) / min(s)
rcond = 1 / cond

print "rcond=", rcond
