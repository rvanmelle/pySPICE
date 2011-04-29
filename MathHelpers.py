import numpy
import numpy.linalg
import Gnuplot
import cmath, math, copy

def step_function(t):
    # This function implements the step function
    if t > 0.0:
	return 1.0
    return 0.0

def pole_residue_solver(t, c, poles, residues):
    # This function solves a time domain function
    # of the form that corresponds to a pole/residue
    # inverse Laplace transform.
    #
    # t : time
    # c : c[0]*delta(t) + c[1]*u(t) + c[2]*t + c[3]*t^2 + ...
    # poles: p1, p2, p3, ...
    # residues: r1, r2, r3, ...
    assert len(poles) == len(residues)
    sum = 0.0

    for i in range(len(c)):
	if i == 0 and t == 0.0:
	    sum += c[0]
	elif i == 1:
	    sum += c[1]
	elif i == 2:
	    sum += c[i]*t
	elif i > 2:
	    assert 0, "%s %s %s" % (c, poles, residues)

    for i in range(len(poles)):
	sum += residues[i] * cmath.exp(poles[i]*t)

    return sum

def integrate(c, p, r, impulse=None, scaling=1):
    # equivalent to multiplying by 1/s
    # f(t) = c[0] + c[1]*t + c[2]*t^2 + r[0]*exp(p[0]*t) + 
    #        r[1]*exp(p[1]*t)
    assert len(p) == len(r)
    new_c = []
    new_r = []
    sum = 0.0
    # special case of impulse
    if impulse is not None:
	sum += impulse

    # loop through the time coefficients
    for i in range(len(c)):
	new_c.append(c[i] / (i+1))

    # loop through the poles and residues
    for i in range(len(p)):
	new_r.append(r[i] / p[i])
	sum -= scaling*new_r[-1]
    
    # add the new scalar to the beginning of the list
    new_c.insert(0, sum)

    return new_c, copy.copy(p), new_r

def zroots(coef, polish=1, sort=1):
    # finds and returns all of the roots for the 
    # given polynominal
    EPS = 2e-6
    MAXM = 100
    m = len(coef) - 1
    roots = range(m)

    ad = copy.copy(coef)
    for j in range(m, 0, -1):
	x = 0+0.0j
	x = laguer(ad, j, x)
	if abs(x.imag) <= 2.0*EPS*abs(x.real):
	    x = x - x.imag
	roots[j-1] = x
	b = ad[j]
	for jj in range(j-1, -1, -1):
	    c = ad[jj]
	    ad[jj] = b
	    b = x*b + c

    if polish:
	for j in range(1, m+1):
	    roots[j-1] = laguer(coef, m, roots[j-1])
	
    if sort:
	for j in range(2, m+1):
	    x = roots[j-1]
	    for i in range(j-1, 0, -1):
		if roots[i-1].real <= x.real:
		    break
		roots[i] = roots[i-1]
	    roots[i] = x

    return roots
	
def laguer(coef, m, x):
    # Returns a single coefficient
    MAXIT = 80
    MT = 10
    EPSS = 1.0e-7
    frac = [0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0]

    for iter in range(1, MAXIT+1):
	b = coef[m]
	err = abs(b)
	d = 0+0j
	f = 0+0j
	abx = abs(x)
	for j in range(m-1, -1, -1):
	    f = d + x*f
	    d = b + x*d
	    b = coef[j] + x*b
	    err = abs(b) + abx*err
	err *= EPSS
	if abs(b) <= err:
	    return x

	g = d / b
	g2 = g * g
	h = g2 - 2.0*(f/b)
	sq = cmath.sqrt(float(m-1)*(m*h - g2))
	gp = g + sq
	gm = g - sq
	abp = abs(gp)
	abm = abs(gm)
	if (abp < abm):
	    gp = gm
	if max(abp, abm) > 0.0:
	    dx = (m+0.0j) / gp
	else:
	    dx = cmath.exp(cmath.log(1+abx)) * complex(math.cos(float(iter)),
						      math.sin(float(iter)))
	x1 = x - dx
	if x.real == x1.real and x.imag == x1.imag:
	    return x
	
	if iter % MT:
	    x = x1
	else:
	    x = x - frac[iter/MT] * dx
	    
    print "too many iterations in laguerre"
    return x

def my_pade(coef, complex=True):
    L = (len(coef) / 2) - 1
    M = L + 1
    assert len(coef) == L + M + 1

    rows = []
    rhs = []
    for i in range(M):
	rows.append([])
	if i >= M:
	    rhs.append(0.0)
	else:
	    rhs.append( -coef[M+i] )
	for j in range(M):
	    if i+j >= len(coef):
		rows[-1].append(0.0)
	    else:
		rows[-1].append( coef[i+j] )

    matrix = numpy.array(rows)
    rhs_matrix = numpy.array(rhs)
    b_coef = list(numpy.linalg.solve_linear_equations( matrix, rhs_matrix ))
    b_coef.reverse()

    if complex == True:
	b_coef.insert(0, 1.0+0j)
    else:
	b_coef.insert(0, 1.0)

    a_coef = []
    for i in range(1, M+1):
	new_coef = 0.0
	for j in range(i):
	    new_coef += b_coef[j]*coef[i-j-1]
	a_coef.append(new_coef)

    return (numpy.array(a_coef), numpy.array(b_coef))

def my_pade2(coef):
    L = len(coef) / 2
    M = L + 1
    assert len(coef) == L + M

    rows = []
    rhs = []
    for i in range(L, M+1):
	rows.append([])
	rhs.append( -coef[i+1] )
	for j in range(0, L):
	    rows[-1].append( coef[i-j] )

    matrix = numpy.array(rows)
    rhs_matrix = numpy.array(rhs)
    b_coef = list(numpy.linalg.solve_linear_equations( matrix, rhs_matrix ))
    b_coef.insert(0, 1.0+0j)

    a_coef = []
    for i in range(1, M+1):
	new_coef = 0.0
	for j in range(i):
	    new_coef += b_coef[j]*coef[i-j-1]
	a_coef.append(new_coef)

    return (numpy.array(a_coef), numpy.array(b_coef))

def pade(coef, complex=True):
    n = len(coef) / 2
    assert len(coef) == 2*n+1
    
    y = []
    rows = []
    for j in range(1, n+1):
	y.append(-coef[n+j])
	rows.append([])
	for k in range(1, n+1):
	    rows[-1].append(coef[j-k+n])

    rhs = numpy.array(y)
    matrix = numpy.array(rows)

    solution = numpy.linalg.solve_linear_equations(matrix, rhs)

    b_coef = solution
    a_coef = [coef[0]]
    for k in range(1,n+1):
	sum = coef[k]
	for j in range(1,k+1):
	    sum += solution[j-1]*coef[k-j]
	a_coef.append(sum)

    b_coef = solution
    b_coef = list(b_coef)
    if complex:
	b_coef.insert(0, 1.0+0j)
    else:
	b_coef.insert(0, 1.0)

    return (numpy.array(a_coef), numpy.array(b_coef))

def solve_poly(c_num, c_den, x):
    num = 0.0
    den = 0.0

    for i in range(len(c_num)):
	num += c_num[i] * pow(x, i)

    for i in range(len(c_den)):
	den += c_den[i] * pow(x, i)

    if den != 0.0:
	return num / den
    return num

from math import *
def func(x):
    return pow( 7 + pow(1+x, 4.0/3.0), 1.0/3.0 )

if __name__=="1__main__":
    print laguer([3,4,1], 2, 100+0.0j)
    print zroots([3,4,1])
    print zroots([-27, -72, -6, 1])

if __name__=="__main__":
    c = []
    p = [-1.0, -2.0]
    r = [2.0,1.0]
    print "c=", c
    print "p=", p
    print "r=", r

    c, p, r = integrate(c,p,r, impulse=3)
    print "\nc=", c
    print "p=", p
    print "r=", r

    c,p,r = integrate(c,p,r)
    print "\nc=", c
    print "p=", p
    print "r=", r

if __name__=="h__main__":

    def my_func(x):
	return math.exp(-x)

    poly = [1.0, -1.0, 0.5, -1.0/6, 1.0/24, -1.0/120]
    a_coef, b_coef = my_pade(poly, complex=False)

    p = Gnuplot.Gnuplot()
    p("set terminal X11")
    p("set logscale y")
    p("set grid")

    res1, res2, res3 = [], [], []
    n = 1000
    for i in range(n):
	x = i*(10.0/n)
	res1.append((x, solve_poly(poly, [], x)))
	res2.append((x, solve_poly(a_coef, b_coef, x)))
	res3.append((x, my_func(x)))
	    
    d = Gnuplot.Data(res1, title='Taylor',
		     with='lines')
    p.plot(d)
    d = Gnuplot.Data(res2, title='pade', with='lines')
    p.replot(d)
    d = Gnuplot.Data(res3, title='func', with='lines')
    p.replot(d)
    #d = Gnuplot.Data(res4, title='my_pade', with='lines')
    #p.replot(d)

    import time
    time.sleep(1000)

if __name__=="t__main__":
    
    poly = [2.0, 1.0/9, 1.0/81, -49.0/8747, 175.0/78732, 0.000356]
    #a_coef, b_coef = pade(poly, complex=False)
    #print "a_coef=", a_coef
    #print "b_coef=", b_coef

    a_coef2, b_coef2 = my_pade(poly, complex=False)
    print "my_a_coef=", a_coef2
    print "my_b_coef=", b_coef2
    
    p = Gnuplot.Gnuplot()
    p("set terminal X11")
    p("set grid")

    res1, res2, res3, res4 = [], [], [], []
    n = 1000
    for i in range(n):
	x = i*(10.0/n)
	res1.append((x, solve_poly(poly, [], x)))
	#res2.append((x, solve_poly(a_coef, b_coef, x)))
	res3.append((x, func(x)))
	res4.append((x, solve_poly(a_coef2, b_coef2, x)))
	    
    d = Gnuplot.Data(res1, title='Taylor',
		     with='lines')
    #p.plot(d)
    #d = Gnuplot.Data(res2, title='pade', with='lines')
    #p.replot(d)
    d = Gnuplot.Data(res3, title='func', with='lines')
    p.replot(d)
    d = Gnuplot.Data(res4, title='my_pade', with='lines')
    p.replot(d)

    import time
    time.sleep(1000)
