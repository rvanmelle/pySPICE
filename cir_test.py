from sim import *

def figure2():
    # Figure 2
    isource(1,0,i=1)
    resistor(1,0,r=50)
    resistor(1,2,r=10)
    capacitor(1,0,c=10e-9)
    resistor(2,3,r=10)
    capacitor(3,0,c=10e-9)
    resistor(3,4,r=10)
    capacitor(4,0,c=10e-9)
    resistor(4,5,r=10)
    capacitor(5,0,c=10e-9)
    resistor(5,0,r=50)
    resistor(3,6,r=10)
    capacitor(6,0,c=10e-9)
    resistor(6,7,r=10)
    capacitor(7,0,c=10e-9)
    resistor(7,0,r=50)

def figure3():
    # Figure 3
    vsource(1,0,v=1)
    resistor(1,2,r=9606) # r1a
    resistor(2,0,r=23280) #r1b
    resistor(2,3,r=6800) # r2
    capacitor(3,0,c=20.5e-9) # c2
    capacitor(2,4,c=94.9e-9) # c1
    opamp(3,4,4,A=50000) # opamp1
    resistor(4,5,r=9304) # rg
    opamp(0,5,6,A=50000) # opamp2
    capacitor(5,6,c=15e-9) # c3
    resistor(5,6,r=52107) # rq
    resistor(5,10,r=9304) # r3
    resistor(6,7,r=9304) # r4
    opamp(0,7,8,A=50000) # opamp3
    capacitor(7,8,c=15e-9) # c4
    resistor(8,9,r=20e3) # r
    resistor(9,10,r=20e3) # r
    opamp(0,9,10,A=50000) # opamp4
    ac(start=1, stop=100e6, steps=10)
    #ac(100)
    plot(AC, v(10), v(5), v(3), v(7))
    #dc()
    #plot(DC)

def figure1():
    # Figure 1
    resistor(1,2,r=50)
    capacitor(2,0,c=0.319e-6)
    inductor(2,0,l=0.3176e-6)
    inductor(2,3,l=1.59e-6)
    capacitor(3,4,c=63.72e-12)
    inductor(4,0,l=0.3176e-6)
    capacitor(4,0,c=0.319e-6)
    resistor(4,0,r=50)
    vsource(1,0,v=1)
    ac(start=1, stop=100e6, steps=10)
    plot(AC, v(4))

#figure3()

import unittest
class TestSequenceFunctions(unittest.TestCase):
    
    def setUp(self):
        reset()

    def assertDoublesEqual(self, d1, d2):
	if d1 == 0.0 or d2 == 0.0:
	    if abs(d1-d2) > 1e-6:
		self.assert_(False)
	else:
	    if (abs(d1 - d2)) / (abs(d1)+abs(d2)) > 1e-6:
		self.assert_(False)

    def test_vcvs(self):
	# testing VCVS
	vsource(1,0,v=1)
	resistor(1,2,r=50)
	resistor(2,0,r=50)
	vcvs(2,0,3,0,A=10)
	resistor(3,0,r=50)
	dc()
	self.assertDoublesEqual( VDC(3), 5.0 )

    def test_spice_voltage_divider(self):
	netlist = """
r1 1 2 50
r2 2 3 50
r3 3 0 50
v1 1 0 5
.DC
"""
	parser = SpiceParser()
	cir = parser.ParseLines(netlist.splitlines())
	self.assertDoublesEqual( cir.VDC("2"), 5*(2.0/3.0) )
	self.assertDoublesEqual( cir.VDC("3"), 5*(1.0/3.0) )

    def test_voltage_divider(self):
	# assert_(x in y)
	# assertEqual(x,y)
	# assertRaises(ValueError, random.sample, self.seq)
	resistor(1,2,r=50)
	resistor(2,3,r=50)
	resistor(3,0,r=50)
	vsource(1,0,v=5)
	dc()
	self.assertDoublesEqual( VDC(2), 5*(2.0/3.0) )
	self.assertDoublesEqual( VDC(3), 5*(1.0/3.0) )

    def test_opamp(self):
	# opamp test
	vsource(1,0,v=1)
	resistor(1,2,r=50)
	resistor(2,3,r=50)
	#capacitor(2,3,c=1e-6)
	#inductor(2,3,l=1e-6)
	#resistor(2,3,r=10e3)
	opamp(0,2,3,A=5000000)
	#ac(start=10, stop=100e6, steps=20)
	#plot(AC, v(1), v(2), v(3))
	dc()
	#plot(DC)
	self.assertDoublesEqual( VDC(2), 0.0)
	self.assertDoublesEqual( VDC(3), -1.0)

if __name__ == '__main__':
    unittest.main()

#list()
#x = solve(freq=0.0)
#print x
#ac(start=10, stop=100e6, steps=20)
#plot(AC, v(1), v(2), v(3))
#plot(AC, v(4), v(6), v(8), v(10))
#plot(AC, v(4), v(3), v(1), v(10))

#dc()
#plot(DC)

#plot(DC, v(1), v(2), v(3))
#plot([v(4)])
#plot(v(3))

