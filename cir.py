import unittest
from sim import SpiceParser

class TestSequenceFunctions(unittest.TestCase):

    def assertDoublesEqual(self, d1, d2):
	if d1 == 0.0 or d2 == 0.0:
	    if abs(d1-d2) > 1e-6:
		self.assert_(False, "%s != %s" % (d1, d2))
	else:
	    if (abs(d1 - d2)) / (abs(d1)+abs(d2)) > 1e-6:
		self.assert_(False, "%s != %s" % (d1, d2))

    def test_inductor(self):
	netlist = """
v1 1 0 AC 2
l1 1 2 1m
l2 2 0 5m
.AC 1k
"""
	parser = SpiceParser()
	cir = parser.ParseLines(netlist.splitlines())
	self.assertDoublesEqual( abs(cir.IF("v1")), 5.30516e-2)
	self.assertDoublesEqual( abs(cir.VF("1")), 2.0 )
	self.assertDoublesEqual( abs(cir.VF("2")), 2.0*(5.0/6.0) )

    def test_capacitor(self):
	# Tests a basic capacitor divider network
	netlist = """
v1 1 0 AC 1
c1 1 2 1u
c2 2 0 5u
.AC 1MEG 
"""
	parser = SpiceParser()
	cir = parser.ParseLines(netlist.splitlines())
	#cir.plot( cir.AC )
	self.assertDoublesEqual( abs(cir.IF("v1")), 5.23598 )
	self.assertDoublesEqual( abs(cir.VF("1")), 1.0 )
	self.assertDoublesEqual( abs(cir.VF("2")), 1.0/6.0)
	#print "***", abs(cir.VF("1"))
	#print "***", abs(cir.VF("2"))
	#print "***", abs(cir.IF("v1"))

    def test_ccvs(self):
	# Tests a current-controlled voltage source
	netlist = """
*v1 1 0 AC 1
i1 1 0 1
v2 1 2 AC 0
r1 2 0 50
r2 1 3 50
h1 3 0 v2 10
r3 4 0 50
.DC 
"""
	parser = SpiceParser()
	cir = parser.ParseLines(netlist.splitlines())
	#cir.plot( cir.DC )
	#self.assertDoublesEqual( cir.VDC("2"), 5*(2.0/3.0) )

    def test_vcvs(self):
	# testing VCVS
	netlist = """
v1 1 0 AC 1
r1 1 2 50
r2 2 0 50
e1 3 0 2 0 10
r3 3 0 50
.DC
"""
	parser = SpiceParser()
	cir = parser.ParseLines(netlist.splitlines())
	self.assertDoublesEqual( cir.VDC("3"), 5.0 )

    def test_voltage_divider(self):
	# Confirms basic resistor operation
	netlist = """
r1 1 2 50
r2 2 3 50
r3 3 0 50
v1 1 0 AC 5
.DC
"""
	parser = SpiceParser()
	cir = parser.ParseLines(netlist.splitlines())
	self.assertDoublesEqual( cir.VDC("1"), 5.0 )
	self.assertDoublesEqual( cir.VDC("2"), 5*(2.0/3.0) )
	self.assertDoublesEqual( cir.VDC("3"), 5*(1.0/3.0) )

    def test_opamp(self):
	# opamp test
	netlist = """
v1 1 0 AC 1
r1 1 2 50
r2 2 3 50
e1 3 0 0 2 5MEG
.DC    
"""
	parser = SpiceParser()
	cir = parser.ParseLines(netlist.splitlines())
	self.assertDoublesEqual( cir.VDC("2"), 0.0)
	self.assertDoublesEqual( cir.VDC("3"), -1.0)

if __name__ == '__main__':
    unittest.main()


