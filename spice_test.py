from commands import *
import math
from SpiceParser import SpiceParser

import unittest
class TestSequenceFunctions(unittest.TestCase):
    
    def assertDoublesEqual(self, d1, d2):
        # helper function to test equality for doubles
        # which may not be "exactly" the same
        if d1 == 0.0 or d2 == 0.0:
            if abs(d1-d2) > 1e-6:
                self.assert_(False, "%s != %s" % (d1, d2))
        else:
            if (abs(d1 - d2)) / (abs(d1)+abs(d2)) > 1e-6:
                self.assert_(False, "%s != %s" % (d1, d2))

    def run_netlist(self, netlist):
        # helper function
        parser = SpiceParser()
        return parser.ParseLines( netlist.splitlines() )

    def test_capacitor(self):
        # TEST: verify basic operation of capacitor using
        #       an AC voltage divider at f=2k Hz
        cir = self.run_netlist("""
v1 1 0 ac 1
c1 1 2 1u
c2 2 0 2u
.AC 2k
""")
        freq = 2e3
        z1 = 1/(2*math.pi*freq*1e-6)
        z2 = 1/(2*math.pi*freq*2e-6)
        self.assertDoublesEqual( cir.VF("2"), 1.0 * (z2 / (z1 + z2)) )
        self.assertDoublesEqual( abs(cir.IF("v1")),
                                 1.0 / (z1+z2) )

    def test_inductor(self):
        # TEST: verify basic operationg of inductor using
        #       an AC voltage divider at f=10k Hz
        cir = self.run_netlist("""
v1 1 0 ac 1
l1 1 2 1n
l2 2 0 5n
.AC 100k
""")
        freq = 100e3
        z1 = 2*math.pi*freq*1e-9
        z2 = 2*math.pi*freq*5e-9
        self.assertDoublesEqual( cir.VF("2"), 1.0 * (z2 / (z1 + z2)) )
        x = cir.IF("v1")
        print x
        print cir.phase(x)
        print cir.cartesianToPolar(x)
        self.assertDoublesEqual( abs(cir.IF("v1")),
                                 1.0 / (z1+z2) )
        
    def test_current_source(self):
        cir = self.run_netlist("""
i1 1 0 dc 1
v1 1 0 dc 0
.DC
""")
        self.assertDoublesEqual( cir.IDC("v1"), 1.0 )

    def test_voltage_source(self):
        # TEST: verify basic operation of voltage source
        cir = self.run_netlist("""
v1 1 0 dc 1
.DC
""")
        self.assertDoublesEqual( cir.VDC("1"), 1.0 )

    def test_ccvs(self):
        # TEST: verify operation of the CCVS using a voltage
        #       divider with a known current value
        cir = self.run_netlist("""
v1 1 0 dc 1
r1 1 2 50
r2 2 0 50
h1 3 4 v1 10
r3 3 0 10
r4 4 0 10
.DC
""")
        gain = 10
        self.assertDoublesEqual( cir.IDC("v1"), -1.0/100.0 )
        self.assertDoublesEqual( cir.VDC("3") - cir.VDC("4"),
                                 gain*cir.IDC("v1") )


    def test_vccs(self):
        # TEST: verify operation of the VCCS
        # CHECK: verify that the correct current is flowing;
        #        NOTE: we need to add a voltage source to measure
        #        this since this is not for free in the matrix
        cir = self.run_netlist("""
v1 1 0 dc 1
r1 1 2 50
r2 2 0 50
v2 3 4 dc 0
g1 4 5 2 1 10
r3 3 0 10
r4 5 0 10
.DC 
""")
        self.assertDoublesEqual( cir.VDC("2"), 0.5 )
        self.assertDoublesEqual( cir.VDC("1"), 1.0 )
        self.assertDoublesEqual( cir.IDC("v2"), -5.0 )
        self.assertDoublesEqual( cir.VDC("3"), 50 )
        self.assertDoublesEqual( cir.VDC("5"), -50 )

    def test_vcvs(self):
        # TEST: verify operation of the VCVS by setting up a simple
        #       simple voltage divider and then a VCVS with gain=10
        # CHECK: verify that the correct gain is achieved
        cir = self.run_netlist("""
v1 1 0 dc 1
r1 1 2 50
r2 2 0 50
e1 3 0 2 0 10
r3 3 0 50
.DC
""")
        self.assertDoublesEqual( cir.VDC("2"), 0.5 )
        self.assertDoublesEqual( cir.VDC("3"), cir.VDC("2")*10 )

    def test_resistor(self):
        # TEST: setup a 3 resistor voltage divider and ensure
        #       that all voltages are calculated properly
        # CHECK: voltage at each node
        cir = self.run_netlist("""
r1 1 2 50
r2 2 3 50
r3 3 0 50
v1 1 0 DC 5
.DC
""")
        self.assertDoublesEqual( cir.IDC("v1"), -5.0/(150.0) )
        self.assertDoublesEqual( cir.VDC("1"), 5.0 )
        self.assertDoublesEqual( cir.VDC("2"), 5*(2.0/3.0) )
        self.assertDoublesEqual( cir.VDC("3"), 5*(1.0/3.0) )


if __name__ == '__main__':
    unittest.main()


