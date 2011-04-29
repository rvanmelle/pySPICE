from SpiceCircuit import SpiceCircuit

###################################################
# What follows is a command-line interface which is more
# convenience for interactive mode.
###################################################
#### GLOBAL DATA ####
x = SpiceCircuit()

#### COMMANDS ####
vcvs = lambda rn1, rn2, n1, n2, A: x.AddVCVS("", rn1, rn2, n1, n2, A) #DONE
opamp = lambda plus, minus, out, A: x.AddVCVS("", plus, minus, out, 0, A)
vsource = lambda n1, n2, vdc=0.0, vac=0.0: x.AddVoltageSource("", n1, n2, vac=vac, vdc=vdc) #DONE
isource = lambda n1, n2, i, idc=0.0, iac=0.0: x.AddCurrentSource("", n1, n2, iac=iac, idc=idc) #DONE
capacitor = lambda n1, n2, c: x.AddCapacitor("", n1, n2, c) #DONE
inductor = lambda n1, n2, l: x.AddInductor("", n1, n2, l) #DONE
list = lambda: x.Print()
resistor = lambda n1, n2, r: x.AddResistor("", n1, n2, r) #DONE
# ANALYSES
ac = lambda start, stop=-1, steps=10, name="ac": x.ac(start, stop, steps, name)
dc = lambda name="dc": x.dc(name)
# DISPLAY
plot = lambda analsis, *args: x.plot(analysis, args)
# RETRIEVAL
VDC = lambda node: x.VDC( node )
IDC = lambda branch: x.IDC( branch )
VF = lambda node: x.VF( node )
IF = lambda branch: x.VF( branch )
#plot = lambda analysis, *args: x.plot(analysis, args)

#### CONSTANTS ####
def reset():
    global x
    x = SpiceCircuit()


