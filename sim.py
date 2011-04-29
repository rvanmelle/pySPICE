#!/usr/bin/env python
import cmd
from SpiceParser import SpiceParser
from SpiceCircuit import SpiceCircuit, SpiceException

class SpiceInteractive(cmd.Cmd):
    def __init__(self, parser, circuit):
        cmd.Cmd.__init__(self)
        self.prompt = "spice% "
        self.__Parser = parser
        self.__Circuit = circuit

    def do_show(self, line):
        print "ANALYSES:", self.__Circuit.GetAnalyses()
        print "BRANCHES:", self.__Circuit.GetBranches()
        print "NODES:", self.__Circuit.GetNodes()

    def do_reset(self, line):
        self.__Circuit.reset_plotter()

    def default(self, line):
        try:
            self.__Parser.ParseLine("."+line)
        except SpiceException, e:
            print str(e)

    def do_EOF(self, line):
        print "Thank you! Exiting.\n"
        sys.exit(0)


if __name__=="__main__":
    # command line version
    import sys
    file_name = sys.argv[1]
    x = SpiceParser(file_name)
    cir = x.GetCircuit()
    interpreter = SpiceInteractive(x, x.GetCircuit())
    interpreter.cmdloop("Welcome to spice!")
    



  
