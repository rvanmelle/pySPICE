#from sim import SpiceException
from SpiceCircuit import SpiceCircuit, SpiceException
import re

class Callable:
    """Simple class to add support for static class methods."""
    def __init__(self, anycallable):
        self.__call__ = anycallable

class NumberConverter:
    """This class converts a number with SPICE suffix into
    a real float.  It is case insensitive.
    If the suffix is not found, it assumed to be UNITS and is ignored."""

    num_re = re.compile("([+-]?[\de\-\.]+)(\w*)")
    suffix_dict = {'': 1,
		   'g':1e9,
                   'meg':1e6,
                   'k':1e3,
                   'm':1e-3,
                   'u':1e-6,
                   'n':1e-9,
                   'p':1e-12} 

    def isnumber(num):
	return NumberConverter.num_re.match(num) is not None
    isnumber = Callable(isnumber)

    def convert(num):
        match_obj = NumberConverter.num_re.match(num)
        assert match_obj is not None
        number = match_obj.group(1)
        suffix = match_obj.group(2)
        if not NumberConverter.suffix_dict.has_key(suffix):
            warning = "Invalid suffix '%s'. IGNORED!" % suffix
            suffix = ''
            print warning
        return float(number) * NumberConverter.suffix_dict[suffix]        
    convert = Callable(convert)

class SpiceParserException(SpiceException):
    pass

class SpiceParser:
    """This class takes a filename or buffer in a Spice-like format
    and constructs a Circuit object.
    """
    func_re = re.compile("(\w+)\(([^()]*)\)")

    def __init__(self, filename=None):
        self.__Filename = filename
        self.__Circuit = None
        if self.__Filename is not None:
            f = open(filename, 'r')
            self.__Circuit = self.ParseLines(f.readlines())

    def GetCircuit(self):
        return self.__Circuit

    def __convert(self, num):
        # Uses the NumberConverter class to perform
        # simple spice style numeric conversions
        return NumberConverter.convert(num)

    def __ProcessParams(self, params, default=None):
        # Takes a list of numeric named parameters in the form:
        #      x=3.3 y=4.4k...
        # and returns a dictionary of the form:
        #      {'x':3.3, 'y':4400}
        # NOTE: If there is only a single parameter and a default
        # parameter name is passed in, then the default name
        # will be used.
        
        if default is not None:
            assert len(params) == 1
        
        param_dict = {}
        for param in params:
            pieces = param.split('=')
            if len(pieces) == 1:
		if NumberConverter.isnumber(pieces[0]):
		    param_dict[default] = self.__convert(pieces[0])
		else:
		    param_dict[default] = pieces[0]
            else:
		if NumberConverter.isnumber(pieces[1]):
		    param_dict[pieces[0]] = self.__convert(pieces[1])
		else:
		    param_dict[pieces[0]] = pieces[1]

        return param_dict

    def __ProcessDotCard(self, pieces):
        # Process any dotcards which are encountered
        # pieces: list[string] ---> the line broken on whitespace
        command = pieces[0]
        if command == ".end":
            # nothing to do
            # FIXME: should force an exit from parser and ignore any
            #        further input
            pass

        elif command == ".dc":
            # form 1: .dc
            assert len(pieces) == 1
            self.__Circuit.dc()

        elif command == ".tran":
            # form 1: .tran <step> <stop>
	    assert len(pieces) >= 3
	    params = {}
            if len(pieces) > 3:
		params = self.__ProcessParams(pieces[3:])
            step = self.__convert(pieces[1])
            stop = self.__convert(pieces[2])
            self.__Circuit.tran(stop, step, options=params)

        elif command == ".ac":
            # form 1: .ac DEC 10 1 100MEG
            # form 2: .ac 1k
	    # form 3: .ac DEC 10 10k 10MEG method=AWE
            if len(pieces) == 2:
                self.__Circuit.ac(start=self.__convert(pieces[1]))
            else:
                assert len(pieces) >= 5
                assert pieces[1] == "dec"
		params = {}
		if len(pieces) > 5:
		    params = self.__ProcessParams(pieces[5:])
                steps_per_decade = self.__convert(pieces[2])
                start = self.__convert(pieces[3])
                stop = self.__convert(pieces[4])
                self.__Circuit.ac(start=start, stop=stop, 
                                  steps=steps_per_decade, 
				  options=params)

        elif command in [".plot", ".replot", ".print"]:
            # .plot ac vm(3) vp(3)
            # .plot tran v(3) v(4)
            # .plot dc v(7)
            assert len(pieces) >= 3
            assert pieces[1] in ['ac', 'dc', 'tran']
            anal_type = {'ac':self.__Circuit.AC, 
                         'dc':self.__Circuit.DC,
                         'tran':self.__Circuit.Transient}[pieces[1]]
            signal_list = []
            logx = False
            logy = False
            file = None
            for signal in pieces[2:]:
                if signal[0] == ">":
                    file = signal[1:]
                elif signal == "logx":
                    logx = True
                elif signal == "logy":
                    logy = True
                elif self.func_re.match(signal):
                    func_name, signal_spec = self.func_re.match(signal).groups()
		    signal_list.append((func_name, signal_spec))
                else:
                    msg = "Must pass a signal retrieval function: %s" % signal
                    raise SpiceParserException(msg)

	    replot = False
	    for func_name, signal_spec in signal_list:
                if command == ".replot" or replot == True:
                    self.__Circuit.replot( anal_type, func_name, signal_spec,
                                           logx, logy, file)
                elif command == ".plot":
                    self.__Circuit.plot( anal_type, func_name, signal_spec,
                                         logx, logy, file)
		    replot = True
                elif command == ".print":
                    val = self.__Circuit.retrieve_signal( anal_type,
                                                          func_name, 
                                                          signal_spec )
                    print "%s = %s" % (signal, val)

        else:
            msg = "Unrecognized dot card: '%s'. IGNORED!" % command
            raise SpiceParserException(msg)

    def ParseLines(self, lines):
        # lines: list[string] which are the lines of the file to be parsed
        # returns the Circuit object
        #
        # NOTE: Spice traditianally ignored the first line,
        #       but I hate this and simply *will not* do this.
        self.__Circuit = SpiceCircuit()
        for line in lines:
            self.ParseLine(line)
        return self.__Circuit

    def ParseLine(self, line):
        line = line.lower() # lowercase the line to be case insensitive
        pieces = line.split() # splits line on whitespace
        if len(pieces) == 0:
            # empty line, let's skip ahead to the next
            return
            
        first_char = pieces[0][0]
        if first_char == '*':
            #comment, do nothing
            pass

        elif first_char == '.':
            # DOT CARD
            # form: .COMMAND [...]
            self.__ProcessDotCard(pieces)

        elif first_char in ["v", "i"]:
            # voltage source or current source
            # form 1: [V|I]XXX N+ N- NUM
            # form 2: [V|I]XXX N+ N- DC NUM
            # form 3: VXXX N+ N- AC NUM
            # form 4: VXXX N+ N- DC NUM AC NUM
            # form 5: VXXX N+ N- PULSE(V1 V2 TD TR TF PW PER)
            # form 6: VXXX N+ N- AC NUM PULSE(V1 V2 TD TR TF PW PER)
            name, n1, n2 = pieces[0:3]
            cv = self.__convert
            params = {'dc':0.0, 'ac':0.0, 'pulse':None}
            if re.search("pulse\(.*\)", line):
                match_obj = re.search("pulse\((.*)\)", line)
                v1, v2, td, tr, tf, pw, per = match_obj.group(1).split()
                params['pulse'] = ( cv(v1), cv(v2), cv(td), 
                                    cv(tr), cv(tf), cv(pw), cv(per) )
                # now, remove the pulse() portion and continue processing
                line = re.sub("pulse\(.*\)", "", line)
                pieces = line.split()

            if len(pieces) == 4: # form 1
                params['dc'] = cv(pieces[3])
            elif len(pieces) == 5: # form 2 or 3
                assert pieces[3] in params, \
                    '%s not in %s' % (pieces[3], params.keys())
                params[pieces[3]] = cv(pieces[4])
            elif len(pieces) == 7: # form 4
                assert pieces[3] in params, \
                    '%s not in %s' % (pieces[3], params.keys())
                assert pieces[5] in params
                assert pieces[3] != pieces[5]
                params[pieces[3]] = cv(pieces[4])
                params[pieces[5]] = cv(pieces[6])

            if first_char == 'v':
                self.__Circuit.AddVoltageSource(\
                    name, n1, n2, vac=params['ac'], vdc=params['dc'],
                    pulse_args=params['pulse'])
            elif first_char == 'i':
                self.__Circuit.AddCurrentSource(\
                    name, n1, n2, iac=params['ac'], idc=params['dc'],
		    pulse_args=params['pulse'])

        elif first_char in ['r', 'l', 'c']:
            # resistor, inductor, capacitor, voltage, current
            # form: RXXX N+ N- [r=]NUM
            func_map = {'r':self.__Circuit.AddResistor,
                        'l':self.__Circuit.AddInductor,
                        'c':self.__Circuit.AddCapacitor,
                        }
            assert len(pieces) == 4
            name, n1, n2 = pieces[0:3]
            params = self.__ProcessParams(pieces[3:], default=first_char)
            assert params.has_key(first_char)
            func_map[first_char](name, n1, n2, params[first_char])

        elif first_char == 'h':
            # Current controlled voltage source
            # form: HXXX N+ N- DEVNAME GAIN
            assert len(pieces) == 5
            name, out_pos, out_neg, in_dev = pieces[0:4]
            params = self.__ProcessParams(pieces[4:], default='gain')
            assert params.has_key('gain')
            self.__Circuit.AddCCVS(name, out_pos, out_neg, in_dev,
                                   params['gain'])

        elif first_char == 'g':
            # voltage controlled current source
            # form: GXXX N+ N- NC+ NC- VALUE
            assert len(pieces) == 6
            name, out_pos, out_neg, in_pos, in_neg, = pieces[0:5]
            params = self.__ProcessParams(pieces[5:], default='gain')
            assert params.has_key('gain')
            self.__Circuit.AddVCCS(name, in_pos, in_neg, out_pos,
                                   out_neg, params['gain'])

        elif first_char == 'e':
            # Voltage controlled voltage source
            # form: EXXX OUT+ OUT- IN+ IN- GAIN
            assert len(pieces) == 6
            name, out_pos, out_neg, in_pos, in_neg = pieces[0:5]
            params = self.__ProcessParams(pieces[5:], default='gain')
            assert params.has_key('gain')
            self.__Circuit.AddVCVS(name, in_pos, in_neg, out_pos, 
                                   out_neg, params['gain'])

