import numpy
import numpy.linalg
import math, cmath
import Gnuplot
import MathHelpers
import copy

class SpiceException(Exception):
    pass

class CircuitException(SpiceException):
    pass

class SpiceCircuit:
    """This is a simple circuit representation.  It performs basic SPICE
    style analyses such as AC, DC, Transient.  It is hooked to Gnuplot
    which it uses for a native plotter.  Matrix solutions are obtained
    using the numpy package.
    """

    #############
    #  CONSTANTS
    #############

    # analysis types
    AC="ac"
    DC="dc"
    Transient="tran"

    # signal types
    Voltage = "Voltage"
    Current = "Current"
    Condition = "Cond"

    # transient solution methods
    Trapezoidal = "trap"
    BackwardEuler = "euler"

    # AC solution methods
    MNA = "mna"
    AWE = "awe"

    # HELPER CLASSES
    class XYFunc:
        def __init__(self, func):
            self.__func = func
        def __call__(self, num):
            return (num[0], self.__func(num[1]))

    def __init__(self):
        self.__Size = 0
        self.__NodeMap = {}
        self.__G = {} # The actual G matrix data
        self.__G_array = None # cache
        self.__C = {} # The actual C matrix data
        self.__C_array = None # cache
        self.__X = [] # The actual X array data
        self.__X_array = None # cache
        self.__B = {} # The actual B array data
        self.__B_array = None # cache
        self.__Branches = {} # The names of branch currents
        self.__Conditions = [] # Conditions numbers

        # Global data during simulations
        self.__Time = None # current time during transient
        self.__CurAnalysis = None #current type of analysis
        
        # Settings
        self.__RecordCondition = True

	self.NUM_AWE = 0
	self.Scaling = 1

        # For plotting
        self.__Plotter = Gnuplot.Gnuplot() 

        # Holds all of the results from simulations
        self.__Results = {}

    def __getattr__(self, name):
        # This provides more convenient access to 
        # various data structures... slists etc
        # are converted into matrices and arrays in
        # this process as well.
        if name == "Size":
            return self.__Size

        elif name == "Time":
            return self.__Time

        elif name == "G":
            # Covert sparse matrix representation into matrix
            if self.__G_array is None:
                rows = []
                for i in range(1,self.__Size+1):
                    rows.append([]) 
                    for j in range(1,self.__Size+1):
                        if self.__G.has_key((i,j)):
                            rows[-1].append(self.__G[(i,j)])
                        else:
                            rows[-1].append(0.0)
                self.__G_array = numpy.array(rows)
            return self.__G_array

        elif name == "C":
            # Covert sparse matrix representation into matrix
            if self.__C_array is None:
                rows = []
                for i in range(1,self.__Size+1):
                    rows.append([]) 
                    for j in range(1,self.__Size+1):
                        if self.__C.has_key((i,j)):
                            rows[-1].append(self.__C[(i,j)])
                        else:
                            rows[-1].append(0.0)
                self.__C_Array = numpy.array(rows)
            return numpy.array(rows)

        elif name == "X":
            if self.__X_array is None:
                x_array = numpy.array(self.__X)
                self.__X_array = numpy.transpose(x_array)
            return self.__X_array

        elif name == "Cond":
            # FIXME: could cache this here
            return numpy.array(self.__Conditions)

        elif name == "B":
            # The B array contains the voltage and current
            # source data.  These may be time dependent sources
            # so the B array needs to be re-evaluated at every
            # point in time.
            # FIXME: Make this far more efficient by only
            # re-evaluating if necessary and certainly only
            # those values which are actually functions
            if self.__B_array is None or self.Time is not None:
                b = []
                for i in range(1, self.__Size+1):
                    if self.__B.has_key(i):
                        # Eval the B function at time self.Time
                        b.append(self.__B[i](self.Time, self.__CurAnalysis))
                    else:
                        b.append(0.0)
                self.__B_array = numpy.transpose(numpy.array([b]))
            return self.__B_array

        else: raise AttributeError, name

    def __setattr__(self, name, value):
        if name == "Time":
            self.__Time = value
            self.__B_array = None
        else:
            self.__dict__[name] = value

    def __Solve(self, freq=0.0):
        # Solve circuit at a particular frequency.
        # This solver uses the built-in lineary equation solver
        # from the numpy package.
        if freq == 0.0:
            total_mna = self.G
        else:
            total_mna = self.G + 2*math.pi*freq*1j*self.C

        solution = numpy.linalg.solve_linear_equations(total_mna, self.B)
        if self.__RecordCondition:
            v, s, wt = numpy.linalg.singular_value_decomposition(total_mna)
            cond = max(s) / min(s)
            self.__Conditions.append((freq,cond))
        return solution

    def reset_plotter(self):
        self.__Plotter = Gnuplot.Gnuplot()

    ##############################
    #### ANALYSIS SECTION #######
    def SetAnalysis(self, newAnalysis):
        self.__CurAnalysis = newAnalysis
        self.__B_array = None

    def dc(self, name="dc"):
        # Run a DC analysis by solving at freq=0.0
        self.__CurAnalysis = self.DC
        self.__Results[name] = self.__Solve(0.0)
        self.__CurAnalysis = None

    def __pade_approximation(self, nodeName, L, M, debug=0):
        assert L == M - 1
        num_moments = L+M+1
        if debug: "DATA L=%s M=%s num_moments=%s" % (L, M, num_moments)

        if debug: print "num_moments", num_moments
        #scaling = 766423.0
        node_index = self.TranslateNode(nodeName)
        node_moments = []
        # step 1: calculate the Moments
        g_inverse = numpy.linalg.inverse(self.G)
        if debug: print "g_inverse", g_inverse
        last_moment = numpy.matrixmultiply(g_inverse, self.B)
        if debug: print "last_moment", last_moment, node_index
        node_moments.append(last_moment[node_index-1][0])

		# test commit
        for i in range(num_moments-1):
            intermediate = -1 * numpy.matrixmultiply(g_inverse, self.C)
            last_moment = numpy.matrixmultiply( intermediate, last_moment )
            moment = self.Scaling * last_moment[node_index-1][0]
            node_moments.append(moment)
            last_moment = self.Scaling * last_moment
            if debug: print "last_moment", last_moment

        print "suggested scaling =", node_moments[0]/node_moments[1]
        if debug: print "node_moments=", node_moments

        # Call the general pade algorithm
        a_coef, b_coef = MathHelpers.my_pade(node_moments, complex=False)

        # Return results
        return a_coef, b_coef, node_moments

    def __ac_awe(self, name, nodeName, L, M, start, stop, steps):

        a_coef, b_coef, node_moments = self.__pade_approximation(nodeName, L, M)

        def solve(a_coef, b_coef, freq):
            numerator = 0+0j
            denominator = 0+0j
            freq = freq / (self.Scaling)

            for i in range(len(a_coef)):
                numerator += a_coef[i] * pow(2*math.pi*freq*1j,i)
            
            for i in range(len(b_coef)):
                denominator += b_coef[i] * pow(2*math.pi*freq*1j, i)

            return numerator / denominator

        # evaluate the function at the specified frequencies
        result_vec = []
        interval = 1.0 / steps
        num_points = int((math.log10(stop) - math.log10(start)) / interval)
        for point in range(num_points):
            freq = math.pow(10, math.log10(start) + point * interval)
            result_vec.append( (freq, abs(solve(a_coef, b_coef, freq))) )

        g = self.__Plotter
        g("set output \"data.ps\"")
        g("set terminal postscript color")
        g("set terminal postscript solid")

        d = Gnuplot.Data(result_vec,
                         title='AWE(%s/%s)' % (L,M),
                         with='lines')
        g.replot(d)

    def __ac_normal(self, name, start, stop, steps):
        if stop == -1:
            # this indicates that we are going to do a single frequency
            # at start
            self.__Results[name].append((start, self.__Solve(start)))
        else:
            interval = 1.0 / steps
            num_points = int((math.log10(stop) - math.log10(start)) / interval)
            for point in range(num_points):
                freq = math.pow(10, math.log10(start) + point * interval)
                self.__Results[name].append((freq, self.__Solve(freq)))

    def ac(self, start, stop=-1, steps=10, name="ac", options={}):
        # start : starting frequency
        # stop : final frequency
        # steps: steps per decade
        self.__Results[name] = []
        self.__CurAnalysis = self.AC
        method = self.MNA
        if options.has_key('method'):
            method = options['method']
        assert method in [self.MNA, self.AWE]

        if method == self.MNA:
            self.__ac_normal(name, start, stop, steps)

        elif method == self.AWE:
            assert options.has_key('node')
            assert options.has_key('order')
            node = str(int(options['node']))
            M = int(options['order'])
            L = int(M-1)
            self.__ac_awe(name, node, L, M, start, stop, steps)

        self.__CurAnalysis = None

    def tran(self, stop, step, name="tran", method=1, options={}):
        self.__Results[name] = []
        self.__CurAnalysis = self.Transient
        method = self.MNA
        if options.has_key('method'):
            method = options['method']
        assert method in [self.MNA, self.AWE]

        if method == self.MNA:
            self.__tran_normal(stop, step, name, method)
        elif method == self.AWE:
            assert options.has_key('node')
            assert options.has_key('order')
            assert options.has_key('scaling')
            assert options.has_key('mode')
            assert options.has_key('size')
            assert options.has_key('tr')
            assert options.has_key('width')
            node = str(int(options['node']))
            M = int(options['order'])
            L = int(M-1)
            self.Scaling = options['scaling']
            mode = options['mode']
            size = options['size']
            tr = options['tr']
            width = options['width']
            self.__tran_awe(name, node, L, M, stop, step, mode, size, 
                            tr, width)

        self.__CurAnalysis = None
        # back to non time-domain
        self.Time = None

    def __tran_awe(self, name, nodeName, L, M, stop, step, mode, size,
                   tr, width):
        # stop : stopping time
        # step: size of the time step
        # L: order of the numerator
        # M: order of the denominator
        # mode: one of "IMPULSE", "IDEAL_STEP", "STEP", "PULSE", "RAMP"
        # size: height of the step or pulse
        # tr: rise and fall time of the STEP or PULSE
        # width: width of the PULSE
        
        self.__CurAnalysis = self.AC
        a_coef, b_coef, node_moments = self.__pade_approximation(\
            nodeName, L, M, debug=0)

        # now, calculate the poles of the sequence
        poles = MathHelpers.zroots(1j*b_coef, sort=0)
        print "L=%s M=%s POLES=%s" % (L, M, poles)

        # calculate the residues of the sequence
        # using the formula given in the paper
        p_rows = []
        for i in range(M+1):
            p_rows.append([])
            for j in range(M):
                p_rows[i].append(poles[j]**-(i+1))
            if i == 0:
                p_rows[i].append(-1+0j)
            else:
                p_rows[i].append(0+0j)

        p_matrix = numpy.array(p_rows)
        rhs = -1*numpy.array(node_moments[0:L+2])
        residues = numpy.linalg.solve_linear_equations( p_matrix, rhs )

        #NOTE: residues are in the form: [k1, k2, ... km, c]
        #    : poles are in the form: [p1, p2, ... pm]
        # i.e. H(s) = c + (k1 / (s+p1)) + (k2 / (s+p2))  + ...
        # L-1{ H(s) } = c*sigma(t) + k1*exp(-p1*t) + k2*exp(-p2*t) = h(t)
        # L-1{ H(s) / s } = INTEGRAL{ h(t) }
            
        # The residues actually contains both the actual residues
        # and the direct coupling c
        impulse = residues[-1]
        residues = residues[0:-1]

        # Now we calculate all of the responses that we may need
        # first, record the impulse response
        c0, poles0, residues0 = [], copy.copy(poles), copy.copy(residues)
        c0.insert(0, impulse)

        # integrate once to get step response and record
        c, poles, residues = MathHelpers.integrate([], poles, residues, impulse)
        c1, poles1, residues1 = copy.copy(c), copy.copy(poles), \
            copy.copy(residues)
        c1.insert(0, 0.0)

        # integrate again to get the ramp response and record
        c2, poles2, residues2 = MathHelpers.integrate(c, poles, residues)
        c2.insert(0, 0.0)

        # Setup some shorthands and calculate required values
        impulse_response = lambda t: MathHelpers.pole_residue_solver(\
            t*self.Scaling, c0, poles0, residues0)
        step_response = lambda t: MathHelpers.pole_residue_solver(\
            t*self.Scaling, c1, poles1, residues1)
        ramp_response = lambda t: MathHelpers.pole_residue_solver(\
            t*self.Scaling, c2, poles2, residues2)
        u = MathHelpers.step_function
        slope = (size / tr) / self.Scaling

        # Calculate the answers by looping through the time
        # steps
        result = []
        for i in range(int(stop/step)):
            t = step*i
            k = slope

            if mode == "ramp":
                answer = k * ramp_response(t)
            elif mode == "impulse":
                answer = impulse_response(t)
            elif mode == "step":
                answer = k * ramp_response(t) \
                    - k * ramp_response(t-tr) * u(t-tr)
            elif mode == "ideal_step":
                answer = size * step_response(t)
            elif mode == "pulse":
                answer = answer3 = k * ramp_response(t) \
                    - k * ramp_response(t-tr) * u(t-tr) \
                    - k * ramp_response(t-tr-width) * u(t-tr-width) \
                    + k * ramp_response(t-2*tr-width) * u(t-2*tr-width)
            else:
                assert 0, "UNRECOGNIZED MODE %s" % mode

            # Add the (time, answer) pair... NOTE: we must
            # extract the real value from the answer
            result.append((t, answer.real))

        # Add the results to the most recent plot
        g = self.__Plotter
        g("set output \"data.ps\"")
        g("set terminal postscript color")
        g("set terminal postscript solid")

        d = Gnuplot.Data(result,
                         title='AWE(%s/%s)' % (L,M),
                         with='lines')
        g.replot(d)            

    def __tran_normal(self, stop, step, name, method):
        # stop : time to stop the analysis
        # step : time step interval
        last_solution = numpy.zeros((self.Size,1))

        if method == 1:
            method = self.Trapezoidal
        else:
            method = self.BackwardEuler
        
        for i in range(0, int(stop/step)):
            #print "tran t=%s" % self.Time
            if method == self.BackwardEuler:
                self.Time = i*step
                mna_total = self.G + (self.C / step)
                intermediate = numpy.dot( (self.C / step), last_solution )
                rhs_total = intermediate + self.B
                last_solution = numpy.linalg.solve_linear_equations(\
                    mna_total, rhs_total)
            elif method == self.Trapezoidal:
                self.Time = i*step
                mna_total = self.G + (2*self.C / step)
                intermediate = numpy.dot( ((2*self.C) / step - self.G), 
                                          last_solution )
                rhs_total = intermediate + self.B
                self.Time = (i+1)*step
                rhs_total += self.B
                last_solution = numpy.linalg.solve_linear_equations(\
                    mna_total, rhs_total)
            self.__Results[name].append((self.Time, last_solution))

    #######################################
    ##### RESULTS RETRIEVAL SECTION #######
    #######################################
    def IDC(self, node):
        # convenience routine for testing
        return self.retrieve_raw( self.DC, self.Current, node )

    def VDC(self, node):
        # convenience routine for testing
        return self.retrieve_raw( self.DC, self.Voltage, node )

    def VF(self, node):
        # convenience routine for testing
        return self.retrieve_raw( self.AC, self.Voltage, node )

    def IF(self, node):
        # convenience routine for testing
        return self.retrieve_raw( self.AC, self.Current, node )

    def retrieve_raw(self, analysis, plot_type, node):
        if analysis in [self.AC, self.Transient]:
            results = self.__Results[analysis]
            # AC/Tran Voltage: returns the (x,y) pair
            if plot_type == self.Condition:
                return self.Cond
            else:
                if plot_type == self.Voltage:
                    index = self.TranslateNode(node)
                elif plot_type == self.Current:
                    index = self.TranslateBranch(node)
                else:
                    assert 0, plot_type

                vec = []
                for freq, solution in results:
                    vec.append((freq, solution[index-1]))
                if len(vec) == 1:
                    #print "1", vec
                    #print "2", vec[0]
                    #print "3", vec[0][1]
                    return vec[0][1][0]
                else:
                    return vec

        elif analysis == self.DC:
            # DC voltage/current: returns the y value
            dc_result = self.__Results[self.DC]
            if plot_type == self.Voltage:
                index = self.TranslateNode(node)
                return dc_result[index-1]
            elif plot_type == self.Current:
                index = self.TranslateBranch(node)
                return dc_result[index-1]

        else:
            assert 0, "Unsupported analysis %s" % analysis

    def vdb(self, num):
        return 20*numpy.log10(abs(num))

    def phase(self, num):
        x = numpy.real(num)
        y = numpy.imag(num)
        r = numpy.sqrt(x**2+y**2)
        phi = numpy.arccos(x/r)
        if y < 0.:
            phi = 2.*numpy.pi-phi
        # FIXME: I added this extra -numpy.pi in order to 
        # match spice... in addition to the source inversion
        return phi - numpy.pi       

    def cartesianToPolar(self, num):
        x = numpy.real(num)
        y = numpy.imag(num)
        r = numpy.sqrt(x**2+y**2)
        phi = numpy.arccos(x/r)
        if y < 0.:
            phi = 2.*numpy.pi-phi
        return r, phi - numpy.pi

    def retrieve_signal(self, analType, funcName, signalSpec):
        spice_func_dict = {
            "v" : None, 
            "i" : None,
            "cond" : None,
            "vm" : self.XYFunc(abs),
            "vp" : self.XYFunc(self.phase),
            "vdb" : self.XYFunc(self.vdb),
            "vr" : self.XYFunc(numpy.real),
            "vi" : self.XYFunc(numpy.imag),
            "im" : self.XYFunc(abs),
            "ip" : self.XYFunc(self.phase),
            "idb" : self.XYFunc(self.vdb),
            "ir" : self.XYFunc(numpy.real),
            "ii" : self.XYFunc(numpy.imag),
            }
        if funcName not in spice_func_dict:
            msg = "Unknown signal retrieval function: %s" % funcName
            raise CircuitException(msg)
        if funcName == "cond":
            plot_type = self.Condition
        elif funcName[0] == "i":
            plot_type = self.Current
        elif funcName[0] == "v":
            plot_type = self.Voltage

        signals = signalSpec.split(',')
        assert len(signals) == 1 or len(signals) == 2
        if len(signals) == 1:
            return map(spice_func_dict[funcName], 
                       self.retrieve_raw(analType, plot_type, signals[0]) )
        else:
            assert 0, "Not supported"

    def plot(self, analysis, funcName, signalSpec, logx, logy, file):
        self.reset_plotter()
        self.__plot(analysis, funcName, signalSpec,
                    logx, logy, file)

    def replot(self, analysis, funcName, signalSpec, logx, logy, file):
        self.__plot(analysis, funcName, signalSpec,
                    logx, logy, file)

    def __plot(self, analysis, funcName, signalSpec, logx, logy, file):
        #print "plot", analysis, funcName, signalSpec
        #print "__plot", funcName, signalSpec
        if analysis not in self.__Results:
            raise CircuitException(\
                "Analyis does not exist '%s': %s" % (analysis, 
                                                     self.GetAnalyses()) )
        if analysis in [self.AC, self.Transient]:
            results = self.__Results[analysis]
            if len(results) == 1:
                # no need to plot anything
                freq, result = results[0]
                #print "***", result
                # FIXME: cartesianToPolar for transient results
                for node_name in self.GetNodes():
                    index = self.TranslateNode(node_name)
                    print "V(%s) = %s" % (node_name, 
                                          cartesianToPolar(result[index-1]))
                for branch_name in self.GetBranches():
                    index = self.TranslateBranch(branch_name)
                    print "I(%s) = %s" % (branch_name,
                                          cartesianToPolar(result[index-1]))
                return
            #g = Gnuplot.Gnuplot()
            g = self.__Plotter
            if logx:
                g("set logscale x")
            if logy:
                g("set logscale y")
            if file is not None:
                g("set output \"%s\"" % file)
                g("set terminal postscript eps")
            else:
                g("set terminal X11")

            g("set grid")
            if analysis == self.AC:
                g.title('Frequency Response')
                g.xlabel('Frequency')
                g.ylabel('Magnitude')
            elif analysis == self.Transient:
                g.title('Transient Response')
                g.xlabel('Time')
                g.ylabel('Voltage')

            vec = self.retrieve_signal( analysis, funcName, signalSpec )
            #print "***vec", vec
            d = Gnuplot.Data(vec,
                             title="%s(%s)" % (funcName, signalSpec),
                             with='lines')
            g.replot(d)

        elif analysis == self.DC:
            dc_result = self.__Results[self.DC]
            if len(args) == 0:
                # plot everything
                #print "NODES:", x.GetNodes()
                for node_name in self.GetNodes():
                    index = self.TranslateNode(node_name)
                    print "V(%s) = %s" % (node_name, dc_result[index-1])
                for branch_name in self.GetBranches():
                    index = self.TranslateBranch(branch_name)
                    print "I(%s) = %s" % (branch_name, dc_result[index-1])
            else:
                for plot_type, node in args:
                    if plot_type == self.Voltage:
                        index = self.TranslateNode(node)
                        print "V(%s) = %s" % (node, dc_result[index-1])
                    elif plot_type == self.Current:
                        index = self.TranslateBranch(node)
                        print "I(%s) = %s" % (node, dc_result[index-1])
                #print dc_result

    def GetAnalyses(self):
        """Returns the name of all of the available analysis results"""
        return self.__Results.keys()

    def GetBranches(self):
        """Returns the names of all the defined branches"""
        return self.__Branches.keys()

    def GetNodes(self):
        return self.__NodeMap.keys()

    def Print(self, name=None):
        # Print out everything to the screen.
        if name is None:
            print "G=", self.G
            print "C=", self.C
            print "X=", self.X
            print "B=", self.B
            print "Branches=", self.Branches
            print "Nodemap=", self.__NodeMap
        elif name == "X":
            print self.X
        elif name == "Branches":
            print self.__Branches

    ##############################
    ### BUILDING THE CIRCUIT #####
    ##############################            
    def TranslateBranch(self, branchName):
        if branchName not in self.__Branches:
            assert 0, "Unknown branch '%s': %s" % (branchName, 
                                                   self.GetBranches())
        return self.__Branches[branchName]

    def TranslateNode(self, nodeName):
        # Return the translation for nodeName
        if nodeName == 0:
            return 0
        if nodeName not in self.__NodeMap:
            raise CircuitException(\
                "Unknown node '%s': %s" % (nodeName, self.GetNodes()) )
        return self.__NodeMap[nodeName]

    def __MapNode(self, nodeName):
        # Returns the internal node index for nodeName.
        # If nodeName does not exist, it is created and
        # added to the map.
        if nodeName == 0 or nodeName == '0':
            return nodeName
        if nodeName not in self.__NodeMap:
            self.__Size += 1
            self.__NodeMap[nodeName] = self.__Size
            self.__X.append((self.Voltage, self.__Size))
        return self.__NodeMap[nodeName]

    def __AddAdmittance(self, node1, node2, x, y):
        # node1: TERMINAL 1 external node name
        # node2: TERMINAL 2 external node name
        # x : dictionary to add equations to
        # y : magnitude of the Admittance
        assert node1 != node2
        n1 = self.__MapNode(node1)
        n2 = self.__MapNode(node2)
        if n1 != 0:
            x.setdefault((n1,n1),0.0)
            x[(n1,n1)] += y
        if n2 != 0:
            x.setdefault((n2,n2),0.0)
            x[(n2,n2)] += y
        if n1 != 0 and n2 != 0:
            x.setdefault((n1,n2),0.0)
            x[(n1,n2)] -= y
            x.setdefault((n2,n1),0.0)
            x[(n2,n1)] -= y

    def AddResistor(self, name, node1, node2, r):
        assert r > 0
        self.__AddAdmittance(node1, node2, self.__G, 1.0/r)
        
    def AddCapacitor(self, name, node1, node2, c):
        assert c > 0
        self.__AddAdmittance(node1, node2, self.__C, c)

    def AddInductor(self, name, node1, node2, l):
        assert l > 0
        self.__AddBranch(name, node1, node2)
        self.__C[(self.__Size, self.__Size)] = -l

    def AddCCVS(self, name, node1, node2, refDevice, A):
        # A : transconductance
        # refDevice: name of the device to measure current through
        # FIXME: Should these be +=?
        n1 = self.__MapNode(node1)
        n2 = self.__MapNode(node2)
        index = self.__Branches[refDevice]
        self.__AddBranch( name, node1, node2 )
        self.__G[(self.__Size, index)] = -1.0
        if n1 != 0:
            self.__G[(n1, self.__Size)] = 1.0
            self.__G[(self.__Size, n1)] = 1.0/A
        if n2 != 0:
            self.__G[(n2, self.__Size)] = -1.0
            self.__G[(self.__Size, n2)] = -1.0/A

    def AddVCCS(self, name, refNode1, refNode2, node1, node2, A):
        # A : gain
        # refNode1, refNode2: external reference nodes
        # (node1 - node2) <= A*(refNode1 - refNode2)
        ref1 = self.__MapNode(refNode1)
        ref2 = self.__MapNode(refNode2)
        n1 = self.__MapNode(node1)
        n2 = self.__MapNode(node2)
        #self.__AddBranch( name, node1, node2 )
        if ref1 != 0:
            if n1 != 0:
                self.__G.setdefault((n1, ref1), 0.0)
                self.__G[(n1, ref1)] += A
            if n2 != 0:
                self.__G.setdefault((n2, ref1), 0.0)
                self.__G[(n2, ref1)] += -A
        if ref2 != 0:
            if n1 != 0:
                self.__G.setdefault((n1, ref2), 0.0)
                self.__G[(n1, ref2)] += -A
            if n2 != 0:
                self.__G.setdefault((n2, ref2), 0.0)
                self.__G[(n2, ref2)] += A

    def AddVCVS(self, name, refNode1, refNode2, node1, node2, A):
        # A : gain
        # refNode1, refNode2: external reference nodes
        # (node1 - node2) <= A*(refNode1 - refNode2)
        ref1 = self.__MapNode(refNode1)
        ref2 = self.__MapNode(refNode2)
        self.__AddBranch( name, node1, node2, 1.0/A )
        if ref1 != 0:
            self.__G[(self.__Size, ref1)] = -1.0
        if ref2 != 0:
            self.__G[(self.__Size, ref2)] = 1.0
        
    def __AddBranch(self, name, node1, node2, mag=1.0):
        # name: label for the branch
        # node1: external node name
        # node2: external node name
        n1 = self.__MapNode(node1)
        n2 = self.__MapNode(node2)
        self.__Size += 1
        # FIXME: this is written badly
        if n1 == 0:
            self.__G[(n2, self.__Size)] = mag
            self.__G[(self.__Size, n2)] = mag
        elif n2 == 0:
            self.__G[(n1, self.__Size)] = mag
            self.__G[(self.__Size, n1)] = mag
        else:
            self.__G[(n1, self.__Size)] = mag
            self.__G[(self.__Size, n1)] = mag
            self.__G[(n2, self.__Size)] = -mag
            self.__G[(self.__Size, n2)] = -mag
        self.__X.append((self.Current, self.__Size))
        self.__Branches[name] = self.__Size

    def __ConstructSourceFunction(self, vdc, vac, pulse_args):
        # vdc: the value to be used during DC and Transient analyses
        # vac: the value to be used during AC analyses
        # pulse: specification tuple for a pulse which consists of the 
        #        following  v1 : voltage before/after pulse
        #                   v2 : voltage during pulse
        #                   td : delay before pulse starts
        #                   tr : rise time
        #                   tf : fall time
        #                   pw : pulse width
        #                   per : period
        # according to nghelp: if a pulse() is specified, that is
        #   the value to be used at time t=0
        def pulse(time, v1, v2, td, tr, tf, pw, per):
            t = time % per
            t1, t2, t3, t4 = td, td+tr, td+tr+pw, td+tr+pw+tf
            if t < t1 or t >= t4:
                # outside the pulse
                val = v1
            if t >= t1 and t < t2:
                # rising edge
                slope = (v2-v1)/(t2-t1)
                val = v1 + slope*(t - t1)
            elif t >=t2 and t < t3:
                # inside the pulse
                val = v2
            elif t >=t3 and t < t4:
                # falling edge
                slope = (v1-v2)/(t4-t3)
                val = v2 + slope*(t-t3)
            return val

        def source(time, analysis, vac, vdc): 
            return {SpiceCircuit.AC : vac,
                    SpiceCircuit.DC : vdc,
                    SpiceCircuit.Transient : vdc}[analysis](time)

        if pulse_args is not None:
            time_func = lambda t: pulse(t, *pulse_args)
        else:
            time_func = lambda t: vdc

        ac_func = lambda t: vac
        source_func = lambda t, anal: source(t, anal, ac_func, time_func)
        return source_func

    def AddVoltageSource(self, name, node1, node2, vdc=0.0, 
                         vac=0.0, pulse_args=None):
        self.__AddBranch(name, node1, node2)
        source_func = self.__ConstructSourceFunction(vdc, vac, pulse_args)
        self.__B[self.__Size] = source_func

    def AddCurrentSource(self, name, node1, node2, iac=0.0, 
                         idc=0.0, pulse_args=None):
        # FIXME: This does not handle not being connect to ground
        n1 = self.__MapNode(node1)
        n2 = self.__MapNode(node2)
	if n1 != 0:
	    source_func = self.__ConstructSourceFunction(idc, iac, pulse_args)
	    self.__B[n1] = source_func
	if n2 != 0:
	    v1, v2, td, tr, tf, pw, per = pulse_args
	    new_args = (-v1, -v2, td, tr, tf, pw, per)
	    source_func = self.__ConstructSourceFunction(-idc, -iac, new_args)
	    self.__B[n2] = source_func
