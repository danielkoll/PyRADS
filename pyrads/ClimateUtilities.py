from __future__ import division, print_function, absolute_import
### THIS SCRIPT IS TAKEN FROM THE COURSEWARE OF
###  Pierrehumbert, 2010, Principles of Planetary Climate
###
###  Downloaded from:
###  https://geosci.uchicago.edu/~rtp1/PrinciplesPlanetaryClimate/Courseware/coursewarePortal.html

#ToDo:  *Check for right column length in setitem and addCurve
#
#       *Implement show,hide for curves (in X() and Y())
#       *Implement missing data coding .
#            self.setMissingDataCode(code) (char or numeric)
#            possibly translate, check and force consistency
#       *Possibly handle missing data for plotting using a
#        fill() or interp() method to fill in missing data using
#        the interpolation routine. The best way to handle
#        missing data in computations is by using Masked Arrays
#
#       *Make and import a PathFinder module, which the instructor
#        will customize for the site.  This module will help the
#        students find locations of datasets and chapter scripts.
#        Could provide "links" to find script that produced a given
#        figure in the text.  Alternately, we could just
#        define directory strings like WorkbookDatasets here in
#        ClimateUtilities.
#
#       *Note: Other courseware modules should avoid
#        importing numpy or Numeric directly. Get it
#        by using "from ClimateUtilities import *" so
#        that the preferred version is imported automatically
#
#       *Looking ahead to Python 3, all print statements should
#        be changed to some kind of "message" function, which
#        can invoke either python 2.xx "print" or python 3 print(...)
#
#

#-----------------------------------------------
#
#Import array package
#
#First try to import numpy, then fall back
#on Numeric.  Arrange things so that the array
#package can be interchangeably referenced as
#either numpy.* or Numeric.* , though it is really
#numpy that is fully supported.  (Almost everything
#will continue to work with Numeric, though, which
#may be useful for older installations)
#------------------------------------------------

try:
    import numpy
    import numpy as Numeric   #For backwards compatibility
    numpy.Float = numpy.float #For backwards compatibility
except:
    try:
        import Numeric
        import Numeric as numpy    #For frontwards compatibility
        numpy.float = numpy.Float  #For frontwards compatibility
        print( "numpy not found. Using Numeric instead")
        print( "Everything should still work, but consider upgrading to numpy")
    except:
        print( "Neither numpy nor Numeric found.")
        print( "Please install numpy (preferred) or Numeric.")

#-----------------------------------------------------
#Import graphics utilities
#
#First try to import ClimateGraphics, which contains the
#implementation of the plotting commands.  The default version
#distributed with the courseware uses PyNgl as a driver, but
#ClimateGraphics can be fairly easily localized to use other
#graphics drivers (e.g. MatPlotLib).  If the import fails, a
#dummy graphics stub module is imported, which allows the courseware
#to be run without the user needing to explicitly comment out
#the plot calls in the Chapter Scripts.
#------------------------------------------------------

try:
    from .ClimateGraphicsMPL import * #Try importing MatPlotLib
#    from ClimateGraphics import * #Try importing Ngl driver
except:
    # If Ngl not found, try importing the graphics interface
    # that uses MatPlotLib.  If you have both installed and
    #prefer MatPlotLib you can change the order of imports here
    #or just do "from ClimateGraphicsMPL import *" after you
    #import ClimateUtilities .
    try:
        from .ClimateGraphicsMPL import * #Try importing MatPlotLib
    except:
        print( "  ")
        print( "Graphics not implemented.")
        print( "Plot routines will not produce graphs.")
        print( "Instead, you can save results from a Curve")
        print( "object c into a text file using c.dump(<FILENAME>)")
        print( "and then plot the data using the graphics program")
        print( "of your choice.")
        from .DummyGraphics import *

#Section 1: -----Data handling utilities-------------------------------

# ToDo: Put documentation on use of Curve object here!
#
#       Add keyword arguments for axes, etc,
#
#       How to handle missing data codes? Note that
#         since lists are converted to Numeric arrays,
#         text codes can't be used.
#       Provide a method to set which data is X.
#
#       Add output option to put out a LaTeX formatted table
#
#       Provide an easy way to add data column long-names
#
class Curve:
    def __init__(self):
         self.Xid = None # id of data to be considered as X
         self.data = {} # Dictionary of data columns
         self.label = {} #Data column label dictionary
         self.scatter = {} #Marker for making a curve a scatter plot
                           #i.e. suppress line drawing and plot symbol
         self.idList = [] #Keeps track of original order of data
         self.NumCurves = 0
         self.description = None # A string providing general information
         self.PlotTitle = '' # Title of the plot
         self.switchXY = 0 # Switches X and Y axes for plotting
         self.reverseX = 0 # Reverses X axis for plotting
         self.reverseY = 0 # Reverses Y axis for plotting
         self.XlogAxis = 0 # Use logarithmic axis for X
         self.YlogAxis = 0 # Use logarithmic axis for Y
         self.Xlabel = '' #X axis label
         self.Ylabel = '' #Y axis label
         #
         # Colors and line titles
         #
    #
    #Install a new curve in the data set, optionally with a variable name and label
    #Any 1D indexable object can be installed here:a Numeric array, a Masked Array,
    #a Masked Variable (as in cdms) or an ordinary list.
    def addCurve(self,data,id = '',label = ''):
        self.NumCurves += 1
        if len(id) == 0:
            id = 'v%d'%(self.NumCurves-1)
        #Transform data from list to a Numeric array here
        if type(data) == type([]):
            data = Numeric.array(data)
        self.data[id] = data
        if self.Xid == None:
            self.Xid = id   #Sets the default id for the X variable
        self.idList.append(id) #Keep track of order in which columns added
        self.label[id] = label
        self.scatter[id] = False
        #ToDo: Add checking for consistent lengths, type, etc. of what's being added
    def listVariables(self):
        return self.idList
    def __getitem__(self,id):
        return self.data[id]
    def __setitem__(self,id,data):
        try:
            n = len(data[:])
        except:
            print( "Object on RHS is not indexable")
            return None
        #Transform data from list to a Numeric array here
        if type(data) == type([]):
            data = Numeric.array(data)
        if id in self.data.keys():
            self.data[id] = data
        else:
            self.addCurve(data,id)

    #Method to return Numeric abcissa array for plotting.
    def X(self):
        #Use of cross section lets us get data from any indexed object
        #However, since Masked variables and Masked Arrays yield
        #their same types as cross sections, we have to check
        #explicitly for a _data component.
        #
        #This method is used mainly for generating arrays used in plotting,
        #and should be streamlined at some point.
        temp = self.data[self.Xid][:]
        if hasattr(temp,'_data'):
            temp = temp._data[:]
        return Numeric.array(temp,Numeric.Float)

    #Method to return Numeric ordinate array for plotting
    def Y(self):
        #Use of cross section lets us get data from any indexed object
        outArray = []
        for id in self.idList:
            if not (id == self.Xid):
              column = self.data[id]
              if hasattr(column,'_data'):
                outArray.append(column._data[:]) #Deals with masked arrays and variables
              else:
                outArray.append(column[:])
        return Numeric.array(outArray,Numeric.Float)

    #Dumps curve to a tab-delimited ascii file with column header
    def dump(self,fileName = 'out.txt'):
        outfile = open(fileName,'w')
        # Write out the data description if it is available.
        if not (self.description == None):
            if not self.description[-1] == '\n':
                self.description += '\n' #Put in a newline if needed
            outfile.write(self.description)
        header = ""
        fmt = ""
        ids = self.idList
        for id in ids:
            header += id+'\t'
            fmt += '%e\t'
        header = header[0:-1]+'\n'  # Replaces last tab with a newline.
        fmt = fmt[0:-1] + '\n'
        outfile.write(header)
        for i in range(len(self.data[ids[0]])):
            out = tuple([(self.data[id])[i] for id in ids])
            outfile.write(fmt%out)
        outfile.close()

    #Extracts a subset of the data and returns it as a new Curve.
    #This is useful if you only want to plot some of the columns.
    #The input argument dataList is a list of column names
    def extract(self,dataList):
        c = Curve()
        for dataName in dataList:
            c.addCurve(self[dataName],dataName)
        return c





#from string import atof
#  deprecated in Python 2.7, gone in Python 3... just use float()

# Scans a list of lines, locates data lines
# and size, splits of column headers and
# splits off general information text.
#
#A header can optionally be specified as input.
#The delimiter need not be specified if the file
#uses any whitespace character (including tabs) as
#column delimiters. Commas are not whitespace characters,
#and so need to be specified if they are used.
#
#ToDo:
#       *Implement missing data coding (with default '-').
#            self.setMissingDataCode(code) (char or numeric)
#            possibly translate, check and force consistency
#       *Replace optional positional arguments with keyword arguments
def scan(buff,inHeader=None,delimiter = None):
    if inHeader == None:
        inHeader = []
    #First delete blank lines
    buff = clean(buff)
    #
    #Now look for patterns that indicate data lines
    #
    startDataLine,endDataLine = findData(buff)
    # Found number of items. Now read in the data
    #
    #Read in the first line. Is it a header?
    header = []
    if delimiter == None:
        line = buff[startDataLine].split()
    else:
        line = buff[startDataLine].split(delimiter)
    #
    try:
        #atof(line[0])
        float(line[0])
    except:
        header = line
    if len(header) == 0:
        header = ['V%d'%i for i in range(len(line))]
        istart = 0
    else:
        istart = 1
    # Replace with inHeader if inHeader has been specified
    #  (Only use the input header if its length is consistent)
    if len(inHeader) == len(line):
        header = inHeader
    #
    # Read in the rest of the lines
    #
    varlist = [[] for i in range(len(header))]
    for line in buff[(startDataLine+istart):endDataLine]:
        if delimiter == None:
            items = line.split()
        else:
            items = line.split(delimiter)
        try:
         for i in range(len(varlist)):
            #varlist[i].append(atof(items[i]))
            varlist[i].append(float(items[i]))
        except:
            print( items)
    vardict = {}
    for name in header:
        vardict[name] = varlist[header.index(name)]
    return vardict,header #header is returned so we can keep cols in orig order

#Eliminates blank lines
def clean(buff):
    buff = [line.strip() for line in buff]
    while(1):
      try:
         buff.remove('')
      except:
         break
    return buff

def findData(buff):
    runStarts = []
    for i in range(len(buff)-1):
        dn = abs(len(buff[i].split())-len(buff[i+1].split()))
        if not dn==0:
             runStarts.append(i+1)
    runStarts.append(len(buff))
    #Find index of run with max length
    nmax = -1
    for i in range(len(runStarts)-1):
        n = runStarts[i+1]-runStarts[i]
        if n > nmax:
            nmax = n
            imax = runStarts[i]
    #Deal with case where entire file is one run
    if nmax == -1:
        nmax = len(buff)
        imax = 0
    return imax,imax+nmax


#
# Function to read space or tab-delimited file into a curve object
# The input header is a list of names of the variables in each
# columns, which can be input optionally, mainly to deal with
# the case in which this information is not in the file being read
import string
def readTable(filename,inHeader = None,delimiter = None):
    f = open(filename)
    buff = f.readlines()
    data,header = scan(buff,inHeader,delimiter)
    c = Curve()
    for key in header:
        c.addCurve(data[key],key)
    return c




#==============================================
#---Section 2: Math utilities-------------------------------------------
#==============================================

#A dummy class useful for passing parameters.
#Just make an instance, then add new members
#at will.
class Dummy:
    pass

#---Polynomial interpolation and extrapolation (adapted
#from Numerical Recipes.
#
#It's used in Romberg extrapolation, but could be useful
#for polynomial OLR fits and so forth as well. Also
#needs online documentation
def polint(xa,ya,x):
    n = len(xa)
    if not (len(xa) == len(ya)):
            print( "Input x and y arrays must be same length")
            return "Error"
        #Set up auxiliary arrays
    c = Numeric.zeros(n,Numeric.Float)
    d = Numeric.zeros(n,Numeric.Float)
    c[:] = ya[:]
    d[:] = ya[:]
    #Find closest table entry
    ns = 0
    diff = abs(xa[0]-x)
    for i in range(n):
            difft = abs(xa[i]-x)
            if difft < diff:
                diff = difft
                ns = i
    y=ya[ns]
    for m in range(1,n):
        for i in range(n-m):
            ho=xa[i]-x
            hp=xa[i+m]-x
            w=c[i+1]-d[i]
            c[i] = ho*w/(ho-hp)
            d[i] = hp*w/(ho-hp)
        if 2*ns < (n-m):
            dy = c[ns]
        else:
            ns -= 1
            dy = d[ns]
        y += dy
        #You can also return dy as an error estimate. Here
        #to keep things simple, we just return y.
    return y

#---------------------------------------------------------
#         Class for doing polynomial interpolation
#         from a table, using polint
#---------------------------------------------------------
#
#Usage:
#         Let xa be a list of independent variable
#         values and ya be a list of the corresponding
#         dependent variable values. Then, to create a function
#         f (actually a callable object, techically) that interpolates
#         or extrapolates to any value x, create f using
#                     f = interp(xa,ya)
#         Then you can get the value you want by invoking f(x)
#         for your desired x.
#
#         By default, the interpolator does fourth-order interpolation
#         using the four nearest neighbors. You can change this by
#         using an optional third argument to the creator. For
#         example
#
#                     f = interp(xa,ya,8)
#         will use the 8 nearest neighbors (if they are available)
#------------------------------------------------------------
class interp:
    def __init__(self,xa,ya,n=4):
        self.xa = Numeric.array(xa)
        self.ya = Numeric.array(ya)
        self.n = n
    def __call__(self,x):
        #Find the closes index to x
        if self.xa[0] < self.xa[-1]:
            i = Numeric.searchsorted(self.xa,x)
        else:
            i = Numeric.searchsorted(-self.xa,-x)
        i1 = max(i-self.n,0)
        i2 = min(i+self.n,len(self.xa))
        return polint(self.xa[i1:i2],self.ya[i1:i2],x)


#----Quadrature (definite integral) by Romberg extrapolation.
#**ToDo: Add documentation and help string

#Before developing a general quadrature class, we'll
#implement a class which efficiently carries out trapezoidal rule
#integration with iterative refinement
class BetterTrap:
    def __init__(self,f,params,interval,nstart):
        self.f = f
        self.n = nstart
        self.interval = interval
        self.params = params
        self.integral = self.dumbTrap(nstart)
    def dumbTrap(self,n):
        a = self.interval[0]
        b = self.interval[1]
        dx = (b-a)/n
        sum = dx*(self.f(a,self.params)+self.f(b,self.params))/2.
        for i in range(1,n):
            x = a+i*dx
            sum = sum + self.f(x,self.params)*dx
        return sum
    def refine(self):
        #Compute the sum of f(x) at the
        #midpoints between the existing intervals.
        #To get the refinement of the trapezoidal
        #rule sum we just add this to half the
        #previous result
        sum = 0.
        a = self.interval[0]
        b = self.interval[1]
        dx = (b-a)/self.n
        #Remember: n is the number of subintervals,
        #not the number of endpoints. Therefore we
        #have one midpoint per subinterval. Keeping that
        #in mind helps us get the range of i right in
        #the following loop
        for i in range(self.n):
            sum = sum + self.f(a+(i+.5)*dx,self.params)*(dx/2.)
        #The old trapezoidal sum was multiplied by
        #the old dx.  To get its correct contribution
        #to the refined sum, we must multiply it by .5,
        #because the new dx is half the old dx
        self.integral = .5*self.integral + sum
        #
        #Update the number of intervals
        self.n = 2*self.n


#Here I define a class called
#romberg, which assists in carrying out evaluation of
#integrals using romberg extrapolation. It assumes polint has
#been imported

class romberg:
    def __init__(self,f,nstart=4):
        self.nstart = nstart
        self.trap = None
        #
        #-------------------------------------------------
        #This snippit of code allows the user to leave the
        #parameter argument out of the definition of f
        #if it isn't needed
        #
        self.fin = f
        #Find the number of arguments of f and append a
        #parameter argument if there isn't any.
        nargs = f.func_code.co_argcount
        if nargs == 2:
            self.f = f
        elif nargs ==1:
            def f1(x,param):
                return self.fin(x)
            self.f = f1
        else:
            name = f.func_name
            print( 'Error: %s has wrong number of arguments'%name)
        #-----------------------------------------------------
        #
        #We keep lists of all our results, for doing
        #Romberg extrapolation. These are re-initialized
        #after each call
        self.nList = []
        self.integralList = []
    def refine(self):
        self.trap.refine()
        self.integralList.append(self.trap.integral)
        self.nList.append(self.trap.n)
        dx = [1./(n*n) for n in self.nList]
        return polint(dx,self.integralList,0.)
    #
    #Use a __call__ method to return the result. The
    #__call__ method takes the interval of integration
    #as its mandatory first argument,takes an optional
    #parameter argument as its second argument, and
    #an optional keyword argument specifying the accuracy
    #desired.
    #**ToDo: Introduce trick to allow parameter argument of
    #integrand to be optional, as in Integrator.  Also, make
    #tolerance into a keyword argument
    #
    def __call__(self,interval,params=None,tolerance=1.e-6):
        self.nList = []
        self.integralList = []
        #Make a trapezoidal rule integrator
        self.trap = BetterTrap(self.f,params,interval,self.nstart)
        self.nList.append(self.nstart)
        self.integralList.append(self.trap.integral)
        #
        #Refine initial evaluation until
        oldval = self.refine()
        newval = self.refine()
        while abs(oldval-newval)>tolerance:
            oldval,newval = newval,self.refine()
        return newval



#-------Runge-Kutta ODE integrator
#**ToDo:
#      * Implement a reset() method which resets to initial conditions.
#         Useful for doing problem over multiple times with different
#         parameters.
#
#      * Referring to the independent variable as 'x' is awful, and
#        confusing in many contexts.  Introduce a variable name
#        dictonary with default names like 'independent' and 'dependent'
#        (and short synonyms) so if fi is the integrator object
#        fi['indep'] is the current value of the independent variable
#        and fi['dep'] is the current (vector) value of the dependent
#        variable. Then allow user to rename or synonym these to
#        the actual user-supplied names.  This is an alternative
#        to using the list returned by fi.next(). Then expunge
#        all references to things like fi.x from the examples and
#        chapter scripts. They are too confusing.  Similarly,
#        change the name of the increment from dx to something else
#
#      * Similarly, we could introduce a dictionary of some sort
#        to make it easier to set up multidimensional systems and
#        refer to the different vector components by name
#        (e.g. refer to v[0] at T, v[1] as dTdy , etc. )
#
#      * Make the integrator object callable. The call can return
#        a list of results for all the intermediate steps, or optionally
#        just the final value.
class integrator:
  '''
  Runge-Kutta integrator, for 1D or multidimensional problems

  Usage:

  First you define a function that returns the
  derivative(s), given the independent and dependent
  variables as arguments. The independent variable (think
  of time) is always a scalar. The dependent variable (think
  of the x and y coordinates of a planet) can be either a scalar
  or a 1D array, according to the problem. For the
  multidimensional case, this integrator will work with any
  kind of array that behaves like a Numeric array, i.e. supports
  arithmetic operations. It will not work with plain Python lists.
  The derivative function should return an array of derivatives, in
  the multidimensional case. The derivative function can have any name.

  The derivative function can optionally have a third argument, to provide
  for passing parameters (e.g. the radius of the Earth) to the
  function.  The "parameter" argument, if present, can be any Python
  entity whatsoever. If you need to pass multiple constants, or
  even tables or functions, to your derivative function, you can
  stuff them all into a list or a Python object.


  Example:
  In the 1D case, to solve the equation
                  dz/dt = -a*t*t/(1.+z*z)
  in which z is the dependent variable and t is the
  independent variable, your derivative function would  be
    def g(t,z):
       return -a*t*t/(1.+z*z)

  treating the parameter a as a global, or perhaps
    def g(t,z,params):
       return -params.a*t*t/(params.b+z*z)

  while in a 2D case, your function might look like:
    def g(t,z):
      return Numeric.array([-z[1],z[0]])

  or perhaps something like:
    def g(t,z):
      return t*Numeric.sin(z)

  or even
    def g(t,z,params):
      return Numeric.matrixmultiply(params(t),z)

  where params(t) in this case is a function returning
  a Numeric square matrix of the right dimension to multiply z.

  BIG WARNING:  Note that all the examples which return a
  Numeric array return a NEW INSTANCE (i.e. copy) of an
  array.  If you try to set up a global array and re-use
  it to return your results from g, you will really be
  just returning a REFERENCE to the same array each time,
  and each call will change the value of all the previous
  results. This will mess up the computation of intermediate
  results in the Runge-Kutta step. An example of the sort of thing
  that will NOT work is:
    zprime = Numeric.zeros(2,Numeric.Float)
    def g(t,z,params):
      zprime[0] = z[1]
      zprime[1] = -z[0]
      return zprime
  Try it out. This defines the harmonic oscillator, and a plot
  of the orbit should give a circle. However, it doesn't  The problem
  reference/value distinction.  The right way to define the function
  would be
    def g(t,z):
      return Numeric.array([z[1],-z[0]])
  Try this one. It should work properly now. Note that any arithmetic
  performed on Numeric array objects returns a new instance of an Array
  object. Hence, a function definition like
    def g(t,z):
      return t*z*z+1.
  will work fine.

  Once you have defined the derivitave function,
  you then proceed as follows.

  First c reate an integrator instance:
    int_g = integrator(g,0.,start,.01)

  where "0." in the argument list is the initial value
  for the independent variable, "start" is the initial
  value for the dependent variable, and ".01" is the
  step size. You then use the integrator as follows:

    int_g.setParams(myParams)
    while int_g.x < 500:
       print( int_g.next())

  The call to setParams is optional. Just use it if your
  function makes use of a parameter object. The next() method
  accepts the integration increment (e.g. dx) as an optional
  argument. This is in case you want to change the step size,
  which you can do at any time.  The integrator continues
  using the most recent step size it knows.

  Each call to int_g.next returns a list, the first of whose
  elements is the new value of the independent variable, and
  the second of whose elements is a scalar or array giving
  the value of the dependent variable(s) at the incremented
  independent variable.

  '''
  def __init__(self, derivs,xstart,ystart,dx=None):
    self.derivsin = derivs
    #
    #The following block checks to see if the derivs
    #function has a parameter argument specified, and
    #writes a new function with a dummy parameter argument
    #appended if necessary. This allows the user to leave
    #out the parameter argument from the function definition,
    #if it isn't needed.
    nargs = derivs.func_code.co_argcount
    if nargs == 3:
        self.derivs = derivs
    elif nargs == 2:
        def derivs1(x,y,param):
            return self.derivsin(x,y)
        self.derivs = derivs1
    else:
        name = derivs.func_name
        print('Error: %s has wrong number of arguments'%name)
    #
    #
    self.x = xstart
    #The next statement is a cheap trick to initialize
    #y with a copy of ystart, which works whether y is
    #a regular scalar or a Numeric array.
    self.y = 0.+ ystart
    self.dx = dx #Can instead be set with the first call to next()
    self.params = None
  #Sets the parameters for the integrator (optional).
  #The argument can be any Python entity at all. It is
  #up to the user to make sure the derivative function can
  #make use of it.
  def setParams(self,params):
      self.params = params
  #Computes next step.  Optionally, takes the increment
  #in the independent variable as an argument.  The
  #increment can be changed at any time, and the most
  #recently used value is remembered, as a default
  def next(self,dx = None):
     if not (dx == None):
         self.dx = dx
     h = self.dx
     hh=h*0.5;
     h6=h/6.0;
     xh=self.x+hh;
     dydx = self.derivs(self.x,self.y,self.params)
     yt = self.y+hh*dydx
     dyt = self.derivs(xh,yt,self.params)
     yt =self.y+hh*dyt
     dym = self.derivs(xh,yt,self.params)
     yt =self.y+h*dym
     dym += dyt
     dyt = self.derivs(self.x+h,yt,self.params)
     self.y += h6*(dydx+dyt+2.0*dym)
     self.x += h
     return self.x,self.y



#**ToDo:
#
#         Store the previous solution for use as the next guess(?)
#
#        Handle arithmetic exceptions in the iteration loop
#
class newtSolve:
    '''
    Newton method solver for function of 1 variable
    A class implementing Newton's method for solving f(x) = 0.

    Usage: solver = newtSolve(f), where f is a function with
    calling sequence f(x,params). Values of x such that
    f(x,params) = 0 are
    then found by invoking solver(guess), where guess
    is the initial guess.  The solver returns the string
    'No Convergence' if convergence fails. The argument
    params allows parameters to be passed to the function.
    It can be left out of the function definition if you don't
    need it. Note that params can be any Python object at all
    (including,e.g.,lists, functions or more complex user-defined objects)

    Optionally, one can specify the derivative function
    in the creator,e.g. solver = newtSolve(f,fp).
    If the derivative function isn't specified, the solver
    computes the derivative approximately using a centered
    difference. Note that in either case you can access
    the derivative function by invoking solver.deriv(x)
    As for f, fp can be optionally defined with a parameter
    argument if you need it. The same parameter object is
    passed to f and fp.

    Use solver.setParams(value) to set the parameter object
    Alternately, the parameter argument can be passed as
    an optional second argument in the solver call. (see
    example below).

    Adjustable constants:
     eps         Increment for computing numerical approximation to
                 the derivative of f
     tolerance   Accuracy criterion for ending the iteration
                 (an approximation to the error in the root)
     nmax        maximum number of iterations

    e.g. to change the maximum number of iterations for an instance
    of the class, set solver.nmax = 10 .
    ----------------Usage Examples-----------------------------

         Example 1: Function without parameters:
          def g(x):
              return x*x - 1.
          roots = newtSolve(g)
          roots(2.)

         Example 2, Function with parameters:
          def g(x,constants):
              return constants.a*x*x - constants.b
          roots = newtSolve(g)
          constants = Dummy()
          constants.a = 1.
          constants.b = 2.
          roots.setParam(constants)
          roots(2.)
          roots(1.)

         Example 2a:
         Instead of using roots.setParam(...) we could do
           roots(2.,constants)
           roots(1.)    the parameters are remembered
           constants.a = 3.
           roots(1.,constants)   We changed the constants

         Example 3, using scan to find initial guesses:
          def g(x):
              return x*x - 1.
          roots = newtSolve(g)
          guesses = roots.scan([-2.,2.],100)
          for guess in guesses:
              print( roots(guess))
    '''
    def __init__(self,f,fprime = None):
        self.fin = f
        #Find the number of arguments of f and append a
        #parameter argument if there isn't any.
        nargs = f.func_code.co_argcount
        if nargs == 2:
            self.f = f
        elif nargs ==1:
            def f1(x,param):
                return self.fin(x)
            self.f = f1
        else:
            name = f.func_name
            print( 'Error: %s has wrong number of arguments'%name)
        self.eps = 1.e-6
        def deriv(x,params):
            return (self.f(x+self.eps,params)- self.f(x-self.eps,params))/(2.*self.eps)
        if fprime == None:
            self.deriv = deriv
        else:
            #A derivative function was explicitly specified
            #Check if it has a parameter argument
            nargs = fprime.func_code.co_argcount
            if nargs == 2:
                self.deriv = fprime #Has a parameter argument
            elif nargs == 1:
                self.fprimein = fprime
                def fprime1(x,param):
                    return self.fprimein(x)
                self.deriv = fprime1
            else:
                name = fprime.func_name
                print( 'Error: %s has wrong number of arguments'%name)
        self.tolerance = 1.e-6
        self.nmax = 100
        self.params = None
    def __call__(self,xGuess,params = None):
        if not (params == None):
            self.setParams(params)
        x = xGuess
        for i in range(self.nmax):
            dx = (self.f(x,self.params)/self.deriv(x,self.params))
            x = x - dx
            if abs(dx) < self.tolerance:
                return x
        return 'No Convergence'
    def setParams(self,params):
        #**ToDo: Check if f1 has a parameter argument
        #defined, and complain if it doesn't
        self.params = params
    def scan(self,interval,n=10):
        #Finds initial guesses to roots in a specified
        #interval, subdivided into n subintervals.
        #e.g. if the instance is called "solver"
        #solver.scan([0.,1.],100) generates a list
        #of guesses between 0. and 1., with a resolution
        #of .01. The larger n is, the less is the chance that
        #a root will be missed, but the longer the search
        #will take.  If n isn't specified, the default value is 10
        #
        #ToDo: Replace this with a bisection search, allowing user
        #to specify the maximum number of distinct guesses that
        #need to be found.
        guessList = []
        dx = (interval[1]-interval[0])/(n-1)
        flast = self.f(interval[0],self.params)
        for x in [interval[0]+ i*dx for i in range(1,n)]:
            fnow = self.f(x,self.params)
            if ((fnow >= 0.)&(flast <=0.)) or ((fnow <= 0.)&(flast >=0.)):
                guessList.append(x)
            flast = fnow
        return guessList
