from __future__ import division, print_function, absolute_import
#----------Section 3: Plotting utilities-------------------------------
#Graphics interface customized for MatPlotLib
#     --------------->UNDER DEVELOPMENT
#
#This is imported into ClimateUtilities. If you want to use a
#different graphics package (e.g. MatPlotLib) as a graphics driver
#in place of Ngl, you only need to rewrite this module. The most
#important thing to rewrite is the plot(...) function, which takes
#a Curve object as input and produces a line plot.  The plot(...)
#function returns a plotObj object to allow further manipulation
#(mostly saving the plot in various formats), but if you have a plotting
#package that provides other ways of saving the plots, or if you just
#want to save the plots by doing screen dumps, you can just have
#the plot(...) function return None, or some dummy object.
#
#The other main function to customize is contour(...) which produces
#contour plots of an array. This function is not very extensively
#used in the courseware, so it is a lower priority for modification.
#This routine has not yet been implemented for MatPlotLib
#-----------------------------------------------------------------------


#Try to import matplotlib plotting routines
try:
    import matplotlib as mpl
    #Configure the backend so things work properly
    #with the Python intepreter and idle. Not needed if you are
    #using ipython. To work properly with idle, idle needs
    #to be started up with the command "idle =n".  Note that
    #this trick limits the options for saving a file; only png
    #format is supported. If you want higher resolution formats
    #(notably eps) you should run using ipython -pylab , and eliminate
    #the following two commands.
    mpl.rcParams['backend'] = 'TkAgg'
    mpl.rcParams['interactive'] = True
    #
    import pylab as pl
except:
    print( 'matplotlib not found on your system')
    print( 'You can still run the courseware, but')
    print( 'will not be able to plot from inside Python')

#A dummy class useful for passing parameters.
#Just make an instance, then add new members
#at will. Not sure why it also has to be defined
#here, since it is also defined in ClimateUtilities,
#but it seems to be necessary
class Dummy:
    pass

# A little class to make resources for Ngl plot options
class resource:
    def __init__(self):
        self.trYReverse = False
        self.trXReverse = False
        self.trXLog = False
        self.trYLog = False
        self.nglFrame = False
        self.nglDraw = False
        #ToDo: Missing data code resource, line styles, line colors

# A little class for use as a return object from plot(), so
# the user has an easy way to delete a window or save a plot
# to an eps file. Plot objects not implemented yet for MPL
# Is there a way to save plots non-interactively in MPL?
class plotObj:
    def __init__(self,workstation,plot,WorkstationResources=None):
        self.workstation = workstation
        self.plot = plot
        self.WorkstationResources = WorkstationResources
    #Deletes a plot window
    def delete(self):
        print( "In MatPlotLib just click the goaway button to dispose a plot")
    def save(self,filename = 'plot'):
        #**ToDo: Can we implement non-interactive save in MPL?
        print( "In MatPlotLib save plots interactively using the file button")
        print( " in the plot window")

#ToDo:
#       *Implement use of missing data coding.
#       *Provide some way to re-use window (e.g. by
#         returning w, and having it be an optional
#         argument. How to clear a window for re-use?
#       *Make some use of dashed lines in line styles
#       *The behavior of axis options like reverseX and
#        XlogAxis is confusing when switchXY = True. We
#        need to think about the semantics of what we
#        mean by the "X" axis, and stick to it.  As
#        currently implemented the "X" axis means the
#        horizontal axis for these options.  (However,
#        the semantics is different for the labeling options!)
#        If we change this (and we probably should), the
#        scripts plotting soundings will also need to be changed.
#        as well as the Workbook intructions and problem sets.
#

#List of line colors and line styles to use
lineColors = ['b','g','r','c','m','y','k']
lineStyles = ['-','--','-.',':']
lineThickness = [2,3,4]
plotSymbols = ['.','o','v','<','s','*','+','x']
def plot(c):
    fig = pl.figure() #Always start a new figure
    #**ToDo: Make sure log-axis options so they work consistently
    #        with Ngl implementation when x/y axes are switched.
    #        (Looks OK, but is that implementation the logical way?)
    #r.trXLog = c.XlogAxis
    #r.trYLog = c.YlogAxis
    if c.XlogAxis:
        pl.semilogx()
    if c.YlogAxis:
        pl.semilogy()
    if c.XlogAxis & c.YlogAxis:
        pl.loglog()
    #
    #Line styles (Not needed in MPL, handled automatically)
    #r.xyLineColors = lineColors
    #r.xyLineThicknesses = lineThickness
    #Plot title
    #r.tiMainString = c.PlotTitle
    pl.title(c.PlotTitle)
    #Axis labels (ToDo: add defaults)
    #X and Y axis labels
    #r.tiXAxisString = c.Xlabel
    #r.tiYAxisString = c.Ylabel
    if c.switchXY:
        pl.ylabel(c.Xlabel)
        pl.xlabel(c.Ylabel)
    else:
        pl.xlabel(c.Xlabel)
        pl.ylabel(c.Ylabel)

    #  Legends, for multicurve plot
    legends = []
    for id in c.listVariables():
        if not id == c.Xid:
            if len(c.label[id]) > 0:
                legends.append(c.label[id])
            else:
                legends.append(id)
    #ToDo: Add option to skip legends
    #
    #Suppress line drawing and just plot symbol for scatter plot curves
    #**ToDo: Implement for MatPlotLib
    #r.xyMarkers = plotSymbols
    #r.xyMarkerColors = lineColors
    #r.xyMarkerSizeF = .01
    #r.xyMarkLineModes = []
    formatList = []
    count = 0
    for id in c.listVariables():
        if not id == c.Xid:
            if c.scatter[id]:
                #r.xyMarkLineModes.append('Markers')
                color = lineColors[count%len(lineColors)]
                symbol = plotSymbols[count%len(plotSymbols)]
                formatList.append(color+symbol)
            else:
                color = lineColors[count%len(lineColors)]
                style = lineStyles[count%len(lineStyles)]
                formatList.append(color+style)
            count += 1
    #
    plotList = [] #Mainly so we can add legends. Could be done with label = ...
    if c.switchXY:
        count = 0
        for data in c.Y():
            plotList.append(pl.plot(data,c.X(),formatList[count]))
            count += 1
    else:
        count = 0
        for data in c.Y():
            plotList.append(pl.plot(c.X(),data,formatList[count]))
            count += 1
    #Do the legends
    pl.legend(plotList,legends)
    #
    #Do the axis reversal
    #We do it here, since we don't know the axis limits until
    #plotting is done
    #r.trYReverse = c.reverseY
    #r.trXReverse = c.reverseX
    axes = pl.gca() #Gets current axis
    if c.reverseX:
        axes.set_xlim(axes.get_xlim()[::-1]) #::-1 reverses the array
    if c.reverseY:
        axes.set_ylim(axes.get_ylim()[::-1])
    #Now re-draw the plot
    pl.draw()
    #(Insert commands needed to show plot, if necessary)
    return plotObj(None,fig) #Eventually we will use this to make subplots and do save option

# A basic contour plotter, which will plot a contour plot
# of a Numeric array. The x and y scales can optionally
# be specified using keyword arguments x and y. For example,
# if we want the x scale to be the array (or list) lat, and
# the y scale to be the array (or list) lon, we would call
# contour as contour(A,x=lat,y=lon).
def contour(A,**kwargs):
    #**ToDo: Add labeled contour lines, option to change contour levels,
    #        and to change palette.
    #
    fig = pl.figure() #Always start a new figure
    # Set axes if they have been specified
    # as keyword arguments
    if 'x' in kwargs.keys():
        x = kwargs['x']
    else:
        x = range(A.shape[1])
    if 'y' in kwargs.keys():
        y = kwargs['y']
    else:
        y = range(A.shape[0])
    cs = pl.contourf(x,y,A)
    cbar = pl.colorbar(cs)
    return plotObj(None,fig)
##    #The following allows an expert user to pass
##    #Ngl options directly to the plotter.
##    #ToDo: Note that
##    #options explicitly specified later will over-ride
##    #what the user specified. This should be fixed,
##    #by checking for things already specified. We should
##    #also allow for using this resource to specify a color
##    #map.
##    if 'resource' in kwargs.keys():
##        r = kwargs['resource']
##    else:
##        r = Dummy()
##
##    rw = Dummy()
##    #Set the color map
##    if 'colors' in kwargs.keys():
##        if (kwargs['colors'] == 'gray') or (kwargs['colors'] == 'grey') :
##            #Set the default greyscale
##            rw.wkColorMap = 'gsdtol'
##        else:
##            rw.wkColorMap = kwargs['colors']
##    else:
##        #Default rainbow color table
##        rw.wkColorMap = "temp1"
