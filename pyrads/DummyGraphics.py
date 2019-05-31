from __future__ import division, print_function, absolute_import
#----------Section 3: Plotting utilities-------------------------------
#This is a dummy graphics routine, to import if a graphics driver
#is not found. It is the fallback import if the import of ClimateGraphics
#fails in ClimateUtilities.  This dummy routine allows courseware
#scripts to be run without the user needing to comment out plot commands.
#It also can be used as a template for customizing the plotter interface
#to work with you locally favored Python graphics driver (e.g. MatPlotLib).
#----------------------------------------------------------------------=

#A dummy class useful for passing parameters.
#Just make an instance, then add new members
#at will.  Not sure why this also has to be defined
#here, since it's defined in ClimateUtilities, but
#it seems to be necessary
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
# to an eps file.
class plotObj:
    def __init__(self,workstation,plot,WorkstationResources = None):
        self.workstation = workstation
        self.plot = plot
        self.WorkstationResources = WorkstationResources
    #Deletes a plot window
    def delete(self):
        #Deletes a plot window, cleans up
        pass
    def save(self,filename = 'plot'):
        #Saves a plot to a file
        pass

#Makes a line pot from a Curve object
def plot(c):
    print( "Plotting not implemented")
    #Set axis options according to information
    #in the curve object c.
    #
    #c.reverseY:
    #   if True, reverse the Y axis
    #c.reverseX:
    #   if True, reverse the X axis
    #c.XlogAxis:
    #   if True, use a logarithmic X axis
    #c.YlogAxis:
    #   if True, use a logarithmic Y axis
    #
    #Customize Line styles and colors here
    #
    #Set thePlot title
    #c.PlotTitle:
    #   String containing the plot title
    #Axis labels
    #X and Y axis labels
    #c.Xlabel:
    #   String containing the X axis label
    #c.Ylabel:
    #   String containing the Y axis label
    #
    #Interchange the X and Y axes
    if c.switchXY:
        pass
        #If True, exchang the axes

    #  Legends, for multicurve plot
    legends = []
    for id in c.listVariables():
        if not id == c.Xid:
            if len(c.label[id]) > 0:
                legends.append(c.label[id])
            else:
                legends.append(id)

    #
    #Suppress line drawing and just plot symbol for scatter plot curves
    #Customize plotting symbols and marker sizes if desired
    for id in c.listVariables():
        if not id == c.Xid:
            if c.scatter[id]:
                #If True, treat this data column as a scatter
                #plot and don't draw lines
                pass
            else:
                #If False, draw lines (default)
                pass
    #
    #Initialize the plot window here, if necessary.
    #w is the handle to the plot window
    w = None
    if c.switchXY:
        #Put in the command for doing the line plot here
        #Do the plot with the X and Y axes in the usual
        #order
        pass
    else:
        #Put in the command for doing the line plot here,
        #but switch the order of the axes.
        pass
    #
    #Now draw the plot
    #Depending on your graphics software, the preceding
    #command may already have drawn the plot. In some graphics
    #packages, after the plot is created, another command needs
    #to be executed in orter to display it. Either way, the
    #variable plotHandle below is a handle referring to the plot.
    #It is used to build a plotObj that can be used for further
    #manipulation of the plot such as saving or deleting. In some
    #graphics packages, which give control over such things from
    #menus in the plot window, the use of a plotObj for this may
    #be unnecessary.
    plotHandle = None
    return plotObj(w,plotHandle) #So that user can delete window or save plot

# A basic contour plotter, which will plot a contour plot
# of a numpy array. The x and y scales can optionally
# be specified using keyword arguments x and y. For example,
# if we want the x scale to be the array (or list) lat, and
# the y scale to be the array (or list) lon, we would call
# contour as contour(A,x=lat,y=lon).
def contour(A,**kwargs):
    print( "Plotting not implemented")
    #The following allows an expert user to pass
    # options directly to the plotter.
    if 'resource' in kwargs.keys():
        r = kwargs['resource']
    else:
        r = Dummy()
    #
    r.cnFillOn = True #Use color fill

    if 'x' in kwargs.keys():
        #Set the X array for the contour plot
        XArray = kwargs['x']
    if 'y' in kwargs.keys():
        #Set the Y array for the contour plot
        YArray = kwargs['y']
    #
    # Now create the plot
    rw = Dummy()
    #Set the color map
    if 'colors' in kwargs.keys():
        if (kwargs['colors'] == 'gray') or (kwargs['colors'] == 'grey') :
            #Set the default greyscale
            #(Substitute the appropriate command for your driver)
            rw.wkColorMap = 'gsdtol'
        else:
            rw.wkColorMap = kwargs['colors']
    else:
        #Default rainbow color table
        rw.wkColorMap = "temp1"

    #Open/initialize a plot window
    w = None
    #Make the plot. plotHandle is the handle returned
    #by the plotter, used for further manipulation of the
    #plot. (Redundant for some kinds of plotting packages)
    plotHandle = None
    #Now draw the plot, if your driver needs this as a separate step
    #(Insert command for drawing the plot here, e.g.
    #ShowPlot(plotHandle).
    #
    #Return a plotObj with the necessary data for further
    #manipulation.
    return plotObj(w,plotHandle,rw) #So user can delete or save plot
