from __future__ import division, print_function, absolute_import
#----------Section 3: Plotting utilities-------------------------------
#These need to be localized, for systems that don't support Ngl
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
#-----------------------------------------------------------------------


#Note: On some older installations, Ngl has to
#be imported as PyNGL instead of Ngl.  If there
#is a trouble with the import, fiddle with that.
try:
    import Ngl
except:
    try:
        import PyNGL as Ngl
    except:
        print( "PyNGL was not found on your system.")
        print( "You will not be able to use graphics functions.")
        print( "PyNGL is available for Mac OSX,Linux and Windows/CygWin")
        print( "Alternately, modify the module ClimateGraphics.py")
        print( "to use some other graphics driver")
        raise('Graphics Import Error')

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
# to an eps file.
class plotObj:
    def __init__(self,workstation,plot,WorkstationResources=None):
        self.workstation = workstation
        self.plot = plot
        self.WorkstationResources = WorkstationResources
    #Deletes a plot window
    def delete(self):
        Ngl.destroy(self.workstation)
    def save(self,filename = 'plot'):
        #'eps' workstation doesn't work reliably, but 'ps' is OK
        weps = Ngl.open_wks('ps',filename,self.WorkstationResources)
        Ngl.change_workstation(self.plot,weps)
        Ngl.draw(self.plot)
        Ngl.change_workstation(self.plot,self.workstation)
        Ngl.destroy(weps)

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
lineColors =range(3,203,20)+range(3,203,20)
#Line colors for old versions of Ngl
#lineColors =range(3,16)+range(3,16)
lineThickness = [2,3,4,2,3,4,2,3,4,2,3,4,2,3,4,2,3,4,2,3,4]
plotSymbols = [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]
def plot(c):
    r = resource()
    r.trYReverse = c.reverseY
    r.trXReverse = c.reverseX
    r.trXLog = c.XlogAxis
    r.trYLog = c.YlogAxis
    #
    #Line styles
    r.xyLineColors = lineColors
    r.xyLineThicknesses = lineThickness
    #Plot title
    r.tiMainString = c.PlotTitle
    #Axis labels (ToDo: add defaults)
    #X and Y axis labels
    r.tiXAxisString = c.Xlabel
    r.tiYAxisString = c.Ylabel
    if c.switchXY:
        r.tiXAxisString,r.tiYAxisString = r.tiYAxisString,r.tiXAxisString

    #  Legends, for multicurve plot
    legends = []
    for id in c.listVariables():
        if not id == c.Xid:
            if len(c.label[id]) > 0:
                legends.append(c.label[id])
            else:
                legends.append(id)
    #ToDo: Add option to skip legends
    #Legends are redundant if there's only one curve
    if len(legends) > 1:
        r.pmLegendDisplayMode = "Always"
        r.xyExplicitLegendLabels = legends
    #
    #Suppress line drawing and just plot symbol for scatter plot curves
    r.xyMarkers = plotSymbols
    r.xyMarkerColors = lineColors
    r.xyMarkerSizeF = .01
    r.xyMarkLineModes = []
    for id in c.listVariables():
        if not id == c.Xid:
            if c.scatter[id]:
                r.xyMarkLineModes.append('Markers')
            else:
                r.xyMarkLineModes.append('Lines')
    #
    w = Ngl.open_wks('x11','Climate Workbook')
    if c.switchXY:
        plot = Ngl.xy(w,c.Y(),c.X(),r)
    else:
        plot = Ngl.xy(w,c.X(),c.Y(),r)
    #
    #Now draw the plot
    r.nglDraw = True
    Ngl.panel(w,[plot],[1,1],r)
    return plotObj(w,plot) #So that user can delete window or save plot

# A basic contour plotter, which will plot a contour plot
# of a Numeric array. The x and y scales can optionally
# be specified using keyword arguments x and y. For example,
# if we want the x scale to be the array (or list) lat, and
# the y scale to be the array (or list) lon, we would call
# contour as contour(A,x=lat,y=lon).
def contour(A,**kwargs):
    #The following allows an expert user to pass
    #Ngl options directly to the plotter.
    #ToDo: Note that
    #options explicitly specified later will over-ride
    #what the user specified. This should be fixed,
    #by checking for things already specified. We should
    #also allow for using this resource to specify a color
    #map.
    if 'resource' in kwargs.keys():
        r = kwargs['resource']
    else:
        r = Dummy()
    #
    r.cnFillOn = True #Uses color fill
    # Set axes if they have been specified
    # as keyword arguments
    if 'x' in kwargs.keys():
        r.sfXArray = kwargs['x']
    if 'y' in kwargs.keys():
        r.sfYArray = kwargs['y']
    #
    # Now create the plot

    rw = Dummy()
    #Set the color map
    if 'colors' in kwargs.keys():
        if (kwargs['colors'] == 'gray') or (kwargs['colors'] == 'grey') :
            #Set the default greyscale
            rw.wkColorMap = 'gsdtol'
        else:
            rw.wkColorMap = kwargs['colors']
    else:
        #Default rainbow color table
        rw.wkColorMap = "temp1"

    w = Ngl.open_wks('x11','Climate Workbook',rw)
    r.nglDraw = False
    r.nglFrame = False
    plot = Ngl.contour(w,A,r)
    #Now draw the plot
    r.nglDraw = True
    Ngl.panel(w,[plot],[1,1],r)
    return plotObj(w,plot,rw) #So user can delete or save plot
