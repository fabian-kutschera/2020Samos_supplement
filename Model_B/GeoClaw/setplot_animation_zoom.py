
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt


#--------------------------
def setplot(plotdata=None):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot
    from clawpack.visclaw.data import ClawPlotData

    from numpy import linspace

    if plotdata is None:
        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items data


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    color_min = -0.4
    color_max = 0.4
    x_min = 23
    x_max = 28
    y_min = 35
    y_max = 40
    def fixup(current_data):
        import pylab
        t = current_data.t
        hours = int(t//3600)
        minutes = (t-hours*3600)/60
        if hours==0:
            time_string = '    '
        else:
            time_string = '%d h ' %(hours)
        time_string = time_string + '%d min' %(minutes)

        pylab.title('')#'Surface at %4.2f hours' % t, fontsize=15)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
        pylab.xticks(ticks = np.arange(x_min, x_max+1, 1), labels=[" ", "24° E", "25° E", "26° E", "27° E", "28° E"])
        pylab.yticks(ticks = np.arange(y_min, y_max+1, 1), labels=["35° N", "36° N", "37° N", "38° N", "39° N", "40° N"])
        pylab.text(26.9, 39.5, time_string, fontsize=15)
        pylab.grid(True)

        #add colorbar
        fig = pylab.gcf()
        ax = pylab.gca()
        a = np.array([[color_min, color_max]])

        from cmcrameri import cm
#        img = pylab.imshow(a, cmap=geoplot.tsunami_colormap)
        img = pylab.imshow(a, cmap=cm.cork)

        #cax = fig.add_axes([0.17, 0.3, 0.01, 0.20])
        cax = fig.add_axes([0.16, 0.18, 0.015, 0.18])
        cb = pylab.colorbar(cax=cax, ticks=np.linspace(color_min, color_max, 5))
        cb.ax.tick_params(which='major', labelsize=15, length=8, width =1, direction='out')
        cb.ax.tick_params(which = 'minor', length = 6, width =1.0)
        pylab.sca(ax)
        pylab.text(23.1, 36.7, 'ssha (m)', va="top", family="sans-serif",  clip_on=True, fontsize=15, zorder=6, rotation='horizontal')#vertical')
        #pylab.text(22.1, 36.35, 'ssha (m)', va="top", family="sans-serif",  clip_on=True, fontsize=15, zorder=6, rotation='horizontal')#vertical')
        ax.tick_params(direction='in')
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
       
        #plot Gauge
        plt.plot(25.152692, 35.348475, marker="^", color="cyan", markersize=12, markeredgecolor="black") #HRAK
        plt.plot(27.423453, 37.032170, marker="^", color="cyan", markersize=12, markeredgecolor="black") #BODR
        plt.plot(27.287792, 36.898362, marker="^", color="cyan", markersize=12, markeredgecolor="black") #KOS1
        plt.plot(27.303632, 36.891013, marker="^", color="cyan", markersize=12, markeredgecolor="black") #KOS2
        plt.plot(24.945808, 37.439969, marker="^", color="cyan", markersize=12, markeredgecolor="black") #SYRO
        plt.plot(26.370550, 38.971880, marker="^", color="cyan", markersize=12, markeredgecolor="black") #PLOM
        plt.plot(26.921840, 35.418600, marker="^", color="cyan", markersize=12, markeredgecolor="black") #NOA-03
        #plot Gauge names
        plt.text(25.152692-0.15, 35.348475-0.17, s="hrak", fontsize=14, bbox=dict(facecolor='white', alpha=0.5, edgecolor='white', pad=0.5))
        plt.text(26.921840-0.24, 35.418600-0.17, s="noa-03", fontsize=14, bbox=dict(facecolor='white', alpha=0.5, edgecolor='white', pad=0.5))
        plt.text(26.370550-0.15, 38.971880+0.1, s="plom", fontsize=14, bbox=dict(facecolor='white', alpha=0.5, edgecolor='white', pad=0.5))
        plt.text(24.945808-0.15, 37.439969-0.17, s="syro", fontsize=14, bbox=dict(facecolor='white', alpha=0.5, edgecolor='white', pad=0.5))
        plt.text(27.423453-0.15, 37.032170+0.1, s="bodr", fontsize=14, bbox=dict(facecolor='white', alpha=0.5, edgecolor='white', pad=0.5))
        plt.text(27.287792-0.15, 36.898362-0.17, s="kos1", fontsize=14, bbox=dict(facecolor='white', alpha=0.5, edgecolor='white', pad=0.5))
        plt.text(27.287792-0.15, 36.891013-0.3, s="kos2", fontsize=14, bbox=dict(facecolor='white', alpha=0.5, edgecolor='white', pad=0.5))


        #plot Jason location
        #import pyproj
        #lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        #myproj=pyproj.Proj('EPSG:5936')
#        Jason = np.genfromtxt('./dart_sta.csv', delimiter=",", skip_header=1)
#        xyz = pyproj.transform(myproj, lla, Jason[:,0],Jason[:,1], radians=False)
#        pylab.plot(Jason[:,0], Jason[:,1], color='k', linewidth=2)
       
    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)
    plotfigure.kwargs = {'figsize':[7.5*1.5*30./35.,7.5*1.5], 'facecolor': 'white'}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.scaled = True
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    from cmcrameri import cm    
#    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmap= cm.cork
    plotitem.pcolor_cmin = color_min
    plotitem.pcolor_cmax = color_max
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0,0]
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land

    from matplotlib import cm as cm1
    from clawpack.visclaw import colormaps
    navajowhite = colormaps.make_colormap({0.:'navajowhite', 1.:'navajowhite'})
    plotitem.pcolor_cmap = navajowhite
    #plotitem.pcolor_cmap = geoplot.land_colors # for ETOPO Land

    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0,0]
    plotitem.patchedges_show = 0
    plotaxes.xlimits = [x_min,x_max]
    plotaxes.ylimits = [y_min,y_max]
    #plotaxes.xlimits = [26,28.5]
    #plotaxes.ylimits = [37,38.5]
   # plotaxes.afteraxes = addgauges

    # Add dashed contour line for shoreline
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [0.]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.amr_contour_show = [1]  # show contours only on finest level
    plotitem.kwargs = {'linewidths': 0.5}


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface at gauges', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    #plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    def add_zeroline(current_data):
        from pylab import plot, legend, xticks, floor, axis, xlabel, ylabel
        t = current_data.t 
        gaugeno = current_data.gaugeno
        if gaugeno==1 or gaugeno==14:
            plt.title("Surface at Hrakleio")
        if gaugeno==2 or gaugeno==8:
            plt.title("Surface at Bodrum")
        if gaugeno==3 or gaugeno==9:
            plt.title("Surface at Kos")
        if gaugeno==4 or gaugeno==10:
            plt.title("Surface at Kos Marina")
        if gaugeno==5 or gaugeno==11:
            plt.title("Surface at Syros")
        if gaugeno==6 or gaugeno==12:
            plt.title("Surface at Plomari")
        if gaugeno==7 or gaugeno==13:
            plt.title("Surface at Kasos")
        if gaugeno==15:
            plt.title("Surface at Syn1")
        if gaugeno==16:
            plt.title("Surface at Syn2")
        if gaugeno==17:
            plt.title("Surface at Syn3")
        if gaugeno==18:
            plt.title("Surface at Syn4")
        if gaugeno==19:
            plt.title("Surface at Syn5")
        if gaugeno==20:
            plt.title("Surface at Syn6")
        if gaugeno==21:
            plt.title("Surface at Syn7")
        plot(t, 0*t, 'k')
        n = int(floor(t.max()/3600.) + 2)
        xticks([3600*i for i in range(n)], ['%i' % i for i in range(n)])
        xlabel('Time [h]')
        ylabel('Ssha [m]')
        spacing = 0.15
        plt.subplots_adjust(left=spacing)
                   
        from matplotlib.ticker import FormatStrFormatter
        ax = plt.gca()
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    plotaxes.afteraxes = add_zeroline


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.parallel = True
    plotdata.print_framenos = [15, 30, 45, 60]          # list of frames to print
#    plotdata.print_framenos = 'all'
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = False                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

# To plot gauge locations on pcolor or contour plot, use this as
# an afteraxis function:

#def addgauges(current_data):
#    from clawpack.visclaw import gaugetools
#    gaugetools.plot_gauge_locations(current_data.plotdata, \
#         gaugenos='all', format_string='ko', add_labels=True)
