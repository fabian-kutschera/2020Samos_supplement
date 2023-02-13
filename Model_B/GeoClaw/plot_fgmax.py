"""
Plot fgmax output from GeoClaw run.

"""

import matplotlib.pyplot as plt
import os
import numpy
from clawpack.geoclaw import fgmax_tools
from clawpack.visclaw import geoplot

def plot_fgmax_grid(outdir,plotdir):

    fg = fgmax_tools.FGmaxGrid()
    fg.outdir = outdir
    data_file = os.path.join(outdir, 'fgmax_grids.data')
    fg.read_fgmax_grids_data(fgno=1, data_file=data_file)
    fg.read_output()

    #clines_zeta = [0.01] + list(numpy.linspace(0.05,0.3,6)) + [0.5,1.0,10.0]
    #clines_zeta = [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 2.0]
    clines_zeta = [0.001, 0.01, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 1.0]
    colors = geoplot.discrete_cmap_1(clines_zeta)
    plt.figure(1)
    plt.clf()
    zeta = numpy.where(fg.B>0, fg.h, fg.h+fg.B)   # surface elevation in ocean
    print(numpy.min(zeta), numpy.max(zeta))
    plt.contourf(fg.X,fg.Y,zeta,clines_zeta,colors=colors)
    cbar = plt.colorbar()
    cbar.ax.set_title('Meters', pad=8)#, rotation=270)
    plt.contour(fg.X,fg.Y,fg.B,[0.],colors='k', linewidths=0.8)  # coastline

    # plot arrival time contours and label:
    #print(fg.arrival_time)
    #arrival_t = fg.arrival_time/3600.  # arrival time in hours
    #print(arrival_t[-1])
    #clines_t = numpy.linspace(0,4,9)  # hours
    #clines_t_label = clines_t[2::2]  # which ones to label 
    #print(clines_t_label)
    #clines_t_colors = ([.5,.5,.5],)
    #con_t = plt.contour(fg.X,fg.Y,arrival_t, clines_t,colors=clines_t_colors, linewidths=0.8) 
    #plt.clabel(con_t, clines_t_label)

    # fix axes:
    #plt.ticklabel_format(style='plain',useOffset=False)
    #plt.xticks(rotation=20)
    #plt.gca().set_aspect(1./numpy.cos(fg.Y.mean()*numpy.pi/180.))
    x_ticks = numpy.arange(23.0,29.0, step=1.0)
    sub_x = numpy.repeat("° E", len(x_ticks))
    plt.xticks(ticks=x_ticks, labels = numpy.char.add(x_ticks.astype(str), sub_x))
    y_ticks = numpy.arange(35.0, 41.0, step=1.0)
    sub_y = numpy.repeat("° N", len(y_ticks))
    plt.yticks(ticks=y_ticks, labels = numpy.char.add(numpy.round(y_ticks, 1).astype(str), sub_y))
    plt.title("Maximum positive tsunami height")
    plt.tight_layout()

    if not os.path.isdir(plotdir): 
        os.mkdir(plotdir)
    fname = os.path.join(plotdir, "amplitude_positive_Fra.png")
    plt.savefig(fname, dpi=300)
    print("Created ",fname)

if __name__=="__main__":
    plot_fgmax_grid(outdir='_output', plotdir='.')
