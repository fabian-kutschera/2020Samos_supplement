# -*- coding: utf-8 -*-

from pylab import *
import os
from clawpack.visclaw import plottools, geoplot
from clawpack.visclaw import animation_tools
from matplotlib import animation, colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from clawpack.geoclaw import fgout_tools
from datetime import timedelta 
from cmcrameri import cm
import matplotlib

rcParams["font.family"] = "sans-serif"
ps = 18
rcParams.update({"font.size": ps})
rcParams["font.family"] = "sans"
matplotlib.rc("xtick", labelsize=ps)
matplotlib.rc("ytick", labelsize=ps)

x1, x2, y1, y2 = 25.4, 27.5, 37.6, 38.3 #Samos

fgno = 1
outdir = '_output'
output_format = 'binary'

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, output_format)

# Plot one frame of fgout data 
fgframe = 4 #7 #4 # frame -1 = minutes
fgout = fgout_grid.read_frame(fgframe)

fig = plt.figure()
fig.set_size_inches(10, 5)
ax1 = plt.gca()

eta_water = where(fgout.h > 0, fgout.eta, nan)
img = imshow(flipud(eta_water.T), extent=fgout.extent_edges, cmap=cm.cork, aspect="auto")
clim(-0.8, 0.8)
plt.contour(fgout.X,fgout.Y,fgout.B,[0.],colors='k', linewidths=0.8)  # coastline

# ticks
plt.xlim(x1,x2)
x_ticks = np.arange(25.5,28.0, step=0.5)
sub_x = np.repeat("° E", len(x_ticks))
plt.xticks(ticks=x_ticks, labels = np.char.add(x_ticks.astype(str), sub_x))
plt.ylim(y1,y2)
y_ticks = np.arange(37.6, 38.2, step=0.2)
sub_y = np.repeat("° N", len(y_ticks))
plt.yticks(ticks=y_ticks, labels = np.char.add(np.round(y_ticks, 1).astype(str), sub_y))

#title('Surface at time %s' % timedelta(seconds=fgout.t))
plt.title("t = {} min".format(int(fgout.t/60)))

# colorbar
plt.text(27.0, 38.15, "ssha (m)", horizontalalignment='center', verticalalignment='center',
             fontsize=ps, rotation="horizontal")
cbaxes = inset_axes(ax1, width="45.2%", height="5%", loc="upper right")
cbar = plt.colorbar(img, orientation="horizontal", cax=cbaxes, ticks=[-0.8,0.0,0.8], extend="both")
cbar.ax.tick_params(which="major", labelsize=ps, length=10, width=2, direction="inout")
cbar.ax.tick_params(which="minor", length=6, width=1.0)
cbar.ax.minorticks_on()

fname = 'fgout_frame%s_bay_Ryo.png' % str(fgframe).zfill(4)
savefig(fname, dpi=300)#, bbox_inches='tight')
print('Created ',fname)

