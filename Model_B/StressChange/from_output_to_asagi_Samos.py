#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:42:41 2022

@author: marcmath
"""

# #!/usr/bin/env python3
import seissolxdmf
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.colors as colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import scipy.linalg
from scipy.interpolate import griddata
from netCDF4 import Dataset
import scipy.ndimage as ndima
import h5py
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import signal

# update default plotting parameters
plt.rcParams["font.family"] = "serif"
plt.rc("xtick", labelsize=14) 
plt.rc("ytick", labelsize=14)

################
## Functions ##
################

def ReadHdf5PosixForBoundaryPlotting(filename,BC):
   sx = seissolxdmf.seissolxdmf(filename)
   xyz = sx.ReadGeometry()
   tetra = sx.ReadConnect()
   nElements=np.shape(tetra)[0]
   boundary = sx.ReadDataChunk('boundary', firstElement=0, nchunk=nElements, idt=0)
   NS=0
   SufaceId= int(BC)
   for faceId in range(0,4):
      boundaryFace = (boundary >> (faceId*8)) & 0xFF;
      if SufaceId==0:
         tetraId = np.where(boundaryFace>0)[0]
      else:
         tetraId = np.where(boundaryFace==SufaceId)[0]
      NS = NS + len(tetraId)
   assert(NS!=0)
   connect=np.zeros((NS,3), dtype=int)
   BC=np.zeros((1,NS))

   s_vert=np.zeros((4,3), int)
   s_vert[0,:] = [0,2,1];   s_vert[1,:] = [0,1,3];   s_vert[3,:] = [0,3,2];   s_vert[2,:] = [1,2,3];

   currentindex = 0
   for faceId in range(0,4):
      boundaryFace = (boundary >> (faceId*8)) & 0xFF;
      if SufaceId==0:
         tetraId = np.where(boundaryFace>0)[0]
      else:
         tetraId = np.where(boundaryFace==SufaceId)[0]
      for idBound in range(0, len(tetraId)):
          trinodes = tetra[ tetraId[idBound], s_vert[faceId,:]]
          connect[currentindex,:] = trinodes
          BC[0,currentindex] = boundaryFace[tetraId[idBound]]
          currentindex = currentindex +1

   return xyz,connect

def plot_field(x, y, z, s=30, xlab=None, ylab=None, title=None,  vmin=None, vmax=None):
    plt.scatter(x, y, s, z, marker='s', vmin=vmin, vmax=vmax, cmap='Spectral')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.colorbar()  

def project_3D_pts_on_plane(p_xyz, in_xyz):
    # perpendicular projection of 3D points onto a plane
    # xyzp = x,y,z coordinates of at least 3 pts belonging to the plane
    # in_xyz = x,y,z coordinates of the points to project
    # out_xyz = x,y,z coordinates of the projected points
    AB = p_xyz[1, :] - p_xyz[0, :]
    AC = p_xyz[2, :] - p_xyz[0, :]
    
    n = np.cross(AB, AC)
    d = np.dot(p_xyz[0, :], n.T)

    t0 = - (n[0]*in_xyz[:, 0] + n[1]*in_xyz[:, 1] + n[2]*in_xyz[:, 2] - d)/(n[0]**2 + n[1]**2 + n[2]**2)
    X = in_xyz[:, 0] + n[0] * t0
    Y = in_xyz[:, 1] + n[1] * t0
    Z = in_xyz[:, 2] + n[2] * t0
    
    out_xyz = np.array([X, Y, Z]).T
    
    return out_xyz
        
def plane_coord(x1, x2, y1, y2, dip, strike, depth1, depth2):
    # depth 1: depth of the top bottom of the fault (= or > 0)
    # depth 2: en profondeur (< 0)
    # (x1,y1), (x2,y2) coordinates of the two top corners of the plane
    # out : xp,yp,zp  coordinates of the four corners of the plane
    
    dx = np.sin(np.deg2rad(strike-90)) * ( abs(depth2-depth1) / np.tan(np.deg2rad(dip)))
    dy = np.cos(np.deg2rad(strike-90)) * ( abs(depth2-depth1) / np.tan(np.deg2rad(dip)))

    xp = [x1, x2, x2 - dx, x1 - dx, x1]
    yp = [y1, y2, y2 + dy, y1 + dy, y1]
    zp = [depth1, depth1, depth2, depth2, depth1]
    
    return xp, yp, zp


def compute_unit_vec_trans_plane(xp, yp, zp):
    # xp,yp,zp coordinate of the 4 corners of the plane (but only 3 pts really needed)
    # hw, hh, th, tw: matrices and scalar needed to go from the flobal coordinates system
    # to the local coordinates system defined by the plane. 
    hw = np.array([xp[1], yp[1], zp[1]] - np.array([xp[2], yp[2], zp[2]])) # fault 3
    hw = hw / np.linalg.norm(hw)
    hh = np.array([xp[3], yp[3], zp[3]]) - np.array([xp[2], yp[2], zp[2]])
    hh = hh / np.linalg.norm(hh)
    th = -np.dot(np.array([xp[2], yp[2], zp[2]]), hh)
    tw = -np.dot(np.array([xp[2], yp[2], zp[2]]), hw)
    
    return hw, hh, th, tw 

def interpolate_2dgrid(xmin,xmax,ymin,ymax,step,data,method='cubic'):    
    # xmin,xmax,ymin,ymax,step : info needed to discretize the plane
    # data (3 x n matrix: x,y coordinates and value of the n data we want to interpolate)   
    l = xmax - xmin
    w = ymax - ymin

    nl = int(l / step)
    nw = int(w / step)

    xr = np.linspace(xmin, xmax, nl)
    yr = np.linspace(ymin, ymax, nw)

    data[:, 0] = data[:, 0]
    data[:, 1] = data[:, 1]

    grid_x, grid_y = np.meshgrid(xr, yr)

    data_interp = griddata(data[:, 0:2], data[:, 2], (grid_x, grid_y), method='cubic')
    data_interp = np.nan_to_num(data_interp)

    return grid_x, grid_y, data_interp

def writeNetcdf4Paraview(sname, x, y, aName, aData):
    "create a netcdf file readable by paraview (but not by ASAGI)"
    fname = sname + "_paraview.nc"
    print("writing " + fname)
    ####Creating the netcdf file
    nx = x.shape[0]
    ny = y.shape[0]
    rootgrp = Dataset(fname, "w", format="NETCDF4")
    rootgrp.createDimension("u", nx)
    rootgrp.createDimension("v", ny)

    vx = rootgrp.createVariable("u", "f4", ("u",))
    vx[:] = x
    vy = rootgrp.createVariable("v", "f4", ("v",))
    vy[:] = y
    for i in range(len(aName)):
        vTd = rootgrp.createVariable(aName[i], "f4", ("v", "u"))
        vTd[:, :] = aData[i][:, :]
    rootgrp.close()


def writeNetcdf4SeisSol(sname, x, y, aName, aData):
    "create a netcdf file readable by ASAGI (but not by paraview)"
    ########## creating the file for SeisSol
    fname = sname + "_ASAGI.nc"
    print("writing " + fname)
    ####Creating the netcdf file
    nx = x.shape[0]
    ny = y.shape[0]

    rootgrp = Dataset(fname, "w", format="NETCDF4")

    rootgrp.createDimension("u", nx)
    rootgrp.createDimension("v", ny)

    vx = rootgrp.createVariable("u", "f4", ("u",))
    vx[:] = x
    vy = rootgrp.createVariable("v", "f4", ("v",))
    vy[:] = y
    ldata4 = [(name, "f4") for name in aName]
    ldata8 = [(name, "f8") for name in aName]
    mattype4 = np.dtype(ldata4)
    mattype8 = np.dtype(ldata8)
    mat_t = rootgrp.createCompoundType(mattype4, "material")

    # this transform the 4 D array into an array of tuples
    arr = np.stack([aData[i] for i in range(len(aName))], axis=2)
    newarr = arr.view(dtype=mattype8)
    newarr = newarr.reshape(newarr.shape[:-1])
    mat = rootgrp.createVariable("data", mat_t, ("v", "u"))
    mat[:] = newarr
    rootgrp.close()

def write_yaml_file(yamlfile, asagifile, segmentID, hh, hw, th, tw):

    nb_seg = len(asagifile)
    template_yaml_begining = f"""!Switch
!Switch
[Ts0,Td0,Pn0]: !EvalModel
    parameters: [tso, tdo, pno]
    model: !Switch
        [tso, tdo, pno]: !Any
          components:"""
    template_yaml_end = f"""
    components: !FunctionMap
       map:
          Ts0: return -tso;
          Td0: return tdo;
          Pn0: return pno;
        """
    fname = f"{yamlfile}.yaml"
    with open(fname, "w") as fid:
        fid.write(template_yaml_begining)
        for i in range(0, nb_seg):
            template_yaml_segment = f"""
            - !GroupFilter
               groups: {segmentID[i]}
               components: !AffineMap
                    matrix:
                      ua: [{hh[i][0]}, {hh[i][1]}, {hh[i][2]}]
                      ub: [{hw[i][0]}, {hw[i][1]}, {hw[i][2]}]
                    translation:
                      ua: {th[i]}
                      ub: {tw[i]}
                    components: !Any
                      - !ASAGI
                          file: {asagifile[i]}_ASAGI.nc
                          parameters: [tso, tdo, pno]
                          var: data
                          interpolation: linear
                      - !ConstantMap
                        map:
                          tso: 0.0
                          tdo: 0.0
                          pno: 0.0 """
            fid.write(template_yaml_segment)
        fid.write(template_yaml_end)
    print(f"done writing {fname}")
###########################
# ## Open SeisSol output ###
# ##########################

model = "Fra" # "Fra"
version = "v4" # "v3"

file='/import/freenas-m-05-seissol/kutschera/MAthesis_fkutschera/simulations/Samos_final_{}_{}/{}_{}_WL-fault.xdmf'.format(model,version,model, version)

#file='/home/marcmath/Bureau/Work/tibet/rupture_sim/slip2stress/supermuc/output/SFS-fault.xdmf'
sx = seissolxdmf.seissolxdmf(file)
nElements = sx.ReadNElements()
connect = sx.ReadConnect()
ndt = sx.ReadNdt()
xyz = sx.ReadGeometry()
tdo = sx.ReadData('Td0',ndt-1)
pno = sx.ReadData('Pn0',ndt-1)
tso = sx.ReadData('Ts0',ndt-1)
bary = (xyz[connect[:,0]] + xyz[connect[:,1]] + xyz[connect[:,2]])/3.
npt = np.shape(bary)[0]


#tso[tso > 20e6] = 20e6;
#tso[tso < -20e6] = -20e6;
#
#tdo[tdo > 20e6] = 20e6; #put the outlier to something close to the average value
#tdo[tdo < -20e6] = -20e6;
#
#pno[pno > 20e6] = 20e6;
#pno[pno < -20e6] = -20e6;

tso[tso > 10e6] = 10e6;
tso[tso < -10e6] = -10e6;

tdo[tdo > 10e6] = 10e6; #put the outlier to something close to the average value
tdo[tdo < -10e6] = -10e6;

pno[pno > 10e6] = 10e6;
pno[pno < -10e6] = -10e6;


#
#filename = "/import/freenas-m-05-seissol/kutschera/MAthesis_fkutschera/mesh/mesh_noWL_Samos_{}.xdmf".format(model)
filename = "/import/freenas-m-05-seissol/kutschera/MAthesis_fkutschera/mesh/mesh_waterlayer_Samos_{}.xdmf".format(model)

#
#with h5py.File(filename, "r") as f:
#    # Print all root level object names (aka keys) 
#    # these can be group or dataset names 
#    print("Keys: %s" % f.keys())
#    # get first object name/key; may or may NOT be a group
#    a_group_key = list(f.keys())[3]
#
#    # get the object type for a_group_key: usually group or dataset
#    print(type(f[a_group_key])) 
#
#    ds_arr = f[a_group_key][()]  # returns as a numpy array
#    print(np.shape(ds_arr))


##########################################
# # ISOLATE THE DIFFERENT FAULT SEGMENTS ##
# #########################################
# Maduo fault is composed of three segments, so we need to isolate the coordinates of each segment in the mesh, in order to produce one ASAGI file per segment

xyz3, connect3 = ReadHdf5PosixForBoundaryPlotting(filename,3)
npt3 = np.shape(xyz3)[0]
bary3 = (xyz3[connect3[:,0]] + xyz3[connect3[:,1]] + xyz3[connect3[:,2]])/3.

#xyz65, connect65 = ReadHdf5PosixForBoundaryPlotting(filename,65)
#npt65 = np.shape(xyz65)[0]
#bary65 = (xyz65[connect65[:,0]] + xyz65[connect65[:,1]] + xyz65[connect65[:,2]])/3.
#
#xyz66, connect66 = ReadHdf5PosixForBoundaryPlotting(filename,66)
#npt66 = np.shape(xyz66)[0]
#bary66 = (xyz66[connect66[:,0]] + xyz66[connect66[:,1]] + xyz66[connect66[:,2]])/3.


# Retrieve the ID corresponding to each segment in the SeisSol output 
ind3 = np.argwhere(np.in1d(bary[:,0],bary3[:,0]))
#ind65 = np.argwhere(np.in1d(bary[:,0],bary65[:,0]))
#ind66 = np.argwhere(np.in1d(bary[:,0],bary66[:,0]))


bary3=bary[ind3[:,0],:]
#bary65=bary[ind65[:,0],:]
#bary66=bary[ind66[:,0],:]

tso3=tso[ind3]
tdo3=tdo[ind3]
pno3=pno[ind3]

#tso65=tso[ind65]
#tdo65=tdo[ind65]
#pno65=pno[ind65]

#tso66=tso[ind66]
#tdo66=tdo[ind66]
#pno66=pno[ind66]


#################
### SEGMENT 1 ###
#################
if model == "Fra":
    if version == "v3":
        ua3 = [-0.9958740467789645, -0.0907462558571351, 0.0]
        ub3 = [0.06963219135524702, -0.762695724677607, 0.6429981255692478]
        ta3 = 79913.50124701673
        tb3 = 26514.39182841975
    elif version == "v4": # actually they are the same
        ua3 = [-0.9958740467789645, -0.0907462558571351, 0.0]
        ub3 = [0.06963219135524702, -0.762695724677607, 0.6429981255692478]
        ta3 = 79913.50124701673
        tb3 = 26514.39182841975
elif model =="Ryo":
    ua3 = [-0.9999802812354174, -0.006279899707447305, 0.0]
    ub3 = [0.005009900725918649, -0.7985471624608284, 0.6019113973173094]
    ta3 = 103931.54917899476
    tb3 = 40612.38129771323    
else:
    print("Sure")

xa  = np.dot(bary3, ua3) + ta3
xb  = np.dot(bary3, ub3) + tb3
xyb3 = np.transpose(np.vstack((xa, xb))) 

step1 = 100
grid_x1, grid_y1, tsoi3 = interpolate_2dgrid(np.min(xyb3[:,0]),np.max(xyb3[:,0]),np.min(xyb3[:,1]),np.max(xyb3[:,1]),step1,np.c_[xyb3,tso3],'cubic')
grid_x1, grid_y1, tdoi3 = interpolate_2dgrid(np.min(xyb3[:,0]),np.max(xyb3[:,0]),np.min(xyb3[:,1]),np.max(xyb3[:,1]),step1,np.c_[xyb3,tdo3],'cubic')
grid_x1, grid_y1, pnoi3 = interpolate_2dgrid(np.min(xyb3[:,0]),np.max(xyb3[:,0]),np.min(xyb3[:,1]),np.max(xyb3[:,1]),step1,np.c_[xyb3,pno3],'cubic')

tsoi3 = ndima.median_filter(tsoi3,size=10)
tdoi3 = ndima.median_filter(tdoi3,size=30)
pnoi3 = ndima.median_filter(pnoi3,size=10)

#wind = signal.windows.hann(np.shape(tsoi3)[1])
wind = signal.windows.tukey(np.shape(tsoi3)[1],alpha=0.5)#0.5) # left and right
tsoi3tuk = np.zeros(np.shape(tsoi3))
tdoi3tuk = np.zeros(np.shape(tdoi3))
pnoi3tuk = np.zeros(np.shape(pnoi3))

for i in range(0, np.shape(tsoi3)[0]):
    tsoi3tuk[i,:] = tsoi3[i,:] * wind
    tdoi3tuk[i,:] = tdoi3[i,:] * wind
    pnoi3tuk[i,:] = pnoi3[i,:] * wind
    
# testing
wind2 = signal.windows.tukey(np.shape(tsoi3)[0],alpha=0.9) #0.3) # top and bottom

#wind2 = signal.windows.hann(np.shape(tsoi3)[0])
for i in range(0, np.shape(tsoi3)[1]):
    tsoi3tuk[:,i] = tsoi3[:,i] * wind2
    tdoi3tuk[:,i] = tdoi3[:,i] * wind2
    pnoi3tuk[:,i] = pnoi3[:,i] * wind2

file_prefix = "{}_{}_tso_tdo_pno".format(model,version)
ldataName = ["Ts0", "Td0","Pn0"]
lgridded_myData = [tsoi3tuk, tdoi3tuk,pnoi3tuk]
writeNetcdf4SeisSol(file_prefix, grid_x1[0,:], grid_y1[:,1], ldataName, lgridded_myData)
writeNetcdf4Paraview(file_prefix, grid_x1[0,:], grid_y1[:,1], ldataName, lgridded_myData)  
 
#################
### SEGMENT 2 ###
#################
#ua65 = [-0.9869469544191454, 0.16104567415107313, 0.0]
#ub65 = [-0.02352603981686094, -0.14417620013200017, 0.9892723329629879]
#ta65 = -1206250.1140838875
#tb65 = 287185.1256218455
#
#xa  = np.dot(bary65, ua65) + ta65
#xb  = np.dot(bary65, ub65) + tb65
#xyb65 = np.transpose(np.vstack((xa, xb))) 
#
#step1 = 100
#grid_x2, grid_y2, tsoi65 = interpolate_2dgrid(np.min(xyb65[:,0]),np.max(xyb65[:,0]),np.min(xyb65[:,1]),np.max(xyb65[:,1]),step1,np.c_[xyb65,tso65],'cubic')
#grid_x2, grid_y2, tdoi65 = interpolate_2dgrid(np.min(xyb65[:,0]),np.max(xyb65[:,0]),np.min(xyb65[:,1]),np.max(xyb65[:,1]),step1,np.c_[xyb65,tdo65],'cubic')
#grid_x2, grid_y2, pnoi65 = interpolate_2dgrid(np.min(xyb65[:,0]),np.max(xyb65[:,0]),np.min(xyb65[:,1]),np.max(xyb65[:,1]),step1,np.c_[xyb65,pno65],'cubic')
#
#tsoi65 = ndima.median_filter(tsoi65,size=10)
#tdoi65 = ndima.median_filter(tdoi65,size=30)
#pnoi65 = ndima.median_filter(pnoi65,size=10)
#
#
#wind = signal.windows.tukey(np.shape(tsoi65)[1],alpha=0.1)
#
#tsoi65tuk = np.zeros(np.shape(tsoi65))
#tdoi65tuk = np.zeros(np.shape(tdoi65))
#pnoi65tuk = np.zeros(np.shape(pnoi65))
#
#for i in range(0, np.shape(tsoi65)[0]):
#    tsoi65tuk[i,:] = tsoi65[i,:] * wind
#    tdoi65tuk[i,:] = tdoi65[i,:] * wind 
#    pnoi65tuk[i,:] = pnoi65[i,:] * wind 
#
#file_prefix = "fault2_tso_tdo_pno"
#ldataName = ["Ts0", "Td0","Pn0"]
#lgridded_myData = [tsoi65tuk, tdoi65tuk, pnoi65tuk]
#writeNetcdf4SeisSol(file_prefix, grid_x2[0,:], grid_y2[:,1], ldataName, lgridded_myData)
#writeNetcdf4Paraview(file_prefix, grid_x2[0,:], grid_y2[:,1], ldataName, lgridded_myData)    
#
##################
#### SEGMENT 3 ###
##################    
#ua66 = [-0.9407628118057804, 0.339065380012885, 0.0]
#ub66 = [-0.04953169759294571, -0.13742948070747152, 0.9892723329629894]
#ta66 = -1542820.231779111
#tb66 = 248900.26132915495
#
#xa  = np.dot(bary66, ua66) + ta66
#xb  = np.dot(bary66, ub66) + tb66
#xyb66 = np.transpose(np.vstack((xa, xb))) 
#
#step1 = 100
#grid_x3, grid_y3, tsoi66 = interpolate_2dgrid(np.min(xyb66[:,0]),np.max(xyb66[:,0]),np.min(xyb66[:,1]),np.max(xyb66[:,1]),step1,np.c_[xyb66,tso66],'cubic')
#grid_x3, grid_y3, tdoi66 = interpolate_2dgrid(np.min(xyb66[:,0]),np.max(xyb66[:,0]),np.min(xyb66[:,1]),np.max(xyb66[:,1]),step1,np.c_[xyb66,tdo66],'cubic')
#grid_x3, grid_y3, pnoi66 = interpolate_2dgrid(np.min(xyb66[:,0]),np.max(xyb66[:,0]),np.min(xyb66[:,1]),np.max(xyb66[:,1]),step1,np.c_[xyb66,pno66],'cubic')
#
#tsoi66 = ndima.median_filter(tsoi66,size=10)
#tdoi66 = ndima.median_filter(tdoi66,size=30)
#pnoi66 = ndima.median_filter(pnoi66,size=30)
#
#wind = signal.windows.tukey(np.shape(tsoi66)[1],alpha=0.1)
#
#tsoi66tuk = np.zeros(np.shape(tsoi66))
#tdoi66tuk = np.zeros(np.shape(tdoi66))
#pnoi66tuk = np.zeros(np.shape(pnoi66))
#
#for i in range(0, np.shape(tsoi66)[0]):
#    tsoi66tuk[i,:] = tsoi66[i,:] * wind
#    tdoi66tuk[i,:] = tdoi66[i,:] * wind    
#    pnoi66tuk[i,:] = pnoi66[i,:] * wind    
#
# Create the asagi file 
#file_prefix = "fault3_tso_tdo_pno"
#ldataName = ["Ts0", "Td0","Pn0"]
#lgridded_myData = [tsoi66tuk, tdoi66tuk, pnoi66tuk]
#writeNetcdf4SeisSol(file_prefix, grid_x3[0,:], grid_y3[:,1], ldataName, lgridded_myData)
#writeNetcdf4Paraview(file_prefix, grid_x3[0,:], grid_y3[:,1], ldataName, lgridded_myData)  

""
print(np.max(tsoi3tuk), np.max(tdoi3tuk), np.max(pnoi3tuk))

""
###############
### FIGURES ###
###############
fig2, axs = plt.subplots(nrows=1, ncols=1, figsize=(12, 15))
plt.subplots_adjust(top=0.962, bottom=0.059, left=0.029, right=0.992, hspace=0.332, wspace=0.200)

p1 = axs.pcolormesh(grid_x1, grid_y1, tsoi3tuk, vmin=-2e6, vmax=2e6)
axs.set_xlabel('Distance along strike [m]', fontsize=14)
axs.set_ylabel('Distance up-dip [m]', fontsize=14)
axs.set_title('Strike component of the shear stress change', fontsize=16)
axs.set_aspect('equal', 'box')
divider = make_axes_locatable(axs)
cax = divider.append_axes('right', size=0.1, pad=0.05)
clb = fig2.colorbar(p1, cax=cax, orientation='vertical')
clb.ax.set_title("[Pa]", fontsize=14, y=-0.1)


axs.text(.01, .97, 'East', ha='left', va='top', fontsize=14, transform=axs.transAxes, bbox=dict(boxstyle="square", facecolor='white'))
axs.text(.94, .97, 'West', ha='left', va='top', fontsize=14, transform=axs.transAxes, bbox=dict(boxstyle="square", facecolor='white'))

#p2 = axs[1].pcolormesh(grid_x2, grid_y2, tsoi65tuk, vmin=-1e7, vmax=1e7)
#axs[1].set_title('Strike component of the shear stress change (Pa) - Segment 2')
#axs[1].set_aspect('equal', 'box')
#axs[1].set_xlabel('Distance along strike (m)')
#axs[1].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[1])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')
#
#p2 = axs[2].pcolormesh(grid_x3, grid_y3, tsoi66tuk, vmin=-1e7, vmax=1e7)
#axs[2].set_title('Strike component of the shear stress change (Pa) - Segment 3')
#axs[2].set_aspect('equal', 'box')
#axs[2].set_xlabel('Distance along strike (m)')
#axs[2].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[2])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')

##
fig2, axs = plt.subplots(nrows=1, ncols=1, figsize=(12, 15))
plt.subplots_adjust(top=0.962, bottom=0.059, left=0.029, right=0.992, hspace=0.332, wspace=0.200)

p1 = axs.pcolormesh(grid_x1, grid_y1, tdoi3tuk, vmin=-4e6, vmax=4e6)
axs.set_xlabel('Distance along strike [m]', fontsize=14)
axs.set_ylabel('Distance up-dip [m]', fontsize=14)
axs.set_title('Dip component of the shear stress change', fontsize=16)
axs.set_aspect('equal', 'box')
divider = make_axes_locatable(axs)
cax = divider.append_axes('right', size=0.1, pad=0.05)
clb = fig2.colorbar(p1, cax=cax, orientation='vertical')
clb.ax.set_title("[Pa]", fontsize=14, y=-0.1)

#p2 = axs[1].pcolormesh(grid_x2, grid_y2, tdoi65tuk, vmin=-2e6, vmax=2e6)
#axs[1].set_title('Dip component of the shear stress change (Pa) - Segment 3')
#axs[1].set_aspect('equal', 'box')
#axs[1].set_xlabel('Distance along strike (m)')
#axs[1].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[1])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')
#
#p2 = axs[2].pcolormesh(grid_x3, grid_y3, tdoi66tuk, vmin=-2e6, vmax=2e6)
#axs[2].set_title('Dip component of the shear stress change (Pa) - Segment 2')
#axs[2].set_aspect('equal', 'box')
#axs[2].set_xlabel('Distance along strike (m)')
#axs[2].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[2])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')

##

fig2, axs = plt.subplots(nrows=1, ncols=1, figsize=(12, 15))
plt.subplots_adjust(top=0.962, bottom=0.059, left=0.029, right=0.992, hspace=0.332, wspace=0.200)

p1 = axs.pcolormesh(grid_x1, grid_y1, pnoi3tuk,  vmin=-2e6, vmax=2e6)
axs.set_xlabel('Distance along strike [m]', fontsize=14)
axs.set_ylabel('Distance up-dip [m]', fontsize=14)
axs.set_title('Normal stress change', fontsize=16)
axs.set_aspect('equal', 'box')
divider = make_axes_locatable(axs)
cax = divider.append_axes('right', size=0.1, pad=0.05)
clb = fig2.colorbar(p1, cax=cax, orientation='vertical')
clb.ax.set_title("[Pa]", fontsize=14, y=-0.1)

#p2 = axs[1].pcolormesh(grid_x2, grid_y2, pnoi65tuk, vmin=-2e7, vmax=2e7)
#axs[1].set_title('Dip component of the shear stress change (Pa) - Segment 3')
#axs[1].set_aspect('equal', 'box')
#axs[1].set_xlabel('Distance along strike (m)')
#axs[1].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[1])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')
#
#p2 = axs[2].pcolormesh(grid_x3, grid_y3, pnoi66tuk, vmin=-2e7, vmax=2e7)
#axs[2].set_title('Dip component of the shear stress change (Pa) - Segment 2')
#axs[2].set_aspect('equal', 'box')
#axs[2].set_xlabel('Distance along strike (m)')
#axs[2].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[2])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')


""
###############
### FIGURES ###
###############
fig2, axs = plt.subplots(nrows=3, ncols=1, figsize=(12, 12))
#plt.subplots_adjust(top=0.962, bottom=0.059, left=0.029, right=0.992)#, hspace=0.1, wspace=0.200)

p1 = axs[0].pcolormesh(grid_x1, grid_y1, tsoi3tuk, vmin=-2e6, vmax=2e6)
#axs[0].set_xlabel('Distance along strike [m]', fontsize=14)
axs[0].set_ylabel('Distance up-dip [m]', fontsize=14)
axs[0].set_title('Strike component of the shear stress change', fontsize=16)
axs[0].set_aspect('equal', 'box')
divider = make_axes_locatable(axs[0])
cax = divider.append_axes('right', size=0.1, pad=0.05)
clb = fig2.colorbar(p1, cax=cax, orientation='vertical')
clb.ax.set_title("[Pa]", fontsize=14, y=-0.14)


axs[0].text(.01, .96, 'East', ha='left', va='top', fontsize=14, transform=axs[0].transAxes, bbox=dict(boxstyle="square", facecolor='white'))
axs[0].text(.938, .96, 'West', ha='left', va='top', fontsize=14, transform=axs[0].transAxes, bbox=dict(boxstyle="square", facecolor='white'))

#p2 = axs[1].pcolormesh(grid_x2, grid_y2, tsoi65tuk, vmin=-1e7, vmax=1e7)
#axs[1].set_title('Strike component of the shear stress change (Pa) - Segment 2')
#axs[1].set_aspect('equal', 'box')
#axs[1].set_xlabel('Distance along strike (m)')
#axs[1].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[1])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')
#
#p2 = axs[2].pcolormesh(grid_x3, grid_y3, tsoi66tuk, vmin=-1e7, vmax=1e7)
#axs[2].set_title('Strike component of the shear stress change (Pa) - Segment 3')
#axs[2].set_aspect('equal', 'box')
#axs[2].set_xlabel('Distance along strike (m)')
#axs[2].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[2])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')

##
#fig2, axs = plt.subplots(nrows=1, ncols=1, figsize=(12, 15))
#plt.subplots_adjust(top=0.962, bottom=0.059, left=0.029, right=0.992, hspace=0.332, wspace=0.200)

p1 = axs[1].pcolormesh(grid_x1, grid_y1, tdoi3tuk, vmin=-4e6, vmax=4e6)
#axs[1].set_xlabel('Distance along strike [m]', fontsize=14)
axs[1].set_ylabel('Distance up-dip [m]', fontsize=14)
axs[1].set_title('Dip component of the shear stress change', fontsize=16)
axs[1].set_aspect('equal', 'box')
divider = make_axes_locatable(axs[1])
cax = divider.append_axes('right', size=0.1, pad=0.05)
clb = fig2.colorbar(p1, cax=cax, orientation='vertical')
clb.ax.set_title("[Pa]", fontsize=14, y=-0.14)

#p2 = axs[1].pcolormesh(grid_x2, grid_y2, tdoi65tuk, vmin=-2e6, vmax=2e6)
#axs[1].set_title('Dip component of the shear stress change (Pa) - Segment 3')
#axs[1].set_aspect('equal', 'box')
#axs[1].set_xlabel('Distance along strike (m)')
#axs[1].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[1])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')
#
#p2 = axs[2].pcolormesh(grid_x3, grid_y3, tdoi66tuk, vmin=-2e6, vmax=2e6)
#axs[2].set_title('Dip component of the shear stress change (Pa) - Segment 2')
#axs[2].set_aspect('equal', 'box')
#axs[2].set_xlabel('Distance along strike (m)')
#axs[2].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[2])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')

##

#fig2, axs = plt.subplots(nrows=1, ncols=1, figsize=(12, 15))
#plt.subplots_adjust(top=0.962, bottom=0.059, left=0.029, right=0.992, hspace=0.332, wspace=0.200)

p1 = axs[2].pcolormesh(grid_x1, grid_y1, pnoi3tuk,  vmin=-2e6, vmax=2e6)
axs[2].set_xlabel('Distance along strike [m]', fontsize=14)
axs[2].set_ylabel('Distance up-dip [m]', fontsize=14)
axs[2].set_title('Normal stress change', fontsize=16)
axs[2].set_aspect('equal', 'box')
divider = make_axes_locatable(axs[2])
cax = divider.append_axes('right', size=0.1, pad=0.05)
clb = fig2.colorbar(p1, cax=cax, orientation='vertical')
clb.ax.set_title("[Pa]", fontsize=14, y=-0.14)

#p2 = axs[1].pcolormesh(grid_x2, grid_y2, pnoi65tuk, vmin=-2e7, vmax=2e7)
#axs[1].set_title('Dip component of the shear stress change (Pa) - Segment 3')
#axs[1].set_aspect('equal', 'box')
#axs[1].set_xlabel('Distance along strike (m)')
#axs[1].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[1])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')
#
#p2 = axs[2].pcolormesh(grid_x3, grid_y3, pnoi66tuk, vmin=-2e7, vmax=2e7)
#axs[2].set_title('Dip component of the shear stress change (Pa) - Segment 2')
#axs[2].set_aspect('equal', 'box')
#axs[2].set_xlabel('Distance along strike (m)')
#axs[2].set_ylabel('Distance along width (m)')
#divider = make_axes_locatable(axs[2])
#cax = divider.append_axes('right', size=0.1, pad=0.05)
#fig2.colorbar(p1, cax=cax, orientation='vertical')

plt.savefig("Stress_change_{}_{}.png".format(model, version), dpi=300)
#######################
# ## Write yaml file ###
# ######################
write_yaml = True
if write_yaml:
    yamlfile = "{}_fault_stress_change_Samos".format(model)
    hh = [ua3]#, ua65, ua66]
    hw = [ub3]#, ub65, ub66]
    th = [ta3]#, ta65, ta66]
    tw = [tb3]#, tb65, tb66]
    asagifile = [file_prefix]#, "fault2_tso_tdo_pno","fault3_tso_tdo_pno"]
    #    asagifile = ["{}_fault1_tso_tdo_pno".format(model)]#, "fault2_tso_tdo_pno","fault3_tso_tdo_pno"]
    segmentID = ["3"]#, "65", "66"]
    write_yaml_file(yamlfile,asagifile,segmentID,hh,hw,th,tw) 

###############################################################################
# # Compute average stress drop

asl = sx.ReadData('ASl',ndt-1)
print(max(asl))
threshold = 0.1*max(asl)
print(threshold)

""
count = 0
stress_drop_accum = 0

for i in range(0,len(asl)):
    if asl[i] < threshold:
        continue
    else:
        count += 1
        stress_drop_accum += np.sqrt(tdo[i]**2 + tso[i]**2)

""
stress_drop_accum/count

""
###############
### FIGURES ###
###############
fig2, axs = plt.subplots(nrows=1, ncols=1, figsize=(12, 15))
plt.subplots_adjust(top=0.962, bottom=0.059, left=0.029, right=0.992, hspace=0.332, wspace=0.200)

p1 = axs.pcolormesh(grid_x1, grid_y1, tsoi3tuk, vmin=-2e6, vmax=2e6)
axs.set_xlabel('Distance along strike (m)')
axs.set_ylabel('Distance along width (m)')
axs.set_title('Strike component of the shear stress change (Pa) - Segment 1')
axs.set_aspect('equal', 'box')
divider = make_axes_locatable(axs)
cax = divider.append_axes('right', size=0.1, pad=0.05)
fig2.colorbar(p1, cax=cax, orientation='vertical')


fig2, axs = plt.subplots(nrows=1, ncols=1, figsize=(12, 15))
plt.subplots_adjust(top=0.962, bottom=0.059, left=0.029, right=0.992, hspace=0.332, wspace=0.200)

p1 = axs.pcolormesh(grid_x1, grid_y1, tdoi3tuk, vmin=-4e6, vmax=4e6)
axs.set_xlabel('Distance along strike (m)')
axs.set_ylabel('Distance along width (m)')
axs.set_title('Dip component of the shear stress change (Pa) - Segment 1')
axs.set_aspect('equal', 'box')
divider = make_axes_locatable(axs)
cax = divider.append_axes('right', size=0.1, pad=0.05)
fig2.colorbar(p1, cax=cax, orientation='vertical')


fig2, axs = plt.subplots(nrows=1, ncols=1, figsize=(12, 15))
plt.subplots_adjust(top=0.962, bottom=0.059, left=0.029, right=0.992, hspace=0.332, wspace=0.200)

p1 = axs.pcolormesh(grid_x1, grid_y1, pnoi3tuk,  vmin=-2e6, vmax=2e6)
axs.set_xlabel('Distance along strike (m)')
axs.set_ylabel('Distance along width (m)')
axs.set_title('Normal stress change (Pa) - Segment 1')
axs.set_aspect('equal', 'box')
divider = make_axes_locatable(axs)
cax = divider.append_axes('right', size=0.1, pad=0.05)
fig2.colorbar(p1, cax=cax, orientation='vertical')


###############################################################################
# # Get forced_rupture_time and hypocenter
#
# [Plicka et al., (2022)](http://geo.mff.cuni.cz/~gallovic/abst/Plicka_etal.Tectonophysics2022.pdf):
#
# The location of the epicenter is very stable, at 37.900◦N and
# 26.817◦E (±1 km) within all the models tested (Fig. 2 and Table 1).
# Although the depth is generally the most challenging parameter to
# constrain (with error for each velocity model of up to 4 km), all models
# indicate hypocenter depths in the upper crust, at ~12 km and shallower.
#
# Centroid location: 37.90 26.59

from pyproj import Transformer
transformer = Transformer.from_crs("epsg:4326", "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=26.25 +lat_0=37.75", 
                                   always_xy=True)
import pandas as pd

""
lon = [26.817]#, 26.59]
lat =  [37.900]#, 37.90]
x, y = transformer.transform(lon, lat)
print(x,y)

""
# initialize data of lists.
df = {'x1': [x[0], x[0]], 'x2': [y[0], y[0]], 'x3': [0, -9000]}
#df = {'x1': [x[0], x[1]], 'x2': [y[0], y[1]], 'x3': [-6000, -6000]}
df = pd.DataFrame(df)
df

""
df.to_csv("potential_hypo.csv", index=False)

###############################################################################
# z = 6350 ?
