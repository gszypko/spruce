#!/usr/bin/env python3

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.image import NonUniformImage
import sys
import argparse

fullnames = {
  'rho': "Density",
  'temp': "Temperature",
  'press': "Pressure",
  'rad': "Radiative Loss Rate",
  'energy': "Energy Density",
  'b': "Magnetic Field",
  'vel': "Velocity",
  'v': "Velocity",
  'vel_x': "X Velocity",
  'vel_y': "Y Velocity",
  'dt': "CFL Timestep",
  'dt_thermal': "Thermal Conduction Timestep",
  'dt_rad': "Radiative Losses Timestep",
  'n': "Number Density",
  'beta': "Plasma Beta"
}
fullunits = {
  'rho': r'g cm$^{-3}$',
  'temp': r'K',
  'press': r'dyn cm$^{-2}$',
  'rad': r'erg cm$^{-3}$ s$^{-1}$',
  'energy': r'erg cm$^{-3}$',
  'b': r'G',
  'vel': r'cm s$^{-1}$',
  'v': r'cm s$^{-1}$',
  'vel_x': r'cm s$^{-1}$',
  'vel_y': r'cm s$^{-1}$',
  'dt': r's',
  'dt_thermal': r's',
  'dt_rad': r's',
  'n': r'cm$^{-3}$',
  'beta': r''
}
for key in fullunits.keys():
  if key != "beta": fullunits[key] = " (" + fullunits[key] + ")"

parser = argparse.ArgumentParser(description='View the output from mhdtoy.')
parser.add_argument('filename', help="the name of the file output from mhdtoy")
parser.add_argument('timestep', type=float, help="the interval (in simulation time units) between frames of output animation")
parser.add_argument('contourvar', help="the simulation variable to display as a contour plot", choices=['rho', 'temp', 'press', 'rad', 'energy', 'vel_x', 'vel_y', 'dt', 'dt_thermal', 'dt_rad', 'n', 'beta'])
parser.add_argument('-v', '--vector', help="designates vector variable to overlay over contour plot", choices=['b','vel','v'])
parser.add_argument('--density', metavar="vec_display_density", type=int, help="set the interval between displayed vectors", default=25)
parser.add_argument('-d', '--diff', metavar='diff_filename', help="filename to difference with original file")
parser.add_argument('-r', '--realtime', action='store_true')
args = parser.parse_args()

vec_interval = args.density #number of vectors in each direction to display(?)

input_file = open(args.filename)
display_interval = float(args.timestep)
output_var = args.contourvar #must be one of rho, temp, press, energy, rad
file_var_name = output_var
if output_var == "beta":
  file_var_name = "press"
vec_var = args.vector

#determine grid size
assert input_file.readline() == "xdim,ydim\n"
dim = input_file.readline().split(',')
xdim = int(dim[0])
ydim = int(dim[1])

xl_ghost = 0
xu_ghost = 0
yl_ghost = 2
yu_ghost = 2

xdim_view = xdim - (xl_ghost + xu_ghost)
ydim_view = ydim - (yl_ghost + yu_ghost)

xl = 0 + xl_ghost
xu = xdim - xu_ghost
yl = 0 + yl_ghost
yu = ydim - yu_ghost

#read in grid cell positions
assert input_file.readline() == "pos_x\n"
X = np.zeros([xdim_view,ydim_view], dtype=float)
rows = input_file.readline().split(';')[xl:xu]
for i in range(xdim_view):
  X[i] = np.asarray(rows[i].split(','))[yl:yu]
  assert len(X[i]) == ydim_view
assert input_file.readline() == "pos_y\n"
Y = np.zeros([xdim_view,ydim_view], dtype=float)
rows = input_file.readline().split(';')[xl:xu]
for i in range(xdim_view):
  Y[i] = np.asarray(rows[i].split(','))[yl:yu]
  assert len(Y[i]) == ydim_view

x_min = X[0][0]
x_max = X[-1][0]
y_min = Y[0][0]
y_max = Y[0][-1]

#read in static magnetic field
bx = np.zeros([xdim_view,ydim_view], dtype=float)
by = np.zeros([xdim_view,ydim_view], dtype=float)
assert input_file.readline() == "b_x\n"
rows_x = input_file.readline().split(';')
assert input_file.readline() == "b_y\n"
rows_y = input_file.readline().split(';')
assert input_file.readline() == "b_z\n"
input_file.readline()

for i in range(xdim_view):
  bx[i] = np.asarray(rows_x[i].split(','))[yl:yu]
  by[i] = np.asarray(rows_y[i].split(','))[yl:yu]

var = []
vec_x = []
vec_y = []
t = []
output_number = 0
line = ""

while True:
  #read through to next time step (or file end)
  line = input_file.readline()
  while line and line[0:2] != "t=":
    line = input_file.readline()
  if not line:
    break
  time = float(line.split('=')[1])

  #read through to correct variable
  #vector quantities must be last in each time step
  line = input_file.readline()
  while line and line.rstrip() != file_var_name:
    if line[0:2] == "t=":
      sys.exit("Specified output variable not found in file")
    input_file.readline()
    line = input_file.readline()
  if not line:
    break

  if output_number == 0 or display_interval == 0 or ((time - t[0])/display_interval >= output_number and not math.isinf(t[-1])):
    output_number += 1
    t.append(time)
    this_var = np.zeros([xdim_view,ydim_view], dtype=float)
    rows = input_file.readline().split(';')[xl:xu]
    if len(rows) != xdim_view: break
    for i in range(xdim_view):
      this_var_row = np.asarray(rows[i].split(','))[yl:yu]
      if len(this_var_row) != ydim_view: break
      this_var[i] = this_var_row
    var.append(this_var)

    if vec_var != None and vec_var != "b":
      line = input_file.readline()
      while line and line.rstrip().split('_')[0] != vec_var:
        if line[0:2] == "t=":
          sys.exit("Specified output vector not found in file")
        input_file.readline()
        line = input_file.readline()
      if not line:
        break

      if output_number == 1 or display_interval == 0 or (time - t[0])/display_interval >= (output_number-1):
        this_vec_x = np.zeros([xdim_view,ydim_view], dtype=float)
        rows = input_file.readline().split(';')[xl:xu]
        if len(rows) != xdim_view: break
        for i in range(xdim_view):
          row_list = rows[i].split(',')[yl:yu]
          if len(row_list) != ydim_view: break
          this_vec_x[i] = np.asarray(row_list)
        this_vec_y = np.zeros([xdim_view,ydim_view], dtype=float)
        input_file.readline()
        rows = input_file.readline().split(';')[xl:xu]
        if len(rows) != xdim_view: break
        for i in range(xdim_view):
          row_list = rows[i].split(',')[yl:yu]
          if len(row_list) != ydim_view: break
          this_vec_y[i] = np.asarray(row_list)
        vec_x.append(this_vec_x)
        vec_y.append(this_vec_y)
  
input_file.close()

if vec_var != None and len(var) > len(vec_x):
  print("pop!")
  var.pop()

if output_var == "rho":
  for i in range(len(var)):
    var[i] = np.ma.masked_where(var[i]<=1.0e-30, var[i])

if output_var == "beta":
  mag_press = (bx*bx + by*by)/(8.0*np.pi)
  for i in range(len(var)):
    var[i] = var[i]/mag_press

fig, ax = plt.subplots()

frame = 0
x = X[:,0]
y = Y[0,:]
if output_var == "rad":
  # im = ax.imshow(np.transpose(var[frame]), animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
  #   interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-5, base=10))
  # ax.set(xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  # im.set_clim(vmin=1e-4)
  # var_colorbar = fig.colorbar(im)
  # var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-5, base=10))
  im.set_data(x,y,np.transpose(var[frame]))
  im.set_clim(vmin=1e-4)
  ax.add_image(im)
  ax.set(xlim=(x_min,x_max), ylim=(y_min,y_max), xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  var_colorbar = fig.colorbar(im)
  var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
elif output_var == "dt" or output_var == "dt_thermal" or output_var == "dt_rad":
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-10, base=10))
  im.set_data(x,y,np.transpose(var[frame]))
  if(np.isnan(np.nanmax(var))): im.set_clim(vmin=1e-4,vmax=1.0)
  else: im.set_clim(vmin=1e-4,vmax=np.nanmax(var))
  ax.add_image(im)
  ax.set(xlim=(x_min,x_max), ylim=(y_min,y_max), xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  var_colorbar = fig.colorbar(im)
  var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
  # im = ax.imshow(np.transpose(var[frame]), animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
  #   interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-10, base=10))
  # ax.set(xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  # if(np.isnan(np.nanmax(var))): im.set_clim(vmin=1e-4,vmax=1.0)
  # else: im.set_clim(vmin=1e-4,vmax=np.nanmax(var))
  # var_colorbar = fig.colorbar(im)
  # var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
elif output_var == "vel_x" or output_var == "vel_y":
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest')
  im.set_data(x,y,np.transpose(var[frame]))
  ax.add_image(im)
  ax.set(xlim=(x_min,x_max), ylim=(y_min,y_max), xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  var_colorbar = fig.colorbar(im)
  var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
  # im = ax.imshow(np.transpose(var[frame]), animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
  #   interpolation='nearest')
  # ax.set(xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  # var_colorbar = fig.colorbar(im)
  # var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
else:
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.LogNorm())
  im.set_data(x,y,np.transpose(var[frame]))
  ax.add_image(im)
  ax.set(xlim=(x_min,x_max), ylim=(y_min,y_max), xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  var_colorbar = fig.colorbar(im)
  var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
  # im = ax.imshow(np.transpose(var[frame]), animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
  #   interpolation='nearest', norm=matplotlib.colors.LogNorm())
  # ax.set(xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  # var_colorbar = fig.colorbar(im)
  # var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")

contour_color_axes = fig.axes[-1]

# pull out correctly spaced vectors
x_interval = (x_max - x_min)/(vec_interval + 1)
y_interval = (y_max - y_min)/(vec_interval + 1)
x_indices = []
y_indices = []
curr_index_num = 1
for i in range(xdim_view):
  if (X[i,0] - x_min)>curr_index_num*x_interval and (X[i-1,0] - x_min)<curr_index_num*x_interval:
    if abs((X[i,0] - x_min)-curr_index_num*x_interval) < abs((X[i-1,0] - x_min)-curr_index_num*x_interval):
      x_indices.append(i)
    else:
      x_indices.append(i-1)
    curr_index_num+=1
curr_index_num = 1
for j in range(ydim_view):
  if (Y[0,j] - y_min)>curr_index_num*y_interval and (Y[0,j-1] - y_min)<curr_index_num*y_interval:
    if abs((Y[0,j] - y_min)-curr_index_num*y_interval) < abs((Y[0,j-1] - y_min)-curr_index_num*y_interval):
      y_indices.append(j)
    else:
      y_indices.append(j-1)
    curr_index_num+=1
x_vec_indices = []
y_vec_indices = []
for i in range(len(x_indices)):
  for j in range(len(y_indices)):
    x_vec_indices.append(x_indices[i])
    y_vec_indices.append(y_indices[j])
x_vec_indices = np.array(x_vec_indices)
y_vec_indices = np.array(y_vec_indices)

if vec_var == "b":
  this_vec_x = bx
  this_vec_y = by
  norm = np.sqrt(this_vec_x**2 + this_vec_y**2)
  np.divide(this_vec_x, norm, out=this_vec_x, where=norm > 0)
  np.divide(this_vec_y, norm, out=this_vec_y, where=norm > 0)
  quiv = ax.quiver(X[vec_interval::vec_interval, vec_interval::vec_interval], \
    Y[vec_interval::vec_interval, vec_interval::vec_interval], \
      this_vec_x[vec_interval::vec_interval, vec_interval::vec_interval], \
        this_vec_y[vec_interval::vec_interval, vec_interval::vec_interval], \
          norm[vec_interval::vec_interval, vec_interval::vec_interval], \
            cmap=plt.cm.plasma, ec='k', lw=0.2, scale_units='inches', angles='xy', scale=20, \
              pivot='mid', norm=matplotlib.colors.LogNorm())
  vec_colorbar = fig.colorbar(quiv)
  vec_colorbar.set_label(fullnames[vec_var]+ " (" + fullunits[vec_var] + ")")
  vector_color_axes = fig.axes[-1]
  # ax.quiver(X[::vec_interval, ::vec_interval], \
  # Y[::vec_interval, ::vec_interval], \
  #   bx[::vec_interval, ::vec_interval], \
  #     by[::vec_interval, ::vec_interval], pivot='mid')
elif vec_var != None:
  this_vec_x = vec_x[frame].copy()
  this_vec_y = vec_y[frame].copy()
  norm = np.sqrt(this_vec_x**2 + this_vec_y**2)
  np.divide(this_vec_x, norm, out=this_vec_x, where=norm > 0)
  np.divide(this_vec_y, norm, out=this_vec_y, where=norm > 0)
  # quiv = ax.quiver(X[vec_interval::vec_interval, vec_interval::vec_interval], \
  #   Y[vec_interval::vec_interval, vec_interval::vec_interval], \
  #     this_vec_x[vec_interval::vec_interval, vec_interval::vec_interval], \
  #       this_vec_y[vec_interval::vec_interval, vec_interval::vec_interval], \
  #         norm[vec_interval::vec_interval, vec_interval::vec_interval], \
  #           cmap=plt.cm.plasma, ec='k', lw=0.2, scale_units='inches', angles='xy', scale=10, \
  #             width = 0.005, headwidth=3, headlength=3, headaxislength=2.5, \
  #               pivot='mid', norm=matplotlib.colors.SymLogNorm(linthresh=1e4, base=10))
  quiv = ax.quiver(X[x_vec_indices, y_vec_indices], \
    Y[x_vec_indices, y_vec_indices], \
      this_vec_x[x_vec_indices, y_vec_indices], \
        this_vec_y[x_vec_indices, y_vec_indices], \
          norm[x_vec_indices, y_vec_indices], \
            cmap=plt.cm.plasma, ec='k', lw=0.2, scale_units='inches', angles='xy', scale=10, \
              width = 0.005, headwidth=3, headlength=3, headaxislength=2.5, \
                pivot='mid', norm=matplotlib.colors.SymLogNorm(linthresh=1e4, base=10))
  vec_colorbar = fig.colorbar(quiv)
  vec_colorbar.set_label(fullnames[vec_var]+ fullunits[vec_var])
  vector_color_axes = fig.axes[-1]
  # ax.quiverkey(quiv, X=0.3, Y=1.1, U=10,
  #   label='Quiver key, length = 10', labelpos='E')
plt.tight_layout()
frame = -2

def updatefig(*args):
    global frame
    # global im
    global ax
    global quiv
    frame = (frame + 1)%len(var)
    im.set_data(x,y,np.transpose(var[frame]))
    # if output_var != "rho": im.autoscale()
    im.autoscale()
    contour_color_axes.cla()
    var_colorbar = fig.colorbar(im, cax=contour_color_axes)
    var_colorbar.set_label(fullnames[output_var]+ fullunits[output_var])
    ax.set(xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
    if vec_var != "b" and vec_var != None:
      vector_color_axes.cla()
      vec_colorbar = fig.colorbar(quiv, cax=vector_color_axes)
      vec_colorbar.set_label(fullnames[vec_var]+ " (" + fullunits[vec_var] + ")")
      this_vec_x = vec_x[frame].copy()
      this_vec_y = vec_y[frame].copy()
      norm = np.sqrt(this_vec_x**2 + this_vec_y**2)
      np.divide(this_vec_x, norm, out=this_vec_x, where=norm > 0)
      np.divide(this_vec_y, norm, out=this_vec_y, where=norm > 0)
      quiv.set_UVC(this_vec_x[x_vec_indices, y_vec_indices], \
          this_vec_y[x_vec_indices, y_vec_indices], \
            norm[x_vec_indices, y_vec_indices])
      quiv.autoscale()
    return im, ax

ani = animation.FuncAnimation(fig, updatefig, frames=(len(var)), repeat=args.realtime, interval=100, blit=False)

if args.realtime:
  plt.show()
else:
  FFwriter = animation.FFMpegWriter(bitrate=2000)
  filename_suffix = "_"+output_var
  if vec_var != None: filename_suffix = filename_suffix + "_" + vec_var
  ani.save(args.filename.split('/')[-2]+filename_suffix+'.mp4', writer = FFwriter, dpi=100)