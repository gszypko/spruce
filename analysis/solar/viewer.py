#!/usr/bin/env python3

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
  'dt_rad': "Radiative Losses Timestep"
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
  'dt_rad': r's'
}

parser = argparse.ArgumentParser(description='View the output from mhdtoy.')
parser.add_argument('filename', help="the name of the file output from mhdtoy")
parser.add_argument('timestep', type=float, help="the interval (in simulation time units) between frames of output animation")
parser.add_argument('contourvar', help="the simulation variable to display as a contour plot", choices=['rho', 'temp', 'press', 'rad', 'energy', 'vel_x', 'vel_y', 'dt', 'dt_thermal', 'dt_rad'])
parser.add_argument('-v', '--vector', help="designates vector variable to overlay over contour plot", choices=['b','vel','v'])
parser.add_argument('--density', metavar="vec_display_density", type=int, help="set the interval between displayed vectors", default=5)
parser.add_argument('-d', '--diff', metavar='diff_filename', help="filename to difference with original file")
parser.add_argument('-r', '--realtime', action='store_true')
args = parser.parse_args()

vec_interval = args.density #interval between B field vectors to display

input_file = open(args.filename)
display_interval = float(args.timestep)
output_var = args.contourvar #must be one of rho, temp, press, energy, rad
vec_var = args.vector

#determine grid size
assert input_file.readline() == "xdim,ydim\n"
dim = input_file.readline().split(',')
xdim = int(dim[0])
ydim = int(dim[1])

#read in grid cell positions
assert input_file.readline() == "pos_x\n"
X = np.zeros([xdim,ydim], dtype=float)
rows = input_file.readline().split(';')
for i in range(xdim):
  X[i] = np.asarray(rows[i].split(','))
  assert len(X[i]) == ydim
assert input_file.readline() == "pos_y\n"
Y = np.zeros([xdim,ydim], dtype=float)
rows = input_file.readline().split(';')
for i in range(xdim):
  Y[i] = np.asarray(rows[i].split(','))
  assert len(Y[i]) == ydim

x_min = X[0][0]
x_max = X[-1][0]
y_min = Y[0][0]
y_max = Y[0][-1]

#read in static magnetic field
bx = np.zeros([xdim,ydim], dtype=float)
by = np.zeros([xdim,ydim], dtype=float)
assert input_file.readline() == "b_x\n"
rows_x = input_file.readline().split(';')
assert input_file.readline() == "b_y\n"
rows_y = input_file.readline().split(';')
assert input_file.readline() == "b_z\n"
input_file.readline()

for i in range(xdim):
  bx[i] = np.asarray(rows_x[i].split(','))
  by[i] = np.asarray(rows_y[i].split(','))

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
  while line and line.rstrip() != output_var:
    if line[0:2] == "t=":
      sys.exit("Specified output variable not found in file")
    input_file.readline()
    line = input_file.readline()
  if not line:
    break

  if output_number == 0 or display_interval == 0 or ((time - t[0])/display_interval >= output_number and not math.isinf(t[-1])):
    output_number += 1
    t.append(time)
    this_var = np.zeros([xdim,ydim], dtype=float)
    rows = input_file.readline().split(';')
    if len(rows) != xdim: break
    for i in range(xdim):
      this_var[i] = np.asarray(rows[i].split(','))
      if len(this_var[i]) != ydim: break
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
        this_vec_x = np.zeros([xdim,ydim], dtype=float)
        rows = input_file.readline().split(';')
        if len(rows) != xdim: break
        for i in range(xdim):
          row_list = rows[i].split(',')
          if len(row_list) != ydim: break
          this_vec_x[i] = np.asarray(row_list)
        this_vec_y = np.zeros([xdim,ydim], dtype=float)
        input_file.readline()
        rows = input_file.readline().split(';')
        if len(rows) != xdim: break
        for i in range(xdim):
          row_list = rows[i].split(',')
          if len(row_list) != ydim: break
          this_vec_y[i] = np.asarray(row_list)
        vec_x.append(this_vec_x)
        vec_y.append(this_vec_y)
  
input_file.close()

if vec_var != None and len(var) > len(vec_x):
  print("pop!")
  var.pop()

fig, ax = plt.subplots()

frame = 0
if output_var == "rad":
  im = ax.imshow(np.transpose(var[frame]), animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-5, base=10))
  ax.set(xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  im.set_clim(vmin=1e-4)
  var_colorbar = fig.colorbar(im)
  var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
elif output_var == "dt" or output_var == "dt_thermal" or output_var == "dt_rad":
  im = ax.imshow(np.transpose(var[frame]), animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-10, base=10))
  ax.set(xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  if(np.isnan(np.nanmax(var))): im.set_clim(vmin=1e-4,vmax=1.0)
  else: im.set_clim(vmin=1e-4,vmax=np.nanmax(var))
  var_colorbar = fig.colorbar(im)
  var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
elif output_var == "vel_x" or output_var == "vel_y":
  im = ax.imshow(np.transpose(var[frame]), animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest')
  ax.set(xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  var_colorbar = fig.colorbar(im)
  var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
else:
  im = ax.imshow(np.transpose(var[frame]), animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.LogNorm())
  ax.set(xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
  var_colorbar = fig.colorbar(im)
  var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")

contour_color_axes = fig.axes[-1]

if vec_var == "b":
  this_vec_x = bx
  this_vec_y = by
  norm = np.sqrt(this_vec_x**2 + this_vec_y**2)
  np.divide(this_vec_x, norm, out=this_vec_x, where=norm > 0)
  np.divide(this_vec_y, norm, out=this_vec_y, where=norm > 0)
  quiv = ax.quiver(X[::vec_interval, ::vec_interval], \
    Y[::vec_interval, ::vec_interval], \
      this_vec_x[::vec_interval, ::vec_interval], \
        this_vec_y[::vec_interval, ::vec_interval], \
          norm[::vec_interval, ::vec_interval], \
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
  quiv = ax.quiver(X[::vec_interval, ::vec_interval], \
    Y[::vec_interval, ::vec_interval], \
      this_vec_x[::vec_interval, ::vec_interval], \
        this_vec_y[::vec_interval, ::vec_interval], \
          norm[::vec_interval, ::vec_interval], \
            cmap=plt.cm.plasma, ec='k', lw=0.2, scale_units='inches', angles='xy', scale=6, \
              width = 0.008, headwidth=3, headlength=3, headaxislength=2.5, \
                pivot='mid', norm=matplotlib.colors.SymLogNorm(linthresh=1e4, base=10))
  vec_colorbar = fig.colorbar(quiv)
  vec_colorbar.set_label(fullnames[vec_var]+ " (" + fullunits[vec_var] + ")")
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
    im.set_data(np.transpose(var[frame]))
    im.autoscale()
    contour_color_axes.cla()
    var_colorbar = fig.colorbar(im, cax=contour_color_axes)
    var_colorbar.set_label(fullnames[output_var]+ " (" + fullunits[output_var] + ")")
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
      quiv.set_UVC(this_vec_x[::vec_interval, ::vec_interval], \
          this_vec_y[::vec_interval, ::vec_interval], \
            norm[::vec_interval, ::vec_interval])
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