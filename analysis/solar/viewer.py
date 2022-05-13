#!/usr/bin/env python3

from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.image import NonUniformImage
import sys
import argparse

# matplotlib.rcParams.update({'font.size': 22})

def read_grid(in_file,xdim,ydim,xl,yl,xu,yu):
  result = []
  for i in range(xdim):
    row = in_file.readline().split(',')[yl:yu]
    if i >= xl and i < xu:
      result.append(row)
  return np.asarray(result).astype(np.float64)

fullnames = {
  'rho': "Density",
  'temp': "Temperature",
  'press': "Pressure",
  'rad': "Radiative Loss Rate",
  'thermal_energy': "Energy Density",
  'be': "Background Magnetic Field",
  'bi': "Induced Magnetic Field",
  'v': "Velocity",
  'v_x': "X Velocity",
  'v_y': "Y Velocity",
  'dt': "CFL Timestep",
  'dt_thermal': "Thermal Conduction Timestep",
  'dt_rad': "Radiative Losses Timestep",
  'n': "Number Density",
  'beta': "Plasma Beta",
  'div_be': "Background Field Divergence",
  'div_bi': "Induced Field Divergence",
  'b_mag': "Field Magnitude"
}
fullunits = {
  'rho': r'g cm$^{-3}$',
  'temp': r'K',
  'press': r'dyn cm$^{-2}$',
  'rad': r'erg cm$^{-3}$ s$^{-1}$',
  'thermal_energy': r'erg cm$^{-3}$',
  'be': r'G',
  'bi': r'G',
  'v': r'cm s$^{-1}$',
  'v_x': r'cm s$^{-1}$',
  'v_y': r'cm s$^{-1}$',
  'dt': r's',
  'dt_thermal': r's',
  'dt_rad': r's',
  'n': r'cm$^{-3}$',
  'beta': r'',
  'div_be': r'G cm$^{-1}$',
  'div_bi': r'G cm$^{-1}$',
  'b_mag': r'G'
}
for key in fullunits.keys():
  if fullunits[key] != r'':
    fullunits[key] = " (" + fullunits[key] + ")"

parser = argparse.ArgumentParser(description='View the output from mhdtoy.')
parser.add_argument('filename', help="the name of the file output from mhdtoy")
parser.add_argument('timestep', type=float, help="the interval (in simulation time units) between frames of output animation")
parser.add_argument('contourvar', help="the simulation variable to display as a contour plot", choices=['rho', 'temp', 'press', 'rad', \
  'thermal_energy', 'v_x', 'v_y', 'dt', 'dt_thermal', 'dt_rad', 'n', 'beta', 'div_be', 'div_bi', 'b_mag'])
parser.add_argument('-v', '--vector', help="designates vector variable to overlay over contour plot", choices=['be','v','bi'])
parser.add_argument('-V', '--vecmode', help="designates mode of display for the chosen vector quantity", choices=['quiver','stream'], default='quiver')
parser.add_argument('--density', metavar="vec_display_density", type=int, help="set the interval between displayed vectors", default=25)
# parser.add_argument('-d', '--diff', metavar='diff_filename', help="filename to difference with original file")
parser.add_argument('-r', '--realtime', action='store_true')
args = parser.parse_args()

vec_interval = args.density #number of vectors in each direction to display(?)
vec_mode = args.vecmode

input_file = open(args.filename)
display_interval = float(args.timestep)
output_var = args.contourvar
file_var_names = np.array([output_var])
if output_var == "beta":
  file_var_names = np.array(["press","bi_x","bi_y"])
if output_var == "b_mag":
  file_var_names = np.array(["bi_x","bi_y"])
vec_var = args.vector

line = input_file.readline()
while line[0] == '#' or len(line) == 0:
  line = input_file.readline()
#determine grid size
assert line == "xdim,ydim\n"
dim = input_file.readline().split(',')
xdim = int(dim[0])
ydim = int(dim[1])

xl_ghost = 2
xu_ghost = 2
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
X = read_grid(input_file,xdim,ydim,xl,yl,xu,yu)
assert input_file.readline() == "pos_y\n"
Y = read_grid(input_file,xdim,ydim,xl,yl,xu,yu)

x_min = X[0][0]
x_max = X[-1][0]
y_min = Y[0][0]
y_max = Y[0][-1]

assert input_file.readline() == "be_x\n"
bx = read_grid(input_file,xdim,ydim,xl,yl,xu,yu)
assert input_file.readline() == "be_y\n"
by = read_grid(input_file,xdim,ydim,xl,yl,xu,yu)

file_vars = [ [] for v in file_var_names ]
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

  # var_found = False
  vars_found = np.array([ False for v in file_var_names ]) #all set to false initially
  vec_x_found = (vec_var == None or vec_var == "be")
  vec_y_found = (vec_var == None or vec_var == "be")
  # collect all necessary data at the current time step
  while not (np.all(vars_found) and vec_x_found and vec_y_found):
    line = input_file.readline()
    while line and not (np.any(line.rstrip() == file_var_names) or line.rstrip() == (str(vec_var)+"_x") or line.rstrip() == (str(vec_var)+"_y")):
      if line[0:2] == "t=" and not np.all(vars_found):
        sys.exit("Specified output variable(s) not found in file")
      if line[0:2] == "t=" and not (vec_x_found and vec_y_found):
        sys.exit("Specified output vector not found in file")
      line = input_file.readline()
    if not line:
      break

    curr_data = []
    if np.any(line.rstrip() == file_var_names): 
      curr_data.append("var")
      var_location = np.argwhere(line.rstrip() == file_var_names)
      assert len(var_location) == 1
      var_index = var_location[0][0]
      assert vars_found[var_index] == False
      vars_found[var_index] = True
    if line.rstrip() == (str(vec_var)+"_x"): 
      curr_data.append("vec_x")
      vec_x_found = True
    if line.rstrip() == (str(vec_var)+"_y"): 
      curr_data.append("vec_y")
      vec_y_found = True
    
    if output_number == 0 or display_interval == 0 or ((time - t[0])/display_interval >= output_number and not math.isinf(t[-1])):
      if (np.all(vars_found) and vec_x_found and vec_y_found):
        output_number += 1
        t.append(time)
      this_var = read_grid(input_file,xdim,ydim,xl,yl,xu,yu)
      if "var" in curr_data:
        file_vars[var_index].append(this_var)
      if "vec_x" in curr_data:
        vec_x.append(this_var)
      if "vec_y" in curr_data:
        vec_y.append(this_var)

  if not line:
    break

input_file.close()

var = []
if output_var == "beta":
  for i in range(len(file_vars[0])):
    mag_press = ((bx+file_vars[1][i])**2 + (by+file_vars[2][i])**2)/(8.0*np.pi)
    var.append(file_vars[0][i]/mag_press)
elif output_var == "b_mag":
  for i in range(len(file_vars[0])):
    field_mag = np.sqrt((bx+file_vars[0][i])**2 + (by+file_vars[1][i])**2)
    var.append(field_mag)
else:
  var = file_vars[0]

if vec_var != None and vec_var != "be" and len(var) > len(vec_x):
  print("pop!")
  var.pop()

if output_var == "rho":
  for i in range(len(var)):
    var[i] = np.ma.masked_where(var[i]<=1.0e-30, var[i])

fig, ax = plt.subplots()
# if (vec_var != None and vec_var != "bi" and vec_var != "be") or output_var == "temp":
#   fig.set_figwidth(7.0/6.0*fig.get_size_inches()[0])

frame = 0
x = X[:,0]
y = Y[0,:]

if vec_mode == "stream":
  avg_space_x = np.abs(np.mean([x[i+1] - x[i] for i in range(len(x)-1)]))
  x_eq = np.asarray([x[0] + i*avg_space_x for i in range(len(x))])
  avg_space_y = np.abs(np.mean([y[i+1] - y[i] for i in range(len(y)-1)]))
  y_eq = np.asarray([y[0] + i*avg_space_y for i in range(len(y))])

max_v = np.finfo(np.float_).min
min_v = np.finfo(np.float_).max
for i in range(len(var)):
  curr_max = np.nanmax(var[i][xl:xu,yl:yu])
  max_v = np.fmax(max_v,curr_max)
  curr_min = np.nanmin(var[i][xl:xu,yl:yu])
  min_v = np.fmin(min_v,curr_min)

if vec_var != None:
  max_v_vec = np.finfo(np.float_).min
  min_v_vec = np.finfo(np.float_).max
  if vec_var == "be":
    magnitude = np.sqrt(bx[xl:xu,yl:yu]**2 + by[xl:xu,yl:yu]**2)
    max_v_vec = np.fmax(np.nanmax(magnitude,where=magnitude>0,initial=max_v_vec),max_v_vec)
    min_v_vec = np.fmin(np.nanmin(magnitude,where=magnitude>0,initial=min_v_vec),min_v_vec)
  else:
    for i in range(len(var)):
      magnitude = np.sqrt(vec_x[i][xl:xu,yl:yu]**2 + vec_y[i][xl:xu,yl:yu]**2)
      curr_max = np.nanmax(magnitude,where=magnitude>0,initial=max_v_vec)
      max_v_vec = np.fmax(max_v_vec,curr_max)
      curr_min = np.nanmin(magnitude,where=magnitude>0,initial=min_v_vec)
      min_v_vec = np.fmin(min_v_vec,curr_min)
  if max_v_vec == np.finfo(np.float_).min:
    max_v_vec = 1.0
  if min_v_vec == np.finfo(np.float_).max:
    min_v_vec = 0.0

if output_var in ["rad","beta","dt","dt_thermal","dt_rad"]: #variables that can go to zero
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-8, base=10))
elif output_var in ["div_bi","div_be"]:
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-22, base=10))
elif output_var in ["vel_x","vel_y"]: #variables with negative and positive values
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-2, base=10))
else:
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.LogNorm())

if output_var == "beta":
  logrange = np.log10(max_v) - np.log10(min_v)
  one_pos = (np.log10(1.0) - np.log10(min_v))/logrange
  small_range = min(abs(one_pos),abs(1.0-one_pos))
  if np.log10(max_v)*np.log10(min_v) < 0.0:
    colors = [(0.0,(0.0,0.1,0.4)),(one_pos-0.5*small_range,(0.0,0.2,0.9)),(one_pos,(1.0,1.0,1.0)),(one_pos+0.5*small_range,(0.9,0.2,0.0)),(1.0,(0.4,0.1,0.0))]
  elif np.log10(max_v) < 0.0:
    middle = 0.5*one_pos
    top_color = [np.interp(1.0,[middle,one_pos],[0.0,1.0]), np.interp(1.0,[middle,one_pos],[0.2,1.0]), np.interp(1.0,[middle,one_pos],[0.9,1.0])]
    colors = [(0.0,(0.0,0.1,0.4)),(min(middle,1.0),(0.0,0.2,0.9)),(1.0,top_color)]
  elif np.log10(max_v) > 0.0:
    middle = 1.0 - 0.5*(1.0-one_pos)
    bott_color = [np.interp(0.0,[one_pos,middle],[1.0,0.9]), np.interp(0.0,[one_pos,middle],[1.0,0.2]), np.interp(0.0,[one_pos,middle],[1.0,0.0])]
    colors = [(0.0,bott_color),(max(middle,0.0),(0.9,0.2,0.0)),(1.0,(0.4,0.1,0.0))]
  beta_cmap = LinearSegmentedColormap.from_list('beta_diverge',colors)
  im.set_cmap(beta_cmap)
elif output_var in ["div_bi","div_be"]:
  colors = [(0.0,(0.5,0.5,1.0)),(0.5,(1.0,1.0,1.0)),(1.0,(1.0,0.5,0.5))]
  div_cmap = LinearSegmentedColormap.from_list('light_rwb',colors)
  im.set_cmap(div_cmap)

im.set_data(x,y,np.transpose(var[frame]))
im.set_clim(vmin=min_v,vmax=max_v)
ax.add_image(im)
ax.set(xlim=(x_min,x_max), ylim=(y_min,y_max), xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame]))
# ax.set_aspect('equal')
var_colorbar = fig.colorbar(im)
var_colorbar.set_label(fullnames[output_var]+ fullunits[output_var])

# contour_color_axes = fig.axes[-1]

if vec_mode == "quiver":
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

if vec_var != None:
  if vec_var == "be":
    this_vec_x = bx
    this_vec_y = by
  else:
    this_vec_x = vec_x[frame].copy()
    this_vec_y = vec_y[frame].copy()
  norm = np.sqrt(this_vec_x**2 + this_vec_y**2)
  np.divide(this_vec_x, norm, out=this_vec_x, where=norm > 0)
  np.divide(this_vec_y, norm, out=this_vec_y, where=norm > 0)
  if vec_mode == "stream":
    num_points = 17
    if this_vec_x[1,1]*this_vec_x[-2,1] > 0.0:
      print("loop detected")
      num_points = int(num_points*2)
      stream_pow = 1
    else:
      stream_pow = 2
    norm_normalized = norm[:,0]**stream_pow/np.trapz(norm[:,0]**stream_pow)
    cdf = [np.trapz(norm_normalized[:(i+1)]) for i in range(len(norm_normalized))]
    stream_points = np.interp(np.linspace(0.05,0.95,num=num_points),cdf,x)
    print(stream_points)
    stream = ax.streamplot(x_eq,y_eq,this_vec_x.transpose(),this_vec_y.transpose(),\
      start_points=np.column_stack((stream_points,y_eq[0]*np.ones_like(stream_points))),\
        color=(0.0,0.0,0.0),density=100,linewidth=1.0,arrowstyle='->',maxlength=1000.0)
  else:
    assert vec_mode == "quiver"
    quiv = ax.quiver(X[x_vec_indices, y_vec_indices], Y[x_vec_indices, y_vec_indices], \
        this_vec_x[x_vec_indices, y_vec_indices], this_vec_y[x_vec_indices, y_vec_indices], norm[x_vec_indices, y_vec_indices], \
              cmap=plt.cm.Greys, ec='k', lw=0.4, scale_units='inches', angles='xy', scale=10, \
                width = 0.005, headwidth=3, headlength=3, headaxislength=2.5, \
                  pivot='mid', norm=matplotlib.colors.SymLogNorm(linthresh=1e-4, base=10))
    vec_colorbar = fig.colorbar(quiv)
    vec_colorbar.set_label(fullnames[vec_var]+ fullunits[vec_var])
    vector_color_axes = fig.axes[-1]
    quiv.set_clim(vmin=min_v_vec,vmax=max_v_vec)

plt.tight_layout()
frame = -2

def updatefig(*args):
    global frame
    global ax
    global quiv
    global stream
    frame = (frame + 1)%len(var)
    im.set_data(x,y,np.transpose(var[frame]))
    ax.set(xlabel="x (cm)",ylabel="y (cm)",title=output_var+", t="+str(t[frame])+" s")
    if vec_var != "be" and vec_var != None:
      this_vec_x = vec_x[frame].copy()
      this_vec_y = vec_y[frame].copy()
      if vec_mode == "quiver":
        vector_color_axes.cla()
        vec_colorbar = fig.colorbar(quiv, cax=vector_color_axes)
        vec_colorbar.set_label(fullnames[vec_var]+ fullunits[vec_var])
        norm = np.sqrt(this_vec_x**2 + this_vec_y**2)
        np.divide(this_vec_x, norm, out=this_vec_x, where=norm > 0)
        np.divide(this_vec_y, norm, out=this_vec_y, where=norm > 0)
        quiv.set_UVC(this_vec_x[x_vec_indices, y_vec_indices], \
            this_vec_y[x_vec_indices, y_vec_indices], \
              norm[x_vec_indices, y_vec_indices])
      else:
        assert vec_mode == "stream"
        stream.lines.remove()
        for art in ax.get_children():
          if not isinstance(art, matplotlib.patches.FancyArrowPatch):
            continue
          art.remove()
        stream = ax.streamplot(x_eq,y_eq,this_vec_x.transpose(),this_vec_y.transpose(),\
          start_points=np.column_stack((stream_points,y_eq[0]*np.ones_like(stream_points))),\
            color=(0.0,0.0,0.0),density=100,linewidth=1.0,arrowstyle='->',maxlength=1000.0)
    plt.tight_layout()
    return im, ax

ani = animation.FuncAnimation(fig, updatefig, frames=(len(var)), repeat=args.realtime, interval=100, blit=False)

if args.realtime:
  plt.show()
else:
  FFwriter = animation.FFMpegWriter(bitrate=2000)
  filename_suffix = "_"+output_var
  if vec_var != None: filename_suffix = filename_suffix + "_" + vec_var
  ani.save(args.filename.split('/')[-2]+filename_suffix+'.mp4', writer = FFwriter, dpi=100)