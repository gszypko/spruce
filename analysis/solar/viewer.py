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
import aia_response as aia
from viewer_lib import *
import os

fullnames = {
  'rho': "Density",
  'temp': "Temperature",
  'press': "Pressure",
  'rad': "Radiative Loss Rate",
  'thermal_energy': "Energy Density",
  'be': "Background Magnetic Field",
  'bi': "Induced Magnetic Field",
  'bi_x': "X Induced Magnetic Field",
  'bi_y': "Y Induced Magnetic Field",
  'bi_z': "Z Induced Magnetic Field",
  'b': "Magnetic Field",
  'v': "Velocity",
  'v_x': "X Velocity",
  'v_y': "Y Velocity",
  'v_z': "Z Velocity",
  'dt': "CFL Timestep",
  'dt_thermal': "Thermal Conduction Timestep",
  'dt_rad': "Radiative Losses Timestep",
  'n': "Number Density",
  'beta': "Plasma Beta",
  'div_be': "Background Field Divergence",
  'div_bi': "Induced Field Divergence",
  'b_mag': "Field Magnitude",
  'field_heating': "Field-Based Heating Rate"
}
fullunits = {
  'rho': r'g cm$^{-3}$',
  'temp': r'K',
  'press': r'dyn cm$^{-2}$',
  'rad': r'erg cm$^{-3}$ s$^{-1}$',
  'thermal_energy': r'erg cm$^{-3}$',
  'be': r'G',
  'bi': r'G',
  'bi_x': r'G',
  'bi_y': r'G',
  'bi_z': r'G',
  'b': r'G',
  'v': r'cm s$^{-1}$',
  'v_x': r'cm s$^{-1}$',
  'v_y': r'cm s$^{-1}$',
  'v_z': r'cm s$^{-1}$',
  'dt': r's',
  'dt_thermal': r's',
  'dt_rad': r's',
  'n': r'cm$^{-3}$',
  'beta': r'',
  'div_b': r'G cm$^{-1}$',
  'div_bi': r'G cm$^{-1}$',
  'b_mag': r'G',
  'field_heating': r'erg cm$^{-3}$ s$^{-1}$'
}
for filter in aia.filter_names():
  fullnames[filter] = filter.split("A")[0] + " " + r'$\mathrm{\AA}$'
  fullunits[filter] = r'Simulated AIA Image'
for key in fullunits.keys():
  if fullunits[key] != r'':
    fullunits[key] = " (" + fullunits[key] + ")"

parser = argparse.ArgumentParser(description='View the output from mhdtoy.')
parser.add_argument('filename', help="the name of the file output from mhdtoy")
parser.add_argument('timestep', type=float, help="the interval (in simulation time units) between frames of output animation")
parser.add_argument('contourvar', help="the simulation variable to display as a contour plot")
parser.add_argument('-D', '--diff_file', help="plot the difference of the contourvar quantity with a different output file")
parser.add_argument('-v', '--vector', help="designates vector variable to overlay over contour plot", choices=['b','v','be','bi'])
parser.add_argument('-s', '--start', help="designates simulation time to begin plotting",type=float,default=0.0)
parser.add_argument('-e', '--end', help="designates simulation time to end plotting",type=float,default=-1.0)
parser.add_argument('-V', '--vecmode', help="designates mode of display for the chosen vector quantity", choices=['quiver','stream'], default='quiver')
parser.add_argument('-so','--streampoint_override', help="filename containing newline-separated x,y source points for streamline vector drawing", type=str)
parser.add_argument('-sr','--streampoint_record', help="mode for easily selecting streampoints; clicking on plot adds or deletes streamlines; \
                                                        press 0 to save the current set of streamlines for later use", action='store_true')
parser.add_argument('-u','--uniform_stream', help="force uniform streamline spacing for vector field", action='store_true')
parser.add_argument('-sp','--stream_power', help="override default power used to scale vector field in stream mode",type=float)
parser.add_argument('-sn','--stream_number', help="override default number of source points for vector field in stream mode",type=int)
parser.add_argument('--density', metavar="vec_display_density", type=int, help="set the interval between displayed vectors", default=25)
parser.add_argument('-t', '--time_label', help="inserts the time into the plot title" , action='store_true')
parser.add_argument('-tr', '--time_label_rounded', help="when -t is specified, rounds the time to the nearest integer", action='store_true')
parser.add_argument('-d', '--dark_mode', help="renders using a black background and white text", action='store_true')
parser.add_argument('-r', '--realtime', help="render in real time, instead of writing out to a video file", action='store_true')
parser.add_argument('-S', '--serif', help="render using a serif font", action='store_true')
parser.add_argument('-fv', '--full_varname', help="use full variable name in title, instead of abbreviation", action='store_true')
parser.add_argument('-nv', '--no_varname', help="omit variable name from title", action='store_true')
parser.add_argument('-nu', '--no_units', help="omit units from title", action='store_true')
parser.add_argument('-L','--low_colorbar', help="set the lower limit of the contour color bar", type=float)
parser.add_argument('-H','--high_colorbar', help="set the upper limit of the contour color bar", type=float)
parser.add_argument('-cl','--colorbar_location',help="set location of colorbar relative to plot (left, right, bottom, or top)",type=str,default='right')
parser.add_argument('-fs', '--font_size', help="set the font size", type=float, default=10)
parser.add_argument('-g', '--ghost_zones', help="includes ghost zones in the contour plot", action='store_true')
parser.add_argument('-nc', '--no_colorbar', help="omit colorbar for contour quantity", action='store_true')
parser.add_argument('-nt','--no_ticks', help="omits labels along the x- and y-axes", action='store_true')
parser.add_argument('-i','--interactive',help="manually control the time stepping of the plot \
                    (with left and right arrow keys) and investigate contour values by mousing over",action='store_true')
args = parser.parse_args()

matplotlib.rcParams.update({'font.size': args.font_size})
if args.dark_mode: plt.style.use('dark_background')
if args.serif: matplotlib.rcParams.update({'font.family': 'serif'})

vec_interval = args.density #number of vectors in each direction to display(?)
vec_mode = args.vecmode

start_time = args.start
end_time = args.end

#construct path to text file containing streamline plotting start points
outpath = args.filename.split("/")[:-1]
startpointfile = ""
for folder in outpath: startpointfile += folder + "/"
startpointfile += f"startpoints_{args.vector}.txt"

display_interval = float(args.timestep)
output_var = args.contourvar

# Define any contour quantities that require multiple file values
file_var_names = np.array([output_var])
if output_var == "beta":
  # file_var_names = np.array(["press","bi_x","bi_y","bi_z"])
  file_var_names = np.array(["press","bi_x","bi_y"])
if output_var == "b_mag":
  # file_var_names = np.array(["bi_x","bi_y","bi_z"])
  file_var_names = np.array(["bi_x","bi_y"])
if output_var in ["div_b","div_bi"]:
  file_var_names = np.array(["bi_x","bi_y"])
if aia.is_filter(output_var):
  file_var_names = np.array(["n","temp"])
vec_var = args.vector
file_vec_name = vec_var
if vec_var == "b":
  file_vec_name = "bi"

num_ghost = 2
if args.ghost_zones: num_ghost = 0

xl_ghost = num_ghost
xu_ghost = num_ghost
yl_ghost = num_ghost
yu_ghost = num_ghost

xdim, ydim, X, Y, file_vars, vec_x, vec_y, bx, by, bz, t =\
  extract_data_from_file(args.filename, file_var_names, display_interval, xl_ghost, xu_ghost, yl_ghost, yu_ghost, start_time, end_time, file_vec_name)

if args.diff_file is not None:
  xdim_diff, ydim_diff, foo, foo, file_vars_diff, foo, foo, foo, foo, foo, t_diff =\
    extract_data_from_file(args.diff_file, file_var_names, display_interval, xl_ghost, xu_ghost, yl_ghost, yu_ghost, start_time, end_time)
  assert xdim_diff == xdim and ydim_diff == ydim
  file_vars_interp = [ [] for v in file_var_names ]
  for i_base in range(len(t)):
    t_base = t[i_base]
    i_diff = np.searchsorted(t_diff,t_base)
    i_diff_lower = min(len(t_diff)-1,max(0,i_diff))
    i_diff_upper = max(0,min(len(t_diff)-1,i_diff))
    for name_idx in range(len(file_var_names)):
      val_lower = file_vars_diff[name_idx][i_diff_lower]
      val_upper = file_vars_diff[name_idx][i_diff_upper]
      if i_diff_lower == i_diff_upper:
        interp = val_lower
      else:
        t_diff_upper = t[i_diff_upper]
        t_diff_lower = t[i_diff_lower]
        interp = (val_lower*(t_diff_upper - t_base) + val_upper*(t_base - t_diff_lower))/(t_diff_upper - t_diff_lower)
      file_vars_interp[name_idx].append(interp)
  file_vars_diff = file_vars_interp

xl = 0 + xl_ghost
xu = xdim - xu_ghost
yl = 0 + yl_ghost
yu = ydim - yu_ghost
xdim_view = xdim - (xl_ghost + xu_ghost)
ydim_view = ydim - (yl_ghost + yu_ghost)

x = X[:,0]
y = Y[0,:]

x_min = X[0][0]
x_max = X[-1][0]
y_min = Y[0][0]
y_max = Y[0][-1]

# Apply any post-reading calculations for contourvars that combine multiple variables
var = apply_contour_computation(output_var, file_vars, x, y, bx, by, bz)
if args.diff_file is not None:
  var_diff = apply_contour_computation(output_var, file_vars_diff, x, y, bx, by, bz)
  for i in range(len(var)):
    var[i] = var[i] - var_diff[i]

if vec_var == "b":
  for i in range(len(vec_x)):
    vec_x[i] = bx + vec_x[i]
    vec_y[i] = by + vec_y[i]

if vec_var != None and vec_var != "be" and len(var) > len(vec_x):
  print("pop!")
  var.pop()

# if output_var in ["rho","b_mag"]:
#   for i in range(len(var)):
#     var[i] = np.ma.masked_where(var[i]<=1.0e-30, var[i])

fig, ax = plt.subplots()
frame = 0

if vec_mode == "stream":
  avg_space_x = np.abs(np.mean([x[i+1] - x[i] for i in range(len(x)-1)]))
  x_eq = np.asarray([x[0] + i*avg_space_x for i in range(len(x))])
  avg_space_y = np.abs(np.mean([y[i+1] - y[i] for i in range(len(y)-1)]))
  y_eq = np.asarray([y[0] + i*avg_space_y for i in range(len(y))])

max_v = np.finfo(np.float_).min
min_v = np.finfo(np.float_).max
for i in range(len(var)):
  curr_max = np.nanmax(var[i])
  if np.isfinite(curr_max): max_v = np.fmax(max_v,curr_max)
  curr_min = np.nanmin(var[i])
  if np.isfinite(curr_min): min_v = np.fmin(min_v,curr_min)

if args.low_colorbar is not None:
  min_v = args.low_colorbar
if args.high_colorbar is not None:
  max_v = args.high_colorbar

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
      if np.isfinite(curr_max): max_v_vec = np.fmax(max_v_vec,curr_max)
      curr_min = np.nanmin(magnitude,where=magnitude>0,initial=min_v_vec)
      if np.isfinite(curr_min): min_v_vec = np.fmin(min_v_vec,curr_min)
  if max_v_vec == np.finfo(np.float_).min:
    max_v_vec = 1.0
  if min_v_vec == np.finfo(np.float_).max:
    min_v_vec = 0.0

# Set contour bar scaling (linear, log, etc)
if args.diff_file is not None or output_var in ["div_bi","div_be"]:
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-20, base=10))
elif output_var in ["rad","dt","dt_thermal","dt_rad","field_heating","b_mag"]: #variables that can go to zero
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-8, base=10))
elif aia.is_filter(output_var):
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest')
else: # output_var in ["v_x","v_y","v_z","bi_x","bi_y"]: #variables with negative and positive values
  im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
    interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-2, base=10))
# else:
#   im = NonUniformImage(ax, animated=True, origin='lower', extent=(x_min,x_max,y_min,y_max),\
#     interpolation='nearest', norm=matplotlib.colors.LogNorm())

# Construct custom colorbar
# if args.diff_file is not None:
#   colors = [(0.0,(0.1,0.1,1.0)),(0.5,(1.0,1.0,1.0)),(1.0,(1.0,0.1,0.1))]
#   div_cmap = LinearSegmentedColormap.from_list('light_rwb',colors)
#   im.set_cmap(div_cmap)
if output_var == "beta" and args.diff_file is None:
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
elif aia.is_filter(output_var):
  im.set_cmap(aia.colormap(output_var))

im.set_data(x,y,np.transpose(var[frame]))
im.set_clim(vmin=min_v,vmax=max_v)
ax.add_image(im)
if args.colorbar_location in ["bottom","top"]:
  ax.set_aspect('equal')
ax.set(xlim=(x_min,x_max), ylim=(y_min,y_max), title=output_var+", t="+str(t[frame]))
if args.no_ticks:
  ax.set(xticks=[],yticks=[])
else:
  ax.set(xlabel="x (Mm)",ylabel="y (Mm)")
if not args.no_colorbar:
  var_colorbar = fig.colorbar(im,location=args.colorbar_location)

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
    if os.path.exists(startpointfile) or args.streampoint_override is not None:
      start_pts=[]
      input_filename = ""
      if args.streampoint_override is not None:
        input_filename = args.streampoint_override
      else:
        input_filename = startpointfile
      with open(input_filename) as f:
        line = f.readline()
        while line:
          if line[0] != "#" and line != '':
            start_pts.append(line.split(","))
            start_pts[-1] = [float(s) for s in start_pts[-1]]
          line = f.readline()
    else:
      if args.stream_number is None:
        num_points = 11
      else:
        num_points = args.stream_number
      if args.uniform_stream:
        stream_points = np.linspace(x[1],x[-2],num_points)
      else:
        if args.stream_power is None:
          if this_vec_y[1,1]*this_vec_y[-2,1] < 0.0:
            print("loop detected")
            num_points = int(num_points*2)
            stream_pow = 1
          else:
            stream_pow = 2
        else:
          stream_pow = args.stream_power
        norm_normalized = norm[:,0]**stream_pow/np.trapz(norm[:,0]**stream_pow)
        cdf = [np.trapz(norm_normalized[:(i+1)]) for i in range(len(norm_normalized))]
        stream_points = np.interp(np.linspace(0.05,0.95,num=num_points),cdf,x)
      start_pts=np.column_stack((stream_points,y_eq[0]*np.ones_like(stream_points)))
    start_pts = np.asarray(start_pts)
    print(start_pts)
    stream = ax.streamplot(x_eq,y_eq,this_vec_x.transpose(),this_vec_y.transpose(),\
      start_points=start_pts,\
        color=(0.0,0.0,0.0,1.0),density=(xdim_view/10,ydim_view/10),linewidth=2.0,arrowstyle='->',arrowsize=1.5,maxlength=10.0)
    stream.lines.set(pickradius=2.0)
  else:
    assert vec_mode == "quiver"
    quiv = ax.quiver(X[x_vec_indices, y_vec_indices], Y[x_vec_indices, y_vec_indices], \
        this_vec_x[x_vec_indices, y_vec_indices], this_vec_y[x_vec_indices, y_vec_indices], norm[x_vec_indices, y_vec_indices], \
              cmap=plt.cm.Greys, ec='k', lw=0.4, scale_units='inches', angles='xy', scale=10, \
                width = 0.005, headwidth=3, headlength=3, headaxislength=2.5, \
                  pivot='mid', norm=matplotlib.colors.SymLogNorm(linthresh=1e-4, base=10))
    vec_colorbar = fig.colorbar(quiv)
    vec_colorbar.set_label(fullnames[vec_var]+ fullunits.get(vec_var,""))
    vector_color_axes = fig.axes[-1]
    quiv.set_clim(vmin=min_v_vec,vmax=max_v_vec)

plt.tight_layout()
frame = -2

def updatefig(*args_arg):
    global frame
    global ax
    global quiv
    global stream
    global start_pts
    if not (args.streampoint_record or args.interactive):
      frame = (frame + 1)%len(var)
    im.set_data(x,y,np.transpose(var[frame]))
    plot_title = ""
    if args.diff_file is not None:
      plot_title += "Difference in "
    if not args.no_varname:
      if args.full_varname: plot_title += fullnames[output_var]
      else: plot_title += output_var
    if not args.no_units:
      plot_title += fullunits.get(output_var,"")
    if args.time_label:
      if len(plot_title) != 0: plot_title += ", "
      if args.time_label_rounded: plot_title += "t="+str(round(t[frame]))+" s"
      else: plot_title += "t="+str(t[frame])+" s"
    ax.set(title=plot_title)
    if not args.no_ticks:
      ax.set(xlabel="x (Mm)",ylabel="y (Mm)")
    if vec_var != "be" and vec_var != None:
      this_vec_x = vec_x[frame].copy()
      this_vec_y = vec_y[frame].copy()
      if vec_mode == "quiver":
        vector_color_axes.cla()
        vec_colorbar = fig.colorbar(quiv, cax=vector_color_axes)
        vec_colorbar.set_label(fullnames[vec_var]+ fullunits.get(vec_var,""))
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
          start_points=start_pts,\
            color=(0.0,0.0,0.0,1.0),density=(xdim_view/10,ydim_view/10),linewidth=2.0,arrowstyle='->',arrowsize=1.5,maxlength=10.0)
        if not args.streampoint_record:
          print("frame " + str(frame) + " plotted")
    plt.tight_layout()
    return im, ax

# Argument: a LineCollection
# Returns: lines, line_segments
# where lines is a list of numpy Mx2 arrays detailing vertices in each discrete line from the LineCollection
# and where line_segments contains the same, except each line is represented by a list of 2x2 arrays indicating individual segments
def extract_separate_lines(segments):
  lines = []
  line_segments = []
  i=0
  while i<len(segments):
    this_line=[segments[i][0]]
    this_line_segments=[segments[i]]
    i+=1
    while i<len(segments) and np.max((abs(segments[i][0,0]-segments[i-1][1,0]),abs(segments[i][0,1]-segments[i-1][1,1])))<1.0:
      this_line.append(segments[i][0])
      this_line_segments.append(segments[i])
      i+= 1
    lines.append(np.asarray(this_line))
    line_segments.append(this_line_segments)
  return lines, line_segments

def modifystreampoints(event):
  ix, iy = event.xdata, event.ydata
  print(f"{ix},{iy} clicked")
  if ix is not None and iy is not None and args.streampoint_record:
    global start_pts
    global stream
    global fig
    global ax
    if stream.lines.contains(event)[0]:
      segments = stream.lines.get_segments()
      lines, line_segments = extract_separate_lines(segments)
      for l in range(len(lines)):
        single_line = matplotlib.collections.LineCollection(line_segments[l],pickradius=5.0,transform=ax.transData)
        if single_line.contains(event)[0]:
          print("Deleting existing line")
          start_pts = np.delete(start_pts,l,axis=0)
    else:
      print("Adding new line")
      start_pts = np.append(start_pts,[[ix, iy]],axis=0)
  updatefig()
  fig.canvas.draw_idle()

def savestreampoints(event):
  if(event.key == "0"):
    global stream
    global args
    global startpointfile
    segments = stream.lines.get_segments()
    lines = extract_separate_lines(segments)[0]
    print(f"Saving {len(lines)} lines to {startpointfile}...")
    with open(startpointfile, 'w') as f:
      for l in lines:
        f.write(f"{l[0][0]},{l[0][1]}\n")

def advancetime(event):
  global frame
  if event.key == "right":
    frame = (frame + 1)%len(var)
    updatefig()
    fig.canvas.draw_idle()
  elif event.key == "left":
    frame = (frame - 1)%len(var)
    updatefig()
    fig.canvas.draw_idle()

def keypresshandler(event):
  if args.interactive:
    advancetime(event)
  if args.streampoint_record:
    savestreampoints(event)

def displayinfo(event):
  global fig
  global x,y
  global var
  global frame
  if event.xdata is None or event.ydata is None:
    fig.suptitle("")
  else:
    xidx = min(np.searchsorted(x,event.xdata),xdim-1)
    if xidx != 0:
      xidx = xidx if (event.xdata-x[xidx-1])>(x[xidx]-event.xdata) else xidx-1
    yidx = min(np.searchsorted(y,event.ydata),ydim-1)
    if yidx != 0:
      yidx = yidx if (event.ydata-y[yidx-1])>(y[yidx]-event.ydata) else yidx-1
    fig.suptitle(f"{var[frame][xidx,yidx]:.3e} at ({xidx},{yidx})",x=0.1,ha='left')
  fig.canvas.draw_idle()


if args.interactive or args.streampoint_record:
  cid = fig.canvas.mpl_connect('key_press_event', keypresshandler)
if args.interactive:
  cid = fig.canvas.mpl_connect('motion_notify_event', displayinfo)
if args.streampoint_record:
  cid = fig.canvas.mpl_connect('button_press_event', modifystreampoints)
else:
  ani = animation.FuncAnimation(fig, updatefig, frames=len(var), repeat=args.realtime, interval=100, blit=False)

if args.realtime or args.streampoint_record or args.interactive:
  plt.show()
else:
  FFwriter = animation.FFMpegWriter(bitrate=2000)
  filename_suffix = "_"+output_var
  if vec_var != None: filename_suffix = filename_suffix + "_" + vec_var
  ani.save(args.filename.split('/')[-2]+filename_suffix+'.mp4', writer = FFwriter, dpi=100)