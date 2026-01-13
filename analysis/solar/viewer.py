#!/usr/bin/env python

from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.image import NonUniformImage
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import LinearNDInterpolator
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from sunpy.map.sources.source_type import source_stretch
import argparse
import aia_response as aia
from viewer_lib import *
import os
from collections import OrderedDict
import bisect

import astropy.units as u
from sunpy.visualization.colormaps import color_tables as ct

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
parser.add_argument('-tpsp','--tracer_particle_streampoint', help="use local tpout file for time-dependent streamline start points", action='store_true')
parser.add_argument('-sr','--streampoint_record', help="mode for easily selecting streampoints; clicking on plot adds or deletes streamlines; \
                                                        press 0 to save the current set of streamlines for later use", action='store_true')
parser.add_argument('-tpr','--tracerparticle_record', help="mode for easily selecting tracer particles; clicking on plot adds or deletes streamlines; \
                                                        press 0 to save the current set of particles to a tpstate file", action='store_true')
parser.add_argument('-u','--uniform_stream', help="force uniform streamline spacing for vector field", action='store_true')
parser.add_argument('-sp','--stream_power', help="override default power used to scale vector field in stream mode",type=float)
parser.add_argument('-sn','--stream_number', help="override default number of source points for vector field in stream mode",type=int)
parser.add_argument('--density', metavar="vec_display_density", type=int, help="set the interval between displayed vectors", default=25)
parser.add_argument('-t', '--time_label', help="inserts the time into the plot title" , action='store_true')
parser.add_argument('-tr', '--time_label_rounded', help="when -t is specified, rounds the time to the nearest integer", action='store_true')
parser.add_argument('-tz', '--time_label_zero', help="when -t is specified, uses the first frame as t=0", action='store_true')
parser.add_argument('-w', '--white_mode', help="renders using a white background and black text", action='store_true', default=False)
parser.add_argument('-d', '--dark_mode', help="nonfunctional option, left in for compatibility", action='store_true')
parser.add_argument('-r', '--realtime', help="render in real time, instead of writing out to a video file", action='store_true')
parser.add_argument('-sa', '--sandbox', help="plot according to the sandbox area in the script", action='store_true')
parser.add_argument('-S', '--serif', help="render using a serif font", action='store_true')
parser.add_argument('-fv', '--full_varname', help="use full variable name in title, instead of abbreviation", action='store_true')
parser.add_argument('-nv', '--no_varname', help="omit variable name from title", action='store_true')
parser.add_argument('-nu', '--no_units', help="omit units from title", action='store_true')
parser.add_argument('-it', '--internal_title', help="places title on interior of plot, instead of above", action='store_true')
parser.add_argument('-L','--low_colorbar', help="set the lower limit of the contour color bar", type=float)
parser.add_argument('-H','--high_colorbar', help="set the upper limit of the contour color bar", type=float)
parser.add_argument('-Lv','--low_colorbar_vec', help="set the lower limit of the quiver vector color bar", type=float)
parser.add_argument('-Hv','--high_colorbar_vec', help="set the upper limit of the quiver vector color bar", type=float)
parser.add_argument('-fsx','--fig_size_x', help="set the horizontal size of the image (inches)", type=float, default=6.4)
parser.add_argument('-fsy','--fig_size_y', help="set the vertical size of the image (inches)", type=float, default=4.8)
parser.add_argument('-mxl','--margin_x_lower', help="set the fractional left margin", type=float, default=-1.0)
parser.add_argument('-mxu','--margin_x_upper', help="set the fractional right margin", type=float, default=-1.0)
parser.add_argument('-myl','--margin_y_lower', help="set the fractional bottom margin", type=float, default=-1.0)
parser.add_argument('-myu','--margin_y_upper', help="set the fractional top margin", type=float, default=-1.0)
parser.add_argument('-cl','--colorbar_location',help="set location of colorbar relative to plot (left, right, bottom, or top)",type=str,default='right')
parser.add_argument('-fs', '--font_size', help="set the font size", type=float, default=10)
parser.add_argument('-to', '--tex_off', help="disables TeX rendering for text", action='store_true')
parser.add_argument('-dpi', '--dpi', help="set the dpi", type=float, default=-1.0)
parser.add_argument('-g', '--ghost_zones', help="includes ghost zones in the contour plot", action='store_true')
parser.add_argument('-gxl','--ghost_xl', help="number of x-lower ghost cells to exclude", type=int, default=2)
parser.add_argument('-gxu','--ghost_xu', help="number of x-upper ghost cells to exclude", type=int, default=2)
parser.add_argument('-gyl','--ghost_yl', help="number of y-lower ghost cells to exclude", type=int, default=2)
parser.add_argument('-gyu','--ghost_yu', help="number of y-upper ghost cells to exclude", type=int, default=2)
parser.add_argument('-nc', '--no_colorbar', help="omit colorbar for contour quantity", action='store_true')
parser.add_argument('-nt','--no_ticks', help="omits labels along the x- and y-axes", action='store_true')
parser.add_argument('-nty','--no_ticks_y', help="omits labels along the y-axis", action='store_true')
parser.add_argument('-ntx','--no_ticks_x', help="omits labels along the x-axis", action='store_true')
parser.add_argument('-nax','--no_axes', help="removes axes, ticks, and labels", action='store_true')
parser.add_argument('-nm','--no_margin', help="removes margins from edges of figure", action='store_true')
parser.add_argument('-i','--interactive',help="manually control the time stepping of the plot \
                    (with left and right arrow keys) and investigate contour values by mousing over",action='store_true')
parser.add_argument('-E','--equal_aspect',help="enforce an equal aspect ratio on the plotted data",action='store_true')
parser.add_argument('-tl', '--tick_list',help="provide comma-separated (no whitespace) list of tick values for colorbar", type=str, default="")
parser.add_argument('-C', '--contour_list',help="provide comma-separated (no whitespace) list of values to include as contour lines", type=str, default="")
parser.add_argument('-te', '--tick_e', help="use abbreviated e notation for colorbar", action="store_true")
parser.add_argument('-slt', '--symlogthreshold', help="override default lower threshold for symlog colorbar scaling", type=float, default=-1.0)
parser.add_argument('-tp', '--tracerparticles', help="display tracer particles", action='store_true')
parser.add_argument('-ltp', '--label_tracerparticles', help="display tracer particle labels", action='store_true')
parser.add_argument('-tps', '--tracerparticles_size', help="set size (in pt) of displayed tracer particles", type=float, default=20)
parser.add_argument('-tpl', '--tracerparticle_list',help="provide comma-separated (no whitespace) list of particle labels to include. Accepts asterisk wildcard", type=str, default="")
parser.add_argument('-tpm', '--tracerparticle_measure',help="provide single particle label to provide measurement for, incl location and contour value at particle location", type=str, default="")
parser.add_argument('-tpts', '--tracerparticle_timeseries',help="provide single particle label to generate time series for", type=str, default="")
parser.add_argument('-chts', '--crosshair_timeseries',help="provide physical x,y coords (example: x50_y100) to generate time series for", type=str, default="")
parser.add_argument('-o', '--out_filename',help="override default naming convention for video file output", type=str, default="")
parser.add_argument('-bt','--backtrace',help="provide physical x_y coords (example: x50_y100) from which to backtrace particles to t=0, comma-separate to do multiple; or, 'custom' for custom trajectory", type=str, default="")
parser.add_argument('-bts', '--backtrace_spacing', help="set time spacing between when new backtraces begin (by default, one per time step)", type=float, default=-1)
parser.add_argument('-bti', '--backtrace_init', help="set label index of first backtraced particle", type=int, default=0)
parser.add_argument('-ps','--particle_synthesize',help="provide physical x,y coords (example: x50_y100) from which to generate multi-particle synthesis; or, 'custom' for custom trajectory", type=str, default="")
parser.add_argument('-pss', '--particle_synthesize_single', help="when -ps is specified, generates particle synthesis instructions using only nearest particle", action='store_true')
parser.add_argument('-ch','--crosshair',help="provide physical x,y coords (example: x50_y100) or comma-separated list thereof from which to generate multi-particle synthesis; or, 'custom' for custom trajectory", type=str, default="")
parser.add_argument('-chs', '--crosshair_size', help="set size (in pt) of displayed crosshair", type=float, default=100)
parser.add_argument('-chl', '--crosshair_labels', help="comma separated list of labels to mark crosshairs with", type=str, default="")
parser.add_argument('-sfs','--still_frames', help="comma-separated list of the indices (based on the time step) to stack into a row of stills", type=str, default="")
parser.add_argument('-cbl','--colorbar_label',help="label contour colorbar with quantity",action='store_true')
parser.add_argument('-stl','--small_ticklabels',help="make tick labels for axes small",action='store_true')

args = parser.parse_args()

line_thickness = 0.5
arrow_size = 1.0

if args.dpi != -1.0:
  matplotlib.rcParams.update({'figure.dpi': args.dpi})
# matplotlib.rcParams.update({'figure.autolayout': True})
matplotlib.rcParams.update({'font.size': args.font_size})
if not args.white_mode: plt.style.use('dark_background')
if args.serif: matplotlib.rcParams.update({'font.family': 'serif'})
if args.contour_list != "":
    contour_vals = [float(v) for v in args.contour_list.split(',')]
if args.tracerparticle_list != "":
    tracerparticle_labels = args.tracerparticle_list.split(',')
else:
    tracerparticle_labels = []

# matplotlib.rcParams.update({'figure.dpi': 300})
# matplotlib.rcParams.update({'font.size': 10})
if not args.tex_off: matplotlib.rcParams.update({'text.usetex': True, 'font.family': 'serif'})

vec_interval = args.density #number of vectors in each direction to display(?)
vec_mode = args.vecmode

start_time = args.start
end_time = args.end

#construct path to text file containing streamline plotting start points
outpath = args.filename.split("/")[:-1]
startpointfile = ""
for folder in outpath: startpointfile += folder + "/"
startpointfile += f"startpoints_{args.vector}.txt"
tracerparticlefile = ""
for folder in outpath: tracerparticlefile += folder + "/"
tracerparticlefile += "particles.tpout"
tracerrecordfile = ""
for folder in outpath: tracerrecordfile += folder + "/"
tracerrecordfile += "tracerpoints.tpstate"

display_interval = float(args.timestep)
output_var = args.contourvar

# Define any contour quantities that require multiple file values
file_var_names = np.array([output_var])
for k in file_vars_dict.keys():
  if output_var in k:
    file_var_names = np.array(file_vars_dict[k])
    break

if args.sandbox:
  file_var_names = np.array(["press","bi_x","bi_y","n"])
if aia.is_filter(output_var):
  file_var_names = np.array(["n","temp"])
vec_var = args.vector
file_vec_name = vec_var
if vec_var == "b":
  file_vec_name = "bi"

xl_ghost = 0 if args.ghost_zones else args.ghost_xl
xu_ghost = 0 if args.ghost_zones else args.ghost_xu
yl_ghost = 0 if args.ghost_zones else args.ghost_yl
yu_ghost = 0 if args.ghost_zones else args.ghost_yu

xdim, ydim, X, Y, file_vars, vec_x, vec_y, bx, by, bz, t =\
  extract_data_from_file(args.filename, file_var_names, display_interval, xl_ghost, xu_ghost, yl_ghost, yu_ghost, start_time, end_time, file_vec_name)

if args.tracerparticles or args.particle_synthesize != "" or args.tracerparticle_timeseries != "":
  # print(tracerparticle_labels)
  tracer_particles, tracer_labels = extract_tracer_particles_from_file(tracerparticlefile, display_interval, start_time, end_time, False, tracerparticle_labels)
  if args.tracerparticle_record:
    print("tp and tpr options not designed to be used at once. Use tpr to create a tpstate file and tp to plot existing particles over time.")
    exit()

if args.tracerparticle_record:
  tracer_particles = ([],[])
  tracer_labels = []
  click_count = 0
  if args.tracerparticles:
    print("tp and tpr options not designed to be used at once. Use tpr to create a tpstate file and tp to plot existing particles over time.")
    exit()

if args.diff_file is not None:
  xdim_diff, ydim_diff, foo, foo, file_vars_diff, foo, foo, foo, foo, foo, t_diff =\
    extract_data_from_file(args.diff_file, file_var_names, display_interval, xl_ghost, xu_ghost, yl_ghost, yu_ghost, start_time, end_time)
  assert xdim_diff == xdim and ydim_diff == ydim
  file_vars_interp = [ [] for v in file_var_names ]
  for i_base in range(len(t)):
    t_base = t[i_base]
    i_diff = np.searchsorted(t_diff,t_base)
    i_diff_lower = min(len(t_diff)-1,max(0,i_diff-1))
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
    var[i] = (var[i] - var_diff[i])/var[i]

if vec_var == "b":
  for i in range(len(vec_x)):
    vec_x[i] = bx + vec_x[i]
    vec_y[i] = by + vec_y[i]

if vec_var != None and vec_var != "be" and len(var) > len(vec_x):
  print("pop!")
  var.pop()

if vec_mode == "stream":
  min_space_x = np.abs(np.min([x[i+1] - x[i] for i in range(len(x)-1)]))
  x_eq = np.asarray([x[0] + i*min_space_x for i in range(math.ceil((x_max - x_min)/min_space_x))])
  min_space_y = np.abs(np.min([y[i+1] - y[i] for i in range(len(y)-1)]))
  y_eq = np.asarray([y[0] + i*min_space_y for i in range(math.ceil((y_max - y_min)/min_space_y))])

# if output_var in ["rho","b_mag"]:
#   for i in range(len(var)):
#     var[i] = np.ma.masked_where(var[i]<=1.0e-30, var[i])

if args.backtrace != "":
  if vec_var != "v":
    print("Backtracing requires the vector variable v to be selected: -v v")
    exit()
  if not args.backtrace_spacing >= args.timestep:
    print("Backtrace spacing must be >= the time step")
    exit()

  print("Generating interpolators...")
  vel_interp_x = RegularGridInterpolator((t,x.copy(order='C'),y.copy(order='C')),np.array(vec_x)/1e8,fill_value=0.0,bounds_error=False)
  vel_interp_y = RegularGridInterpolator((t,x.copy(order='C'),y.copy(order='C')),np.array(vec_y)/1e8,fill_value=0.0,bounds_error=False)

  total_msnum = args.backtrace_init

  if args.backtrace == "custom":
    print("Using custom probe trajectory (not yet implemented)...")
    exit()
  else:
    bt_locs = args.backtrace.split(',')

    for bt_loc in bt_locs:
      if len(bt_loc.split('_')) != 2 or bt_loc.split('_')[0][0] != "x" or bt_loc.split('_')[1][0] != "y":
        print("Expected backtrace point syntax is underscore-separated, with x and y coordinates prepended by x and y respectively. Example: x50_y100")
        exit()
      probe_x = float(bt_loc.split('_')[0][1:])
      probe_y = float(bt_loc.split('_')[1][1:])
      traced_particles = []

      print("Allocating traced particle memory...")
      for i in range(len(t)):
        traced_particles.append(([],[]))

      print("Tracing...")
      if args.backtrace_spacing < 0.0: interval = args.timestep
      else: interval =  args.backtrace_spacing
      num_ples = int((t[-1] - t[0]) // interval) + 1

      for ple in range(num_ples):
        print("ple ",ple)
        ple_t_idx = np.argmin(np.abs(np.array(t)-(ple*interval)))
        traced_particles[ple_t_idx][0].append(probe_x)
        traced_particles[ple_t_idx][1].append(probe_y)
        # Forward trace
        for t_j in range(ple_t_idx+1,len(t)):
          old_pos_x = traced_particles[t_j-1][0][ple]
          old_pos_y = traced_particles[t_j-1][1][ple]
          old_time = t[t_j-1]
          vel_x = vel_interp_x([old_time,old_pos_x,old_pos_y])[0]
          vel_y = vel_interp_y([old_time,old_pos_x,old_pos_y])[0]
          delta_t = t[t_j] - old_time
          traced_particles[t_j][0].append(old_pos_x + delta_t*vel_x)
          traced_particles[t_j][1].append(old_pos_y + delta_t*vel_y)
        # Backward trace
        for t_j in range(ple_t_idx-1,-1,-1):
          old_pos_x = traced_particles[t_j+1][0][ple]
          old_pos_y = traced_particles[t_j+1][1][ple]
          old_time = t[t_j+1]
          vel_x = vel_interp_x([old_time,old_pos_x,old_pos_y])[0]
          vel_y = vel_interp_y([old_time,old_pos_x,old_pos_y])[0]
          delta_t = t[t_j] - old_time
          traced_particles[t_j][0].append(old_pos_x + delta_t*vel_x)
          traced_particles[t_j][1].append(old_pos_y + delta_t*vel_y)
      bt_outpath = '/'.join(args.filename.split('/')[:-1])+f"/backtraced_x{probe_x:.0f}_y{probe_y:.0f}.tpstate"
      print("Backtrace complete, writing starting positions to "+bt_outpath+"...")
      with open(bt_outpath, 'w') as f:
        for ple in range(num_ples):
          f.write(f"{traced_particles[0][0][ple]*1.0e8},{traced_particles[0][1][ple]*1.0e8}#ms{total_msnum}")
          total_msnum += 1
          if ple != num_ples - 1: f.write('\n')
    exit()

if args.particle_synthesize != "":
  if args.particle_synthesize == "custom":
    print("Using custom probe trajectory for synthesis instructions (not yet implemented)...")
    exit()
  else:
    if len(args.particle_synthesize.split('_')) != 2 or args.particle_synthesize.split('_')[0][0] != "x" or args.particle_synthesize.split('_')[1][0] != "y":
      print("Expected synthesize point syntax is underscore-separated, with x and y coordinates prepended by x and y respectively. Example: x50_y100")
      exit()
    probe_x = float(args.particle_synthesize.split('_')[0][1:])
    probe_y = float(args.particle_synthesize.split('_')[1][1:])
    ps_outpath = '/'.join(args.filename.split('/')[:-1])+f"/instructions_x{probe_x:.0f}_y{probe_y:.0f}.plesyn"
    print("Generating particle synthesis instructions at "+ps_outpath+"...")
    max_ple_dist = 20.0
    with open(ps_outpath, 'w') as f:
      for t_i in range(len(t)):
        curr_time = t[t_i]
        curr_closest_ple = -1
        curr_closest_dist = np.inf
        curr_secondclosest_ple = -1
        curr_secondclosest_dist = np.inf
        ples = []
        dists = []
        for ple in range(len(tracer_particles[t_i][0])):
          curr_dist = np.sqrt((tracer_particles[t_i][0][ple] - probe_x)**2 + (tracer_particles[t_i][1][ple] - probe_y)**2)
          if "ms" in tracer_labels[t_i][ple] and (curr_dist < max_ple_dist):
            ples.append(ple)
            dists.append(curr_dist)
            # if curr_dist < curr_closest_dist:
            #   curr_closest_ple = ple
            #   curr_closest_dist = curr_dist
        zipped = sorted(zip(dists, ples))
        ples, dists = [t[1] for t in zipped], [t[0] for t in zipped]
        if len(ples)>=1:
          curr_closest_ple = ples[0]
          curr_closest_dist = dists[0]
          curr_secondclosest_ple = curr_closest_ple
          curr_secondclosest_dist = curr_closest_dist
          if len(ples)>1 and not args.particle_synthesize_single:
            for i in range(1,len(ples)):
              # check that the currently considered particle is not on the *other side* of the closest particle (to the probe)
              closest_to_curr_x = (tracer_particles[t_i][0][ples[i]] - tracer_particles[t_i][0][ples[0]])
              closest_to_curr_y = (tracer_particles[t_i][1][ples[i]] - tracer_particles[t_i][1][ples[0]])
              closest_to_probe_x = (probe_x - tracer_particles[t_i][0][ples[0]])
              closest_to_probe_y = (probe_y - tracer_particles[t_i][1][ples[0]])
              closest_to_curr_mag = np.sqrt(closest_to_curr_x**2 + closest_to_curr_y**2)
              if closest_to_curr_mag == 0.0:
                continue
              projected_distance = closest_to_curr_x/closest_to_curr_mag*closest_to_probe_x + closest_to_curr_y/closest_to_curr_mag*closest_to_probe_y
              if projected_distance > 0.0:
                curr_secondclosest_ple = ples[i]
                curr_secondclosest_dist = closest_to_curr_mag - projected_distance
                curr_closest_dist = projected_distance
                break
        f.write(f"t={t[t_i]},{tracer_labels[t_i][curr_closest_ple]},{curr_closest_dist}")
        if curr_closest_ple != curr_secondclosest_ple: f.write(f",{tracer_labels[t_i][curr_secondclosest_ple]},{curr_secondclosest_dist}")
        if t_i != len(t)-1: f.write('\n')
    exit()
    

# ################## SANDBOX AREA ###################
if args.sandbox:
  print("Starting sandbox mode")

  assert vec_var == "v"
  traced_particles = []
  print("Generating interpolators...")
  vel_interp_x = RegularGridInterpolator((t,x.copy(order='C'),y.copy(order='C')),np.array(vec_x)/1e8,fill_value=0.0,bounds_error=False)
  vel_interp_y = RegularGridInterpolator((t,x.copy(order='C'),y.copy(order='C')),np.array(vec_y)/1e8,fill_value=0.0,bounds_error=False)
  probe_x = -50
  probe_y = 350

  print("Allocating traced particle memory...")

  for i in range(len(t)):
    traced_particles.append(([],[]))
  print("Tracing...")
  for ple in range(len(t)):
    print("ple ",ple)
    traced_particles[ple][0].append(probe_x)
    traced_particles[ple][1].append(probe_y)
    # Forward trace
    for t_j in range(ple+1,len(t)):
      old_pos_x = traced_particles[t_j-1][0][ple]
      old_pos_y = traced_particles[t_j-1][1][ple]
      old_time = t[t_j-1]
      vel_x = vel_interp_x([old_time,old_pos_x,old_pos_y])[0]
      vel_y = vel_interp_y([old_time,old_pos_x,old_pos_y])[0]
      delta_t = t[t_j] - old_time
      traced_particles[t_j][0].append(old_pos_x + delta_t*vel_x)
      traced_particles[t_j][1].append(old_pos_y + delta_t*vel_y)
    # Backward trace
    for t_j in range(ple-1,-1,-1):
      old_pos_x = traced_particles[t_j+1][0][ple]
      old_pos_y = traced_particles[t_j+1][1][ple]
      old_time = t[t_j+1]
      vel_x = vel_interp_x([old_time,old_pos_x,old_pos_y])[0]
      vel_y = vel_interp_y([old_time,old_pos_x,old_pos_y])[0]
      delta_t = t[t_j] - old_time
      traced_particles[t_j][0].append(old_pos_x + delta_t*vel_x)
      traced_particles[t_j][1].append(old_pos_y + delta_t*vel_y)
  # ple_start_t_idx = []
  # for i in range(len(vec_x)):
  #   traced_particles.append(([],[]))
  # start_t_idx = 0
  # for i in range(len(tracer_particles[0][0])):
  #   # CHANGE THIS NEXT
  #   # traced_particles[start_t_idx][0].append(tracer_particles[0][0][i])
  #   # traced_particles[start_t_idx][1].append(tracer_particles[0][1][i])
  #   ple_start_t_idx.append(start_t_idx)
  #   start_t_idx = (start_t_idx + 1)%len(t)
  # print("Tracing...")
  # for ple in range(len(tracer_particles[0][0])):
  #   print("ple ",ple)
  #   traced_particles[ple_start_t_idx[ple]][0].append(tracer_particles[ple_start_t_idx[ple]][0][ple])
  #   traced_particles[ple_start_t_idx[ple]][1].append(tracer_particles[ple_start_t_idx[ple]][1][ple])
  #   # Forward trace
  #   for t_j in range(ple_start_t_idx[ple]+1,len(t)):
  #     old_pos_x = traced_particles[t_j-1][0][ple]
  #     old_pos_y = traced_particles[t_j-1][1][ple]
  #     old_time = t[t_j-1]
  #     vel_x = vel_interp_x([old_time,old_pos_x,old_pos_y])[0]
  #     vel_y = vel_interp_y([old_time,old_pos_x,old_pos_y])[0]
  #     delta_t = t[t_j] - old_time
  #     traced_particles[t_j][0].append(old_pos_x + delta_t*vel_x)
  #     traced_particles[t_j][1].append(old_pos_y + delta_t*vel_y)
  #   # Backward trace
  #   for t_j in range(ple_start_t_idx[ple]-1,-1,-1):
  #     old_pos_x = traced_particles[t_j+1][0][ple]
  #     old_pos_y = traced_particles[t_j+1][1][ple]
  #     old_time = t[t_j+1]
  #     vel_x = vel_interp_x([old_time,old_pos_x,old_pos_y])[0]
  #     vel_y = vel_interp_y([old_time,old_pos_x,old_pos_y])[0]
  #     delta_t = t[t_j] - old_time
  #     traced_particles[t_j][0].append(old_pos_x + delta_t*vel_x)
  #     traced_particles[t_j][1].append(old_pos_y + delta_t*vel_y)
  
  with open('/'.join(args.filename.split('/')[:-1])+"/backtraced.tpstate", 'w') as f:
    for ple in range(len(t)):
      f.write(f"{traced_particles[0][0][ple]*1.0e8},{traced_particles[0][1][ple]*1.0e8}#ms{ple}")
      if ple != len(t) - 1: f.write('\n')
  
  tracer_particles = traced_particles
  tracer_labels = [[f"t{i}" for i in range(len(traced_particles))] for j in range(len(t))]
  # tracer_labels = [["" for i in range(len(traced_particles))] for j in range(len(t))]
  vec_var = None

  with open('/'.join(args.filename.split('/')[:-1])+"/instructions.plesyn", 'w') as f:
    for t_i in range(len(t)):
      curr_time = t[t_i]
      curr_closest_ple = 0
      curr_closest_dist = np.inf
      for ple in range(len(t)):
        curr_dist = np.sqrt((tracer_particles[t_i][0][ple] - probe_x)**2 + (tracer_particles[t_i][1][ple] - probe_y)**2)
        if curr_dist < curr_closest_dist:
          curr_closest_ple = ple
          curr_closest_dist = curr_dist
      f.write(f"t={t[t_i]},{tracer_labels[t_i][curr_closest_ple]},{curr_closest_dist}")
      if t_i != len(t)-1: f.write('\n')
  # tracer_particles[frame][0][i]/1e8,tracer_particles[frame][1][i]/1e8

#   line_thickness = 0.6
#   arrow_size = 0.6
#   matplotlib.rcParams.update({'font.size': 10})
#   matplotlib.rcParams.update({'figure.dpi': 580.0})
#   fig, axs = plt.subplots(2,6,figsize=(8.0,3.0))
#   fig.subplots_adjust(hspace=-0.2, wspace=0.15)
#   # fig.subplots_adjust(hspace=0.0, wspace=0.0)
#   index = 0
#   fig.supxlabel("x (Mm)")
#   fig.supylabel("y (Mm)",x=0.04)

#   beta = []
#   n = []
#   for i in range(len(var)):
#     mag_press = ((bx+file_vars[1][i])**2 + (by+file_vars[2][i])**2)/(8.0*np.pi)
#     beta.append(file_vars[0][i]/mag_press)
#     n.append(file_vars[3][i])
#   # plasma beta
#   max_beta = np.finfo(np.float_).min
#   min_beta = np.finfo(np.float_).max
#   max_n = np.finfo(np.float_).min
#   min_n = np.finfo(np.float_).max
#   for i in range(len(var)):
#     curr_max_beta = np.nanmax(beta[i])
#     if np.isfinite(curr_max_beta): max_beta = np.fmax(max_beta,curr_max_beta)
#     curr_min_beta = np.nanmin(beta[i])
#     if np.isfinite(curr_min_beta): min_beta = np.fmin(min_beta,curr_min_beta)
#     curr_max_n = np.nanmax(n[i])
#     if np.isfinite(curr_max_n): max_n = np.fmax(max_n,curr_max_n)
#     curr_min_n = np.nanmin(n[i])
#     if np.isfinite(curr_min_n): min_n = np.fmin(min_n,curr_min_n)
#   logrange = np.log10(max_beta) - np.log10(min_beta)
#   one_pos = (np.log10(1.0) - np.log10(min_beta))/logrange
#   small_range = min(abs(one_pos),abs(1.0-one_pos))
#   large_range = max(abs(one_pos),abs(1.0-one_pos))
#   if np.log10(max_beta)*np.log10(min_beta) < 0.0:
#     one_pos_color = (1.0,1.0,1.0)
#     zero_color = (0.0,0.1,0.4)
#     one_color = (0.4,0.1,0.0)
#     half_color_low = (0.0,0.2,0.9)
#     half_color_high = (0.9,0.2,0.0)
#     xp = [0.0,0.5*large_range]
#     if one_pos < 0.5:
#       if one_pos < 1.0/3.0:
#         s = one_pos
#         fp = [[1.0,half_color_low[0]],[1.0,half_color_low[1]],[1.0,half_color_low[2]]]
#         colors = [(0.0,[np.interp(s,xp,fp[0]),np.interp(s,xp,fp[1]),np.interp(s,xp,fp[2])]),\
#                   (one_pos,one_pos_color),(one_pos+0.5*large_range,half_color_high),(1.0,one_color)]
#       else:
#         s = one_pos-0.5*large_range
#         fp = [[half_color_low[0],zero_color[0]],[half_color_low[1],zero_color[1]],[half_color_low[2],zero_color[2]]]
#         colors = [(0.0,[np.interp(s,xp,fp[0]),np.interp(s,xp,fp[1]),np.interp(s,xp,fp[2])]),\
#                   (one_pos-0.5*large_range,half_color_low),(one_pos,one_pos_color),(one_pos+0.5*large_range,half_color_high),(1.0,one_color)]
#     elif one_pos > 0.5:
#       if one_pos > 2.0/3.0:
#         s = 1.0-one_pos
#         fp = [[1.0,half_color_high[0]],[1.0,half_color_high[1]],[1.0,half_color_high[2]]]
#         colors = [(0.0,zero_color),(one_pos-0.5*large_range,half_color_low),(one_pos,one_pos_color),\
#                   (1.0,[np.interp(s,xp,fp[0]),np.interp(s,xp,fp[1]),np.interp(s,xp,fp[2])])]
#       else:
#         s = 1.0-one_pos-0.5*large_range
#         fp = [[half_color_high[0],one_color[0]],[half_color_high[1],one_color[1]],[half_color_high[2],one_color[2]]]
#         colors = [(0.0,zero_color),(one_pos-0.5*large_range,half_color_low),(one_pos,one_pos_color),(one_pos+0.5*large_range,half_color_high),\
#                   (1.0,[np.interp(s,xp,fp[0]),np.interp(s,xp,fp[1]),np.interp(s,xp,fp[2])])]
#     else:
#       colors = [(0.0,zero_color),(one_pos-0.5*small_range,half_color_low),(one_pos,one_pos_color),(one_pos+0.5*small_range,half_color_high),(1.0,one_color)]
#   elif np.log10(max_beta) < 0.0:
#     middle = 0.5*one_pos
#     top_color = [np.interp(1.0,[middle,one_pos],[0.0,1.0]), np.interp(1.0,[middle,one_pos],[0.2,1.0]), np.interp(1.0,[middle,one_pos],[0.9,1.0])]
#     colors = [(0.0,(0.0,0.1,0.4)),(min(middle,1.0),(0.0,0.2,0.9)),(1.0,top_color)]
#   elif np.log10(max_beta) > 0.0:
#     middle = 1.0 - 0.5*(1.0-one_pos)
#     bott_color = [np.interp(0.0,[one_pos,middle],[1.0,0.9]), np.interp(0.0,[one_pos,middle],[1.0,0.2]), np.interp(0.0,[one_pos,middle],[1.0,0.0])]
#     colors = [(0.0,bott_color),(max(middle,0.0),(0.9,0.2,0.0)),(1.0,(0.4,0.1,0.0))]
#   beta_cmap = LinearSegmentedColormap.from_list('beta_diverge',colors)

#   im1 = None
#   for ax in axs[0]:
#     im = NonUniformImage(ax, animated=False, origin='lower', extent=(x_min,x_max,y_min,y_max),\
#       interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-20, base=10))
#     if index == 0: im1 = im
#     im.set_cmap(beta_cmap)
#     im.set_data(x,y,np.transpose(beta[index]))
#     im.set_clim(vmin=min_beta,vmax=max_beta)
#     ax.add_image(im)
#     ax.set_xlim(x_min, x_max)
#     ax.set_ylim(y_min, y_max)
#     ax.set_title("t="+str(round(t[index] - t[0]))+" s")
#     if index > 0:
#       ax.yaxis.set_ticklabels([])
#     # if index == 5:
#     #   fig.colorbar(im)
#     # ax.xaxis.set_ticklabels([])
#     ax.tick_params(axis='both', which='both', labelsize=6)
#     ax.label_outer()

#     if vec_var != None:
#       if vec_var == "be":
#         this_vec_x = bx
#         this_vec_y = by
#       else:
#         this_vec_x = vec_x[index].copy()
#         this_vec_y = vec_y[index].copy()
#       norm = np.sqrt(this_vec_x**2 + this_vec_y**2)
#       np.divide(this_vec_x, norm, out=this_vec_x, where=norm > 0)
#       np.divide(this_vec_y, norm, out=this_vec_y, where=norm > 0)
#       if vec_mode == "stream":
#         if os.path.exists(startpointfile) or args.streampoint_override is not None:
#           start_pts=[]
#           input_filename = ""
#           if args.streampoint_override is not None:
#             input_filename = args.streampoint_override
#           else:
#             input_filename = startpointfile
#           with open(input_filename) as f:
#             line = f.readline()
#             while line:
#               if line[0] != "#" and line != '':
#                 start_pts.append(line.split(","))
#                 start_pts[-1] = [float(s) for s in start_pts[-1]]
#               line = f.readline()
#         else:
#           if args.stream_number is None:
#             num_points = 11
#           else:
#             num_points = args.stream_number
#           if args.uniform_stream:
#             stream_points = np.linspace(x[1],x[-2],num_points)
#           else:
#             if args.stream_power is None:
#               if this_vec_y[1,1]*this_vec_y[-2,1] < 0.0:
#                 print("loop detected")
#                 num_points = int(num_points*2)
#                 stream_pow = 1
#               else:
#                 stream_pow = 2
#             else:
#               stream_pow = args.stream_power
#             norm_normalized = norm[:,0]**stream_pow/np.trapz(norm[:,0]**stream_pow)
#             cdf = [np.trapz(norm_normalized[:(i+1)]) for i in range(len(norm_normalized))]
#             stream_points = np.interp(np.linspace(0.05,0.95,num=num_points),cdf,x)
#           start_pts=np.column_stack((stream_points,y_eq[0]*np.ones_like(stream_points)))
#         start_pts = np.asarray(start_pts)
#         stream = ax.streamplot(x_eq,y_eq,this_vec_x.transpose(),this_vec_y.transpose(),\
#           start_points=start_pts,\
#             color=(0.0,0.0,0.0),broken_streamlines=False,linewidth=line_thickness,arrowstyle='->',arrowsize=arrow_size,maxlength=10.0)
#         stream.lines.set(pickradius=2.0)
#     ax.set_aspect('equal')

#     index += 1

#   index = 0
#   im2 = None
#   # density
#   for ax in axs[1]:
#     im = NonUniformImage(ax, animated=False, origin='lower', extent=(x_min,x_max,y_min,y_max),\
#       interpolation='nearest', norm=matplotlib.colors.SymLogNorm(linthresh=1e-20, base=10))
#     if index == 0: im2 = im
#     im.set_data(x,y,np.transpose(n[index]))
#     im.set_clim(vmin=min_n,vmax=max_n)
#     ax.add_image(im)
#     ax.set_xlim(x_min, x_max)
#     ax.set_ylim(y_min, y_max)
#     # ax.set_title(str(index))
#     if index > 0:
#       ax.yaxis.set_ticklabels([])
#     ax.tick_params(axis='both', which='both', labelsize=6)
#     ax.label_outer()

#     if vec_var != None:
#       if vec_var == "be":
#         this_vec_x = bx
#         this_vec_y = by
#       else:
#         this_vec_x = vec_x[index].copy()
#         this_vec_y = vec_y[index].copy()
#       norm = np.sqrt(this_vec_x**2 + this_vec_y**2)
#       np.divide(this_vec_x, norm, out=this_vec_x, where=norm > 0)
#       np.divide(this_vec_y, norm, out=this_vec_y, where=norm > 0)
#       if vec_mode == "stream":
#         if os.path.exists(startpointfile) or args.streampoint_override is not None:
#           start_pts=[]
#           input_filename = ""
#           if args.streampoint_override is not None:
#             input_filename = args.streampoint_override
#           else:
#             input_filename = startpointfile
#           with open(input_filename) as f:
#             line = f.readline()
#             while line:
#               if line[0] != "#" and line != '':
#                 start_pts.append(line.split(","))
#                 start_pts[-1] = [float(s) for s in start_pts[-1]]
#               line = f.readline()
#         else:
#           if args.stream_number is None:
#             num_points = 11
#           else:
#             num_points = args.stream_number
#           if args.uniform_stream:
#             stream_points = np.linspace(x[1],x[-2],num_points)
#           else:
#             if args.stream_power is None:
#               if this_vec_y[1,1]*this_vec_y[-2,1] < 0.0:
#                 print("loop detected")
#                 num_points = int(num_points*2)
#                 stream_pow = 1
#               else:
#                 stream_pow = 2
#             else:
#               stream_pow = args.stream_power
#             norm_normalized = norm[:,0]**stream_pow/np.trapz(norm[:,0]**stream_pow)
#             cdf = [np.trapz(norm_normalized[:(i+1)]) for i in range(len(norm_normalized))]
#             stream_points = np.interp(np.linspace(0.05,0.95,num=num_points),cdf,x)
#           start_pts=np.column_stack((stream_points,y_eq[0]*np.ones_like(stream_points)))
#         start_pts = np.asarray(start_pts)
#         stream = ax.streamplot(x_eq,y_eq,this_vec_x.transpose(),this_vec_y.transpose(),\
#           start_points=start_pts,\
#             color=(0.0,0.0,0.0),broken_streamlines=False,linewidth=line_thickness,arrowstyle='->',arrowsize=arrow_size,maxlength=10.0)
#         stream.lines.set(pickradius=2.0)
#     index += 1
#     ax.set_aspect('equal')
#   # plt.tight_layout()
#   fig.subplots_adjust(left=0.1,right=0.895)
#   cax1 = plt.axes([0.91, 0.514, 0.01, 0.305])
#   cbar1 = fig.colorbar(im1,cax=cax1,label="Plasma Beta")
#   for label in cbar1.ax.yaxis.get_ticklabels()[::2]:
#     label.set_visible(False)
#   cax2 = plt.axes([0.91, 0.172, 0.01, 0.305])
#   cbar2 = fig.colorbar(im2,cax=cax2,label="Num. Density")
#   plt.savefig("testing_200heat.png")
#   exit()
# ###############################################

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

if args.low_colorbar_vec is not None:
  min_v_vec = args.low_colorbar_vec
if args.high_colorbar_vec is not None:
  max_v_vec = args.high_colorbar_vec


frames = [0]
if args.still_frames != "": frames = [int(k) for k in args.still_frames.split(",")]
num_frames = len(frames)
i_frame = 0

# wr = [1.0 for f in frames]
# if num_frames > 1:
#   wr.append(0.3)

fig, axis = plt.subplots(nrows=1,ncols=num_frames,sharey=True,figsize=(args.fig_size_x, args.fig_size_y)) #6.4, 4.8
# if num_frames > 1:
#   axs = axis
# else:
#   axs = [axis]

for ax in fig.axes:

  def construct_plot_title(f):
    global args
    global t
    plot_title = ""
    if args.diff_file is not None:
      plot_title += "Difference in "
    if not args.no_varname:
      if args.full_varname: plot_title += fullnames[output_var]
      else: plot_title += output_var
    if not args.no_units:
      if args.diff_file is not None:
        plot_title += " (Fraction)"
      else:
        plot_title += fullunits.get(output_var,"")
    if args.time_label:
      if len(plot_title) != 0: plot_title += ", "
      if args.time_label_rounded:
        if args.time_label_zero:
          plot_title += "t="+str(round(t[f]-t[0]))+" s"
        else:
          plot_title += "t="+str(round(t[f]))+" s"
      else:
        if args.time_label_zero:
          plot_title += "t="+str(t[f]-t[0])+" s"
        else:
          plot_title += "t="+str(t[f])+" s"
    return plot_title

  if args.no_axes: ax.axis('off')

  im = NonUniformImage(ax, origin='lower')
  im.set(animated=True, extent=(x_min,x_max,y_min,y_max), interpolation=('bilinear' if aia.is_filter(output_var) or output_var == "aia_composite" else 'nearest'))

  # Set colorbar scaling, with any other particular settings
  for k in symlogthresholds.keys():
    if output_var in k:
      thresh = (symlogthresholds[k] if args.symlogthreshold == -1.0 else args.symlogthreshold)
      im.set(norm=matplotlib.colors.SymLogNorm(linthresh=thresh, base=10))
      if output_var in ["cumulative_electron_heating","cumulative_ion_heating","cumulative_electron_heating_temp","cumulative_ion_heating_temp","div_b"]:
        im.set(cmap='seismic')
        biggest_mag = max(abs(min_v),abs(max_v))
        im.set_clim(vmin=-biggest_mag,vmax=biggest_mag)
      break
  else:
    if args.diff_file is not None:
      thresh = (1e-15 if args.symlogthreshold == -1.0 else args.symlogthreshold)
      im.set(norm=matplotlib.colors.SymLogNorm(linthresh=thresh, base=10))
    elif output_var in ["beta"]:
      im.set(norm=matplotlib.colors.CenteredNorm(),cmap='seismic')
    elif aia.is_filter(output_var):
      # im.set(norm=matplotlib.colors.PowerNorm(gamma=0.4),cmap=aia.colormap(output_var))
      im.set(norm=ImageNormalize(stretch=source_stretch(OrderedDict(), AsinhStretch(0.001)), clip=False),cmap=ct.aia_color_table(int((output_var).replace("A",""))*u.angstrom))
      # for f in range(len(var)):
      #   var[f] = np.square(var[f])
    else:
      thresh = (1e-2 if args.symlogthreshold == -1.0 else args.symlogthreshold)
      im.set(norm=matplotlib.colors.SymLogNorm(linthresh=thresh, base=10))

  if output_var not in ["cumulative_electron_heating","cumulative_ion_heating","cumulative_electron_heating_temp","cumulative_ion_heating_temp","div_b"]:
    im.set_clim(vmin=min_v,vmax=max_v)

  if output_var in ["anomalous_template"]:
    im.set(cmap='Reds',norm=matplotlib.colors.Normalize(vmin=0,vmax=1))
    

  frame = frames[i_frame]
  if output_var == "aia_composite":
    im.set_data(x,y,var[frame].transpose(1,0,2))
  else:
    im.set_data(x,y,np.transpose(var[frame]))
  ax.add_image(im)
  if args.contour_list != "":
    cs = ax.contour(x,y,np.transpose(var[frame]), colors='k',levels=contour_vals)
    csl = ax.clabel(cs,inline=True)
  if args.equal_aspect:
    ax.set_aspect('equal')

  ax.set(xlim=(x_min,x_max), ylim=(y_min,y_max))
  if args.internal_title:
    if output_var == "aia_composite":
      ann = [ax.annotate(output_var,(0.05,0.95),xycoords='axes fraction',ha='left',va='top'),
            ax.annotate(output_var,(0.05,0.9),xycoords='axes fraction',ha='left',va='top'),
            ax.annotate(output_var,(0.05,0.85),xycoords='axes fraction',ha='left',va='top')]
    else:
      ann = ax.annotate(output_var+", t="+str(t[frame]),(0.05,0.95),xycoords='axes fraction',ha='left',va='top')
  else:
    ax.set(title=construct_plot_title(frame))
  
  if args.small_ticklabels:
    ax.tick_params(axis='both', which='major', labelsize=0.6*args.font_size)
    ax.tick_params(axis='both', which='minor', labelsize=0.4*args.font_size)
  if args.no_ticks_x or args.no_ticks:
    ax.xaxis.set_ticklabels([])
  else:
    ax.set(xlabel="x (Mm)")
  if args.no_ticks_y or args.no_ticks:
    ax.yaxis.set_ticklabels([])
  elif i_frame == 0:
    ax.set(ylabel="y (Mm)")
  if not args.no_colorbar and output_var != "aia_composite" and not (num_frames > 1 and i_frame == num_frames-1):
    if num_frames == 1:
      var_colorbar = fig.colorbar(im,location=args.colorbar_location,label=fullnames[output_var]+" "+fullunits.get(output_var,""))
      if args.colorbar_label: var_colorbar.set_label(shortnames[output_var]+ fullunits.get(output_var,""))
    if args.tick_list != "":
      tick_vals = [float(v) for v in args.tick_list.split(',')]
      if args.tick_e:
        var_colorbar.set_ticks(tick_vals,\
          labels=[f'{v:.1e}'.replace("e-0","e-").replace("e+0","e") for v in tick_vals])
      else:
        var_colorbar.set_ticks(tick_vals)

  if args.tracerparticle_timeseries or args.crosshair_timeseries:
    tp_measurements = []
    tp_measure_interp = RegularGridInterpolator((t,x.copy(order='C'),y.copy(order='C')),np.array(var),fill_value=0.0,bounds_error=False)
    if args.crosshair_timeseries:
      chx = float(args.crosshair_timeseries.split("_")[0].split("x")[1])
      chy = float(args.crosshair_timeseries.split("_")[1].split("y")[1])
    for i in range(len(t)):
      try:
        if args.crosshair_timeseries:
          tp_measurements.append(tp_measure_interp([t[i],chx,chy])[0])
        else:
          tpm_idx = tracer_labels[i].index(args.tracerparticle_timeseries)
          tp_measurements.append(tp_measure_interp([t[i],tracer_particles[i][0][tpm_idx],tracer_particles[i][1][tpm_idx]])[0])
        # print(f"{tp_measurements[-1]},",end="")
      except ValueError:
        tp_measurements.append(np.NaN)
    del tp_measure_interp
    figts, axts = plt.subplots(figsize=(args.fig_size_x, args.fig_size_y))
    axts.plot(t,np.array(tp_measurements))
    # axts.plot(t,np.array(tp_measurements)/1.0e6 - from_ms_elec[:len(tp_measurements)],label='elec')
    # axts.plot(t,np.array(tp_measurements)/1.0e6 - from_ms_ion[:len(tp_measurements)],label='ion')
    axts.set_yscale('log')
    axts.set_facecolor('white')
    axts.set_xlabel("Time (s)")
    axts.set_ylabel(f"{fullnames[output_var]} {fullunits.get(output_var,'')}")
    # axts.set_ylabel("Density. (g cm^-3)")
    # axts.legend()
    # axts.set_ylim([0.04,15])
    # axts.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axts.set_yticks([0.1e6,1.0e6,10.0e6])
    plt.tight_layout()
    if args.crosshair_timeseries:
      figts.savefig(args.filename.split('/')[-2]+"_"+output_var+"_"+args.crosshair_timeseries+"_timeseries.png")
    else:
      figts.savefig(args.filename.split('/')[-2]+"_"+output_var+"_"+args.tracerparticle_timeseries+"_timeseries.png")
    exit()

  if args.tracerparticles: 
    tp = ax.scatter(np.array(tracer_particles[frame][0]),np.array(tracer_particles[frame][1]),marker=".",c=('k' if output_var == "beta" else 'w'),s=args.tracerparticles_size,zorder=5,edgecolor='w',linewidth=0.025*args.tracerparticles_size)
    if args.label_tracerparticles:
      tp_labels = []
      for i in range(len(tracer_labels[frame])):
        tp_labels.append(ax.annotate(tracer_labels[frame][i],(tracer_particles[frame][0][i],tracer_particles[frame][1][i]),ha='center',va='bottom',c='w',fontsize='x-small'))
    if args.tracerparticle_measure:
      tp_measurements = []
      tp_measure_interp = RegularGridInterpolator((t,x.copy(order='C'),y.copy(order='C')),np.array(var),fill_value=0.0,bounds_error=False)
      for i in range(len(t)):
        try:
          tpm_idx = tracer_labels[i].index(args.tracerparticle_measure)
          tp_measurements.append(f"{tp_measure_interp([t[i],tracer_particles[i][0][tpm_idx],tracer_particles[i][1][tpm_idx]])[0]:.3e}"
                                          +fullunits.get(output_var,"").replace("(","").replace(")","")+f" at ({tracer_particles[i][0][tpm_idx]:.1f},{tracer_particles[i][1][tpm_idx]:.1f})")
        except ValueError:
          tp_measurements.append("null")
      del tp_measure_interp
      tp_measure = ax.annotate(f"{args.tracerparticle_measure}: {tp_measurements[frame]}",(0.01,0.01),xycoords='figure fraction',ha='left',va='bottom')


  if args.crosshair != "":
    ch_strs = args.crosshair.split(",")
    chs_x = [float(c.split("_")[0][1:]) for c in ch_strs]
    chs_y = [float(c.split("_")[1][1:]) for c in ch_strs]
    ch = ax.scatter(chs_x,chs_y,c='k',edgecolors='w',zorder=10,marker='*',linewidth=0.005*args.crosshair_size,s=args.crosshair_size)
    if args.crosshair_labels != "":
      ch_ls = args.crosshair_labels.split(",")
      if len(ch_ls) != len(chs_x):
        print("Crosshair labels should be one per crosshair")
        exit()
      for i in range(len(ch_ls)):
        x_offset = 12
        y_offset = 0
        y_offset_increment = 8
        if i < len(ch_ls)-1 and chs_x[i+1] > chs_x[i]:
          if chs_y[i+1] < chs_y[i]: y_offset += y_offset_increment
          elif chs_y[i+1] > chs_y[i]: y_offset -= y_offset_increment
        if i > 0 and chs_x[i-1] > chs_x[i]:
          if chs_y[i-1] < chs_y[i]: y_offset += y_offset_increment
          elif chs_y[i-1] > chs_y[i]: y_offset -= y_offset_increment
        ax.annotate(ch_ls[i],(chs_x[i]+x_offset,chs_y[i]+y_offset),va='center',ha='left',c=('k' if output_var == "beta" else 'w'),size='x-small')

  if args.tracerparticle_record:
    tp = ax.scatter([],[],c='k',s=3)
    tp_labels = []
    for i in range(len(tracer_labels)):
      tp_labels.append(ax.annotate(tracer_labels[0][i],(tracer_particles[frame][0][i],tracer_particles[frame][1][i]),ha='center',va='bottom'))

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

    if vec_mode == "stream":
      interp_x = LinearNDInterpolator(list(zip(X.flatten(), Y.flatten())), this_vec_x.flatten(),fill_value=0.0)
      interp_y = LinearNDInterpolator(list(zip(X.flatten(), Y.flatten())), this_vec_y.flatten(),fill_value=0.0)
      x_eq_grid, y_eq_grid = np.meshgrid(x_eq,y_eq)
      this_vec_x = interp_x(x_eq_grid, y_eq_grid)
      this_vec_y = interp_y(x_eq_grid, y_eq_grid)
      
      norm = np.sqrt(this_vec_x**2 + this_vec_y**2)
      np.divide(this_vec_x, norm, out=this_vec_x, where=norm > 0)
      np.divide(this_vec_y, norm, out=this_vec_y, where=norm > 0)
      if args.tracer_particle_streampoint:
        start_pts_series = extract_tracer_particles_from_file(tracerparticlefile, display_interval, start_time, end_time, True, [])
        start_pts_series = [np.transpose(np.asarray(s)) for s in start_pts_series[0]]
      elif os.path.exists(startpointfile) or args.streampoint_override is not None:
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
        start_pts=[[x[0],y[0]]]
        # if args.stream_number is None:
        #   num_points = 11
        # else:
        #   num_points = args.stream_number
        # if args.uniform_stream:
        #   stream_points = np.linspace(x[1],x[-2],num_points)
        # else:
        #   if args.stream_power is None:
        #     if this_vec_y[1,1]*this_vec_y[-2,1] < 0.0:
        #       print("loop detected")
        #       num_points = int(num_points*2)
        #       stream_pow = 1
        #     else:
        #       stream_pow = 2
        #   else:
        #     stream_pow = args.stream_power
        #   norm_normalized = norm[:,0]**stream_pow/np.trapz(norm[:,0]**stream_pow)
        #   cdf = [np.trapz(norm_normalized[:(i+1)]) for i in range(len(norm_normalized))]
        #   stream_points = np.interp(np.linspace(0.05,0.95,num=num_points),cdf,x)
        # start_pts=np.column_stack((stream_points,y_eq[0]*np.ones_like(stream_points)))
      if args.tracer_particle_streampoint:
        start_pts = start_pts_series[frame]
        cleaned = []
        for s in start_pts:
          if s[0] > x[1] and s[0] < x[-2] and s[1] > y[1] and s[1] < y[-2]:
            cleaned.append(s)
        start_pts = cleaned
      start_pts = np.asarray(start_pts)
      # print(start_pts)
      stream = ax.streamplot(x_eq,y_eq,this_vec_x,this_vec_y,start_points=start_pts,\
          color=(0.0,0.0,0.0),broken_streamlines=False,linewidth=line_thickness,arrowstyle='->',arrowsize=arrow_size,maxlength=10.0)
      stream.lines.set(pickradius=2.0)
    else:
      norm = np.sqrt(this_vec_x**2 + this_vec_y**2)
      np.divide(this_vec_x, norm, out=this_vec_x, where=norm > 0)
      np.divide(this_vec_y, norm, out=this_vec_y, where=norm > 0)
      assert vec_mode == "quiver"
      quiv = ax.quiver(X[x_vec_indices, y_vec_indices], Y[x_vec_indices, y_vec_indices], \
          this_vec_x[x_vec_indices, y_vec_indices], this_vec_y[x_vec_indices, y_vec_indices], norm[x_vec_indices, y_vec_indices], \
                cmap=plt.cm.Greys, ec='k', lw=0.0, scale_units='inches', angles='xy', scale=10, \
                  width = 0.005, headwidth=3, headlength=3, headaxislength=2.5, \
                    pivot='mid', norm=matplotlib.colors.SymLogNorm(linthresh=1e-4, base=10))
      quiv.set_clim(vmin=min_v_vec,vmax=max_v_vec)
      vec_colorbar = fig.colorbar(quiv)
      vec_colorbar.set_label(label=shortnames[output_var]+ fullunits.get(output_var,""))
      vector_color_axes = fig.axes[-1]

  plt.tight_layout()
  if args.no_margin:
    fig.subplots_adjust(left=(0.0 if args.margin_x_lower < 0.0 else args.margin_x_lower),
                        right=(1.0 if args.margin_x_upper < 0.0 else (1.0 - args.margin_x_upper)),
                        bottom=(0.0 if args.margin_y_lower < 0.0 else args.margin_y_lower),
                        top=(1.0 if args.margin_y_upper < 0.0 else (1.0 - args.margin_y_upper)))
  else:
    if args.margin_x_lower >= 0.0 : fig.subplots_adjust(left=args.margin_x_lower)
    if args.margin_x_upper >= 0.0 : fig.subplots_adjust(right=(1.0 - args.margin_x_upper))
    if args.margin_y_lower >= 0.0 : fig.subplots_adjust(bottom=args.margin_y_lower)
    if args.margin_y_upper >= 0.0 : fig.subplots_adjust(top=(1.0 - args.margin_y_upper))

  i_frame += 1
  if num_frames > 1:
    # if i_frame != 1:
    #   ax.set(ylabel="")
    if i_frame == num_frames:
      fig.subplots_adjust(wspace=0.05,right=0.85,bottom=0.25,top=0.8)
      cbar_ax = fig.add_axes([0.86, 0.29, 0.02, 0.47])
      var_colorbar = fig.colorbar(im,cax=cbar_ax,location=args.colorbar_location)
      if args.colorbar_label: var_colorbar.set_label(label=shortnames[output_var]+" "+fullunits.get(output_var,""))
      fig.savefig(args.filename.split('/')[-2]+"_"+output_var+"_stills.png")
      exit()
    else:
      # var_colorbar.remove()
      if vec_var != None and vec_mode == "quiver":
        vec_colorbar.remove()

  
  frame = -2

  def updatefig(*args_arg):
      global frame
      global ax
      global quiv
      global stream
      global start_pts
      global im
      global tp, tp_labels, tp_measure, tp_measurements
      global tracer_particles, tracer_labels
      global cs, csl
      global ann
      if not (args.streampoint_record or args.tracerparticle_record or args.interactive):
        frame = (frame + 1)%len(var)
      if output_var == "aia_composite":
        im.set_data(x,y,var[frame].transpose(1,0,2))
      else:
        im.set_data(x,y,np.transpose(var[frame]))
      if args.contour_list != "":
        for c in cs.collections:
          c.remove()
        for c in csl:
          c.remove()
        cs = ax.contour(x,y,np.transpose(var[frame]), colors='k',levels=contour_vals)
        csl = ax.clabel(cs,inline=True)
      if args.tracerparticles:
        tp.set_offsets(np.transpose(np.array(tracer_particles[frame])))
        if args.label_tracerparticles:
          for i in range(len(tp_labels)):
            tp_labels.pop().remove()
          for i in range(len(tracer_labels[frame])):
            tp_labels.append(ax.annotate(tracer_labels[frame][i],(tracer_particles[frame][0][i],tracer_particles[frame][1][i]),ha='center',va='bottom',c='w',fontsize='xx-small'))
        if args.tracerparticle_measure:
          tp_measure.remove()
          tp_measure = ax.annotate(f"{args.tracerparticle_measure}: {tp_measurements[frame]}",(0.01,0.01),xycoords='figure fraction',ha='left',va='bottom')
      if args.tracerparticle_record:
        tp.set_offsets(np.transpose(np.array(tracer_particles)))
        for i in range(len(tp_labels)):
          tp_labels.pop().remove()
        for i in range(len(tracer_labels)):
          tp_labels.append(ax.annotate(tracer_labels[i],(tracer_particles[0][i],tracer_particles[1][i]),ha='center',va='bottom'))
      plot_title = construct_plot_title(frame)
      if args.internal_title:
        if output_var == "aia_composite":
          if ann is not None:
            for a in ann: a.remove()
          filters = fullnames[output_var].split(",")
          ann = [ax.annotate(filters[0],(0.05,0.95),c='r',xycoords='axes fraction',ha='left',va='top'),
                ax.annotate(filters[1],(0.05,0.9),c='g',xycoords='axes fraction',ha='left',va='top'),
                ax.annotate(filters[2],(0.05,0.85),c='b',xycoords='axes fraction',ha='left',va='top')]
        else:
          if ann is not None: ann.remove()
          ann = ax.annotate(plot_title,(0.05,0.95),xycoords='axes fraction',ha='left',va='top')
      else:
        ax.set(title=plot_title)
      if not (args.no_ticks or args.no_ticks_x):
        ax.set(xlabel="x (Mm)")
      if not (args.no_ticks or args.no_ticks_y or i_frame > 0):
        ax.set(ylabel="y (Mm)")
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
          interp_x = LinearNDInterpolator(list(zip(X.flatten(), Y.flatten())), this_vec_x.flatten(),fill_value=0.0)
          interp_y = LinearNDInterpolator(list(zip(X.flatten(), Y.flatten())), this_vec_y.flatten(),fill_value=0.0)
          x_eq_grid, y_eq_grid = np.meshgrid(x_eq,y_eq)
          this_vec_x = interp_x(x_eq_grid, y_eq_grid)
          this_vec_y = interp_y(x_eq_grid, y_eq_grid)
          stream.lines.remove()
          for art in ax.get_children():
            if not isinstance(art, matplotlib.patches.FancyArrowPatch):
              continue
            art.remove()
          if args.tracer_particle_streampoint:
            start_pts = start_pts_series[frame]
            cleaned = []
            for s in start_pts:
              if s[0] > x[1] and s[0] < x[-2] and s[1] > y[1] and s[1] < y[-2]:
                cleaned.append(s)
            start_pts = cleaned
          stream = ax.streamplot(x_eq,y_eq,this_vec_x,this_vec_y,start_points=start_pts,\
            color=(0.0,0.0,0.0),broken_streamlines=False,linewidth=line_thickness,arrowstyle='->',arrowsize=arrow_size,maxlength=10.0)
          if not (args.streampoint_record or args.tracerparticle_record):
            print('\x1b[1K\r'+"frame " + str(frame) + " plotted", end='')
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
      while i<len(segments) and np.sqrt(((segments[i][0,0]-segments[i-1][1,0])**2+(segments[i][0,1]-segments[i-1][1,1])**2))<0.01:
        this_line.append(segments[i][0])
        this_line_segments.append(segments[i])
        i+= 1
      lines.append(np.asarray(this_line))
      line_segments.append(this_line_segments)
    return lines, line_segments

  def modifytracerpoints(event):
    ix, iy = event.xdata, event.ydata
    print(f"{ix},{iy} clicked")
    if ix is not None and iy is not None and args.tracerparticle_record:
      global tp, tp_labels
      global fig
      global ax
      global tracer_particles, tracer_labels
      global click_count
      tracer_clicked = False
      clicked_i = -1
      for i in range(len(tracer_labels)):
        if np.sqrt((tracer_particles[0][i] - ix)**2 + (tracer_particles[1][i] - iy)**2) < 3:
          print(f"{tracer_particles[0][i]},{tracer_particles[1][i]} clicked (particle {i}, {tracer_labels[i]})")
          tracer_clicked = True
          clicked_i = i
          break
      if tracer_clicked:
        print("Removing existing tracer...")
        tracer_particles[0].pop(clicked_i)
        tracer_particles[1].pop(clicked_i)
        tracer_labels.pop(clicked_i)
      else:
        print("Adding new tracer...")
        tracer_particles[0].append(ix)
        tracer_particles[1].append(iy)
        click_count += 1
        tracer_labels.append("t"+str(click_count))
    updatefig()
    fig.canvas.draw_idle()

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
          single_line = matplotlib.collections.LineCollection(line_segments[l],pickradius=2.0,transform=ax.transData)
          if single_line.contains(event)[0]:
            print("Deleting existing line")
            segs = single_line.get_segments()
            for j in range(len(segs)):
              segs[j] = segs[j].tolist()
            segs_flat = [sublist[0] for sublist in segs]
            curr_min_dist = np.finfo(np.float_).max
            curr_best_k = -1
            for k in range(len(start_pts)):
              curr_pt = start_pts[k]
              closest_dist = min([np.sqrt((curr_pt[0]-line_pt[0])**2 + (curr_pt[1]-line_pt[1])**2) for line_pt in segs_flat])
              if closest_dist < curr_min_dist:
                curr_min_dist = closest_dist
                curr_best_k = k
            start_pts = np.delete(start_pts,curr_best_k,axis=0)
      else:
        print("Adding new line")
        start_pts = np.append(start_pts,[[ix, iy]],axis=0)
    updatefig()
    fig.canvas.draw_idle()

  def savetracerpoints(event):
    global stream
    global args
    global tracerrecordfile
    global tracer_particles, tracer_labels
    print(f"Saving {len(tracer_labels)} tracer points to {tracerrecordfile}...")
    scale_factor = 1.0e8
    with open(tracerrecordfile, 'w') as f:
      for i in range(len(tracer_labels)):
        f.write(f"{tracer_particles[0][i]*scale_factor:.15e},{tracer_particles[1][i]*scale_factor:.15e}#{tracer_labels[i]}\n")

  def savestreampoints(event):
    global stream
    global args
    global startpointfile
    y_min = Y[0][50]
    segments = stream.lines.get_segments()
    lines = extract_separate_lines(segments)[0]
    print(f"Saving {len(lines)} lines to {startpointfile}...")
    scale_factor = 1.0e8
    with open(startpointfile, 'w') as f:
      line_num = 1
      for l in lines:
        if l[-1][1] > l[0][1]:
          if l[0][1] < y_min:
            for i in range(1,len(l)):
              if l[i][1] > y_min:
                print(i)
                f.write(f"{l[i][0]*scale_factor:.15e},{l[i][1]*scale_factor:.15e}#s{line_num}\n")
                break
            else:
              print("Streamline always below y_min. Using middle point...")
              middle = int(round(0.5*len(l)))
              f.write(f"{l[middle][0]*scale_factor:.15e},{l[middle][1]*scale_factor:.15e}#s{line_num}\n")
          else:
            f.write(f"{l[0][0]*scale_factor:.15e},{l[0][1]*scale_factor:.15e}#s{line_num}\n")
        else:
          if l[-1][1] < y_min:
            for i in range(2,len(l)):
              if l[-i][1] > y_min:
                print(i)
                f.write(f"{l[-i][0]*scale_factor:.15e},{l[-i][1]*scale_factor:.15e}#s{line_num}\n")
                break
            else:
              print("Streamline always below y_min. Using middle point...")
              middle = int(round(0.5*len(l)))
              f.write(f"{l[middle][0]*scale_factor:.15e},{l[middle][1]*scale_factor:.15e}#s{line_num}\n")
          else:
            f.write(f"{l[-1][0]*scale_factor:.15e},{l[-1][1]*scale_factor:.15e}#s{line_num}\n")
        line_num += 1

  def setmeasurepoint(event):
    global measurepoint
    global mp_scatter
    global ax
    global x,y
    if event.xdata is None or event.ydata is None:
      measurepoint = None
    else:
      xidx = min(np.searchsorted(x,event.xdata),xdim-1)
      if xidx != 0:
        xidx = xidx if (event.xdata-x[xidx-1])>(x[xidx]-event.xdata) else xidx-1
      yidx = min(np.searchsorted(y,event.ydata),ydim-1)
      if yidx != 0:
        yidx = yidx if (event.ydata-y[yidx-1])>(y[yidx]-event.ydata) else yidx-1
      measurepoint = (xidx,yidx)
    try:
      mp_scatter
    except NameError:
        real_coords = (x[measurepoint[0]],y[measurepoint[1]])
        mp_scatter = ax.scatter((real_coords[0],),(real_coords[1],),c='r',s=10)
    else:
        if measurepoint is None:
          mp_scatter.set_offsets(np.ma.masked_array([0, 0], mask=True))
        else:
          mp_scatter.set_offsets(((x[measurepoint[0]],y[measurepoint[1]]),))
    # fig.canvas.draw_idle()

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
      if(event.key == "0"):
        savestreampoints(event)
      elif(event.key == "9"):
        drawautofieldlines()
    if args.tracerparticle_record:
      if(event.key == "0"):
        savetracerpoints(event)

  def displayinfo(event):
    global fig
    global x,y
    global var
    global frame
    global measurepoint
    global mp_plot
    if event.xdata is None or event.ydata is None:
      fig.suptitle("")
      try:
        mp_plot
      except NameError:
        pass
      else:
        if measurepoint is None:
          mp_plot.set_data([[],[]])
    else:
      xidx = min(np.searchsorted(x,event.xdata),xdim-1)
      if xidx != 0:
        xidx = xidx if (event.xdata-x[xidx-1])>(x[xidx]-event.xdata) else xidx-1
      yidx = min(np.searchsorted(y,event.ydata),ydim-1)
      if yidx != 0:
        yidx = yidx if (event.ydata-y[yidx-1])>(y[yidx]-event.ydata) else yidx-1
      suptext = f"{var[frame][xidx,yidx]:.3e} at ({xidx},{yidx})"
      if measurepoint is not None:
        measuredist = np.sqrt((measurepoint[0] - xidx)**2 + (measurepoint[1] - yidx)**2)
        suptext += f"\n{measuredist:.2f} cells from ({measurepoint[0]},{measurepoint[1]})"
        try:
          mp_plot
        except NameError:
          mp_plot = ax.plot((x[measurepoint[0]],event.xdata),(y[measurepoint[1]],event.ydata),'r--')[0]
        else:
          mp_plot.set_data([[x[measurepoint[0]],event.xdata],[y[measurepoint[1]],event.ydata]])
      fig.suptitle(suptext,x=0.99,y=0.99,ha='right',va='top',size='small')
    fig.canvas.draw_idle()

  def drawautofieldlines():
    global start_pts
    global stream
    global fig
    global ax
    global x_min, x_max
    global y_min, y_max
    num_cells_x = 3
    num_cells_y = 6
    dx = (x_max - x_min)/num_cells_x
    dy = (y_max - y_min)/num_cells_y
    curr_x = x_min + 0.5*dx
    curr_y = y_max - 0.5*dy
    while curr_y > y_min:
      curr_x = x_min + 0.5*dx
      while curr_x < x_max:
        min_d = (ax.transData.transform((x_min + min(dx,dy), 0))[0] - ax.transData.transform((x_min, 0))[0])*0.5
        pix_pos = ax.transData.transform((curr_x, curr_y))
        faux_me = matplotlib.backend_bases.MouseEvent('button_press_event', fig.canvas, pix_pos[0], pix_pos[1])
        segments = stream.lines.get_segments()
        lines, line_segments = extract_separate_lines(segments)
        covered = False
        for l in range(len(lines)):
          single_line = matplotlib.collections.LineCollection(line_segments[l],pickradius=0.5*min_d,transform=ax.transData)
          if single_line.contains(faux_me)[0]:
            covered = True
            break
        if not covered:
          start_pts = np.append(start_pts,[[curr_x, curr_y]],axis=0)
          updatefig()
          fig.canvas.draw_idle()
        curr_x += dx
      curr_y -= dy
    
  if num_frames == 1:
    if args.interactive or args.streampoint_record or args.tracerparticle_record:
      cid = fig.canvas.mpl_connect('key_press_event', keypresshandler)
    if args.interactive:
      cid = fig.canvas.mpl_connect('motion_notify_event', displayinfo)
      cid2 = fig.canvas.mpl_connect('button_press_event', setmeasurepoint)
    if args.streampoint_record:
      cid = fig.canvas.mpl_connect('button_press_event', modifystreampoints)
    if args.tracerparticle_record:
      cid = fig.canvas.mpl_connect('button_press_event', modifytracerpoints)
    else:
      ani = animation.FuncAnimation(fig, updatefig, frames=len(var), repeat=args.realtime, interval=100, blit=False)

    if args.realtime or args.streampoint_record or args.tracerparticle_record or args.interactive:
      measurepoint = None
      plt.show()
    else:
      FFwriter = animation.FFMpegWriter(bitrate=2000*2*4*10,fps=20,metadata={'comment': " ".join(sys.argv[:])})
      if args.out_filename != "":
        fname = args.out_filename + '.mp4'
      else:
        filename_suffix = "_"+output_var
        if vec_var != None: filename_suffix = filename_suffix + "_" + vec_var
        if args.sandbox: filename_suffix = filename_suffix + "_SANDBOX"
        fname = args.filename.split('/')[-2]+filename_suffix+'.mp4'
      ani.save(fname, writer = FFwriter, dpi=200)
      print()