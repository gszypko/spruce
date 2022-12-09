import numpy as np
import math
import sys
import aia_response as aia

def read_grid(in_file,xdim,ydim,xl,yl,xu,yu):
  result = []
  for i in range(xdim):
    row = in_file.readline().split(',')[yl:yu]
    if i >= xl and i < xu:
      result.append(row)
  return np.asarray(result).astype(np.float64)

def div(a_x,a_y,x,y):
  result = np.zeros_like(a_x)
  for i in range(1,len(x)-1):
    for j in range(1,len(y)-1):
      result[i,j] = (a_x[i+1,j] - a_x[i-1,j])/(x[i+1]-x[i-1]) + (a_y[i,j+1] - a_y[i,j-1])/(y[j+1]-y[j-1])
  return result

def apply_contour_computation(output_var, file_vars, x, y, bx, by, bz):
  var = []
  if output_var == "beta":
    for i in range(len(file_vars[0])):
        # mag_press = ((bx+file_vars[1][i])**2 + (by+file_vars[2][i])**2 + (bz+file_vars[3][i])**2)/(8.0*np.pi)
        mag_press = ((bx+file_vars[1][i])**2 + (by+file_vars[2][i])**2)/(8.0*np.pi)
        var.append(file_vars[0][i]/mag_press)
  elif output_var == "b_mag":
    for i in range(len(file_vars[0])):
        # field_mag = np.sqrt((bx+file_vars[0][i])**2 + (by+file_vars[1][i])**2 + (bz+file_vars[2][i])**2)
        field_mag = np.sqrt((bx+file_vars[0][i])**2 + (by+file_vars[1][i])**2)
        var.append(field_mag)
  elif output_var == "div_bi":
    for i in range(len(file_vars[0])):
        divergence = np.abs(div(file_vars[0][i],file_vars[1][i],x,y))
        var.append(divergence)
  elif aia.is_filter(output_var):
    for i in range(len(file_vars[0])):
        resp = aia.response(output_var,file_vars[1][i])
        var.append(np.sqrt((file_vars[0][i])**2*resp))
  else:
    var = file_vars[0]
  return var

def extract_data_from_file(filename, file_var_names, display_interval, xl_ghost, xu_ghost, yl_ghost, yu_ghost, start_time, end_time, file_vec_name=None):
  input_file = open(filename)
  line = input_file.readline()
  while line[0] == '#' or len(line) == 0:
    line = input_file.readline()
  #determine grid size
  assert line == "xdim,ydim\n"
  dim = input_file.readline().split(',')
  xdim = int(dim[0])
  ydim = int(dim[1])

  xl = 0 + xl_ghost
  xu = xdim - xu_ghost
  yl = 0 + yl_ghost
  yu = ydim - yu_ghost

  #read in grid cell positions
  assert input_file.readline() == "pos_x\n"
  X = read_grid(input_file,xdim,ydim,xl,yl,xu,yu)/1.0e8
  assert input_file.readline() == "pos_y\n"
  Y = read_grid(input_file,xdim,ydim,xl,yl,xu,yu)/1.0e8

  assert input_file.readline() == "be_x\n"
  bx = read_grid(input_file,xdim,ydim,xl,yl,xu,yu)
  assert input_file.readline() == "be_y\n"
  by = read_grid(input_file,xdim,ydim,xl,yl,xu,yu)
  assert input_file.readline() == "be_z\n"
  bz = read_grid(input_file,xdim,ydim,xl,yl,xu,yu)
  
  file_vars = [ [] for v in file_var_names ]
  vec_x = []
  vec_y = []
  t = []
  output_number = 0
  line = ""

  time = -1.0
  while True:
    #read through to next time step (or file end)
    line = input_file.readline()
    while line and line[0:2] != "t=":
      line = input_file.readline()
    if not line:
      break
    time = float(line.split('=')[1])

    while time < start_time:
      line = input_file.readline()
      while line and line[0:2] != "t=":
        line = input_file.readline()
      if not line:
        break
      time = float(line.split('=')[1])
    if not line:
      break

    if end_time > 0.0 and time > end_time+1.0:
      break

    # var_found = False
    vars_found = np.array([ False for v in file_var_names ]) #all set to false initially
    vec_x_found = (file_vec_name == None or file_vec_name == "be")
    vec_y_found = (file_vec_name == None or file_vec_name == "be")
    # collect all necessary data at the current time step
    while not (np.all(vars_found) and vec_x_found and vec_y_found):
      line = input_file.readline()
      while line and not (np.any(line.rstrip() == file_var_names) or line.rstrip() == (str(file_vec_name)+"_x") or line.rstrip() == (str(file_vec_name)+"_y")):
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
      if line.rstrip() == (str(file_vec_name)+"_x"): 
        curr_data.append("vec_x")
        vec_x_found = True
      if line.rstrip() == (str(file_vec_name)+"_y"): 
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
  return xdim, ydim, X, Y, file_vars, vec_x, vec_y, bx, by, bz, t
