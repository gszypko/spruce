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
      result[i,j] = (a_x[i+1,j] - a_x[i-1,j])/((x[i+1]-x[i-1])*1.0e8) + (a_y[i,j+1] - a_y[i,j-1])/((y[j+1]-y[j-1])*1.0e8)
  return result

def grad_x(quantity,x,y):
  result = np.zeros_like(quantity)
  for i in range(1,len(x)-1):
    for j in range(1,len(y)-1):
      result[i,j] = (quantity[i+1,j] - quantity[i-1,j])/((x[i+1]-x[i-1])*1.0e8)
  return result

def grad_y(quantity,x,y):
  result = np.zeros_like(quantity)
  for i in range(1,len(x)-1):
    for j in range(1,len(y)-1):
      result[i,j] = (quantity[i,j+1] - quantity[i,j-1])/((y[j+1]-y[j-1])*1.0e8)
  return result

def apply_contour_computation(output_var, file_vars, x, y, bx, by, bz):
  var = []
  if output_var == "current_density_z":
    for i in range(len(file_vars[0])):
        curl_z = grad_x(by+file_vars[1][i],x,y) - grad_y(bx+file_vars[0][i],x,y)
        var.append(3.0e10/(4.0*np.pi)*curl_z)
  elif output_var == "current_density":
    for i in range(len(file_vars[0])):
        curl_z = grad_x(by+file_vars[1][i],x,y) - grad_y(bx+file_vars[0][i],x,y)
        curl_y = -grad_x(bz+file_vars[2][i],x,y)
        curl_x = grad_y(bz+file_vars[2][i],x,y)
        var.append(3.0e10/(4.0*np.pi)*np.sqrt(curl_x**2 + curl_y**2 + curl_z**2))
  elif output_var == "xia_comp":
    for i in range(len(file_vars[0])):
        mag_2d = np.sqrt((bx+file_vars[0][i])**2 + (by+file_vars[1][i])**2)
        b_hat_x = (bx+file_vars[0][i])/mag_2d
        b_hat_y = (by+file_vars[1][i])/mag_2d
        curv_x = b_hat_x*grad_x(b_hat_x,x,y) + b_hat_y*grad_y(b_hat_x,x,y)
        curv_y = b_hat_x*grad_x(b_hat_y,x,y) + b_hat_y*grad_y(b_hat_y,x,y)
        curv = (np.sqrt(curv_x*curv_x + curv_y*curv_y))
        with np.errstate(divide='ignore'):
          roc = 1.0/curv
        b_mag = np.sqrt((bx+file_vars[0][i])**2 + (by+file_vars[1][i])**2 + (bz+file_vars[2][i])**2)
        # b_mag = mag_2d
        n = file_vars[3][i]
        #0.29557738 0.75       1.085      1.215
        #8.9262941e-17 3.2500000e-01 1.7450000e+00 3.3500000e-01
        c_1 = 8.9262941e-17
        b_pow = 3.2500000e-01
        n_pow = 1.7450000e+00
        roc_pow = 3.3500000e-01
        # c_1 = 4.54*3.0e-11
        # b_pow = 1.75-0.25
        # n_pow = 0.125+0.825
        # roc_pow = 0.75 - 0.5
        # c_1 = 299
        # b_pow = 1.75
        # n_pow = 0.125
        # roc_pow = 0.75
        # c_1 = 2.0e-6
        # b_pow = 1.5
        # n_pow = 1.0
        # roc_pow = 0.75
        # c_1 = 8.7e-9
        # b_pow = 0.5
        # n_pow = 1.5
        # roc_pow = 1.0
        heat = np.ma.masked_invalid(c_1*(b_mag**b_pow)*(n**n_pow)/(roc**roc_pow))
        rad = np.ma.masked_invalid(file_vars[4][i])
        rad = np.ma.masked_where(rad==0.0,rad)
        with np.errstate(divide='ignore',invalid='ignore'):
          # ratio = np.ma.masked_invalid((heat - rad)/(1.38e-16*n))
          ratio = np.ma.masked_invalid((heat - rad)/(rad))
          var.append(ratio)
  elif output_var == "xia_heat":
    for i in range(len(file_vars[0])):
        mag_2d = np.sqrt((bx+file_vars[0][i])**2 + (by+file_vars[1][i])**2)
        b_hat_x = (bx+file_vars[0][i])/mag_2d
        b_hat_y = (by+file_vars[1][i])/mag_2d
        curv_x = b_hat_x*grad_x(b_hat_x,x,y) + b_hat_y*grad_y(b_hat_x,x,y)
        curv_y = b_hat_x*grad_x(b_hat_y,x,y) + b_hat_y*grad_y(b_hat_y,x,y)
        curv = (np.sqrt(curv_x*curv_x + curv_y*curv_y))
        with np.errstate(divide='ignore'):
          roc = 1.0/curv
        b_mag = np.sqrt((bx+file_vars[0][i])**2 + (by+file_vars[1][i])**2 + (bz+file_vars[2][i])**2)
        n = file_vars[3][i]
        c_1 = 4.54
        var.append(c_1*(b_mag**1.75)*(n**0.125)/(roc**0.75))
  elif output_var == "roc":
    for i in range(len(file_vars[0])):
        mag = np.sqrt((bx+file_vars[0][i])**2 + (by+file_vars[1][i])**2)
        b_hat_x = (bx+file_vars[0][i])/mag
        b_hat_y = (by+file_vars[1][i])/mag
        curv_x = b_hat_x*grad_x(b_hat_x,x,y) + b_hat_y*grad_y(b_hat_x,x,y)
        curv_y = b_hat_x*grad_x(b_hat_y,x,y) + b_hat_y*grad_y(b_hat_y,x,y)
        curv = (np.sqrt(curv_x*curv_x + curv_y*curv_y))
        with np.errstate(divide='ignore'):
          var.append(np.ma.masked_equal(np.where(curv > 0.0, 1.0/curv/1.0e8, 0.0),0.0))
  elif output_var == "beta":
    for i in range(len(file_vars[0])):
        # mag_press = ((bx+file_vars[1][i])**2 + (by+file_vars[2][i])**2 + (bz+file_vars[3][i])**2)/(8.0*np.pi)
        mag_press = ((bx+file_vars[1][i])**2 + (by+file_vars[2][i])**2)/(8.0*np.pi)
        var.append(np.log10(file_vars[0][i]/mag_press))
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
        var.append((file_vars[0][i])**2*resp)
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
