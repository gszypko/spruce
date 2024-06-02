import numpy as np
import math
import sys
import aia_response as aia
import scipy.ndimage

fullnames = {
  'rho': "Density",
  'temp': "Temperature",
  'press': "Pressure",
  'rad': "Radiative Loss Rate",
  'thermal_energy': "Thermal Energy Density",
  'be': "Background Magnetic Field",
  'bi': "Induced Magnetic Field",
  'bi_x': "X Induced Magnetic Field",
  'bi_y': "Y Induced Magnetic Field",
  'bi_z': "Z Induced Magnetic Field",
  'b_x': "X Magnetic Field",
  'b_y': "Y Magnetic Field",
  'b_z': "Z Magnetic Field",
  'b': "Magnetic Field",
  'v': "Velocity",
  'v_x': "X Velocity",
  'v_y': "Y Velocity",
  'v_z': "Z Velocity",
  'dt': "CFL Timestep",
  'dt_thermal': "Thermal Conduction Timestep",
  'dt_rad': "Radiative Losses Timestep",
  'n': "Density",
  'beta': "Plasma Beta",
  'div_be': "Background Field Divergence",
  'div_bi': "Induced Field Divergence",
  'b_mag': "Field Magnitude",
  'field_heating': "Field-Based Heating Rate",
  'roc': "Radius of Curvature",
  'xia_heat': "Heating Heuristic from Xia et al. 2014",
  'xia_comp': "Xia et al. 2014 Heating",
  'current_density': "Current Density Magnitude",
  'current_density_z': "Z Current Density",
  'anomalous_factor': "Anomalous Resistivity Factor",
  'thermal_conduction': "Thermal Conduction",
  'thermal_conduction_temp': "Thermal Conduction",
  'heating': "Heating Required for Balance",
  'field_heating_vs_heating': "Field Heating vs. Heating Required",
  'anomalous_diffusivity': "Anomalous Magnetic Diffusivity",
  'diffusivity': "Magnetic Diffusivity (Before Template)",
  'anomalous_template': "Anomalous Diffusivity Template",
  'lambda': "Lambda(T) Suppression Factor",
  'losses_over_rad': "Radiative and Conductive Losses",
  'joule_heating': "Joule Heating",
  'joule_heating_temp': "Joule Heating",
  'diff_time_scale': "Time Scale of Magnetic Diffusion"
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
  'b_x': r'G',
  'b_y': r'G',
  'b_z': r'G',
  'b': r'G',
  'v': r'cm s$^{-1}$',
  'v_x': r'cm s$^{-1}$',
  'v_y': r'cm s$^{-1}$',
  'v_z': r'cm s$^{-1}$',
  'dt': r's',
  'dt_thermal': r's',
  'dt_rad': r's',
  'n': r'cm$^{-3}$',
  'beta': r'Log',
  'div_b': r'G cm$^{-1}$',
  'div_bi': r'G cm$^{-1}$',
  'b_mag': r'G',
  'field_heating': r'erg cm$^{-3}$ s$^{-1}$',
  'roc': r'Mm',
  'xia_heat': r'erg cm$^{-3}$ s$^{-1}$',
  'xia_comp': r'Fractional Excess wrt Rad. Loss',
  'current_density': r'G cm$^{-1}$',
  'current_density_z': r'G cm$^{-1}$',
  'anomalous_factor': r'cm$^{-4}$',
  'thermal_conduction': r'erg cm$^{-3}$ s$^{-1}$',
  'thermal_conduction_temp': r'K s$^{-1}$',
  'heating': r'erg cm$^{-3}$ s$^{-1}$',
  'field_heating_vs_heating': r'Fractional Excess wrt. Heating Required',
  'anomalous_diffusivity': r'cm$^2$ s$^{-1}$',
  'diffusivity': r'cm$^2$ s$^{-1}$',
  'lambda': r'',
  'losses_over_rad': r'Fraction of Rad. Loss',
  'joule_heating': r'erg cm$^{-3}$ s$^{-1}$',
  'joule_heating_temp': r'K s$^{-1}$',
  'diff_time_scale': r's'
}

symlogthresholds = {
  ("div_bi","div_be","diffusivity"): 1e-25,
  ("dt","dt_thermal","dt_rad","b_mag","xia_heat"): 1e-8,
  ("xia_comp","field_heating_vs_heating"): 1e-2,
  ("current_density",): 1e-1,
  ("thermal_conduction",): 1e-6,
  ("thermal_conduction_temp",): 1e0,
  ("rad","heating"): 1e-8,
  ("lambda",): 1e-1,
  ("rho",): 1e-15,
  ("losses_over_rad",): 1e0,
  ("joule_heating",): 1e-4
}

file_vars_dict = {
  ("field_heating_vs_heating",): ("rad","thermal_conduction","field_heating"),
  ("heating",): ("rad","thermal_conduction"),
  ("anomalous_factor",): ("bi_x","bi_y","n"),
  ("current_density_z",): ("bi_x","bi_y"),
  ("current_density",): ("bi_x","bi_y","bi_z"),
  ("xia_comp",): ("bi_x","bi_y","bi_z","n","rad","thermal_conduction"),
  ("xia_heat",): ("bi_x","bi_y","bi_z","n"),
  ("roc",): ("bi_x","bi_y"),
  ("beta",): ("press","bi_x","bi_y","bi_z"),
  ("b_mag",): ("bi_x","bi_y","bi_z"),
  ("div_b","div_bi"): ("bi_x","bi_y"),
  ("lambda",): ("rad","n","temp"),
  ("losses_over_rad",): ("rad","thermal_conduction"),
  ("joule_heating_temp",): ("joule_heating","n"),
  ("thermal_conduction_temp",): ("thermal_conduction","n"),
  ("b_x",): ("bi_x",),
  ("b_y",): ("bi_y",),
  ("b_z",): ("bi_z",),
  ("diff_time_scale",): ("anomalous_diffusivity",)
}


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

def grad2_x(quantity,x,y):
  result = np.zeros_like(quantity)
  for i in range(1,len(x)-1):
    for j in range(1,len(y)-1):
      result[i,j] = (quantity[i+1,j] - 2.0*quantity[i,j] + quantity[i-1,j])/((0.5*(x[i+1]-x[i-1])*1.0e8)**2.0)
  return result

def grad2_y(quantity,x,y):
  result = np.zeros_like(quantity)
  for i in range(1,len(x)-1):
    for j in range(1,len(y)-1):
      result[i,j] = (quantity[i,j+1] - 2.0*quantity[i,j] + quantity[i,j-1])/((0.5*(y[j+1]-y[j-1])*1.0e8)**2.0)
  return result

def laplacian(quantity,x,y):
  return grad2_x(quantity,x,y) + grad2_y(quantity,x,y)

def lambda_opticallythin(t):
  regimes = {
    4.49: (0.0,1.0),
    4.97: (1.09e-31,2.0),
    5.67: (8.87e-17,-1.0),
    6.18: (1.90e-22,0.0),
    6.55: (3.53e-13,-1.5),
    6.90: (3.46e-25,1.0/3.0),
    7.63: (5.49e-16,-1.0),
    100.0: (1.96e-27,0.5)
  }
  for k in sorted(regimes.keys()):
    if np.log10(t) <= k:
      return regimes[k][0]*pow(t,regimes[k][1])
  assert False

def lambda_array(t):
  vfunc = np.vectorize(lambda t : lambda_opticallythin(t))
  return vfunc(t)

def apply_contour_computation(output_var, file_vars, x, y, bx, by, bz):
  var = []
  if output_var == "diff_time_scale":
    for i in range(len(file_vars[0])):
        diffusivity = file_vars[0][i]
        dx = np.maximum(x - np.roll(x,1), np.roll(x,-1) - x)*1.0e8
        dy = np.maximum(y - np.roll(y,1), np.roll(y,-1) - y)*1.0e8
        var.append((np.ma.array(np.outer(dx,dy))/np.ma.masked_less_equal(diffusivity,0.0)))
  elif output_var == "field_heating_vs_heating":
    for i in range(len(file_vars[0])):
        heating = - (file_vars[0][i] + file_vars[1][i])
        field_heating = file_vars[2][i]
        var.append((field_heating - heating)/heating)
  elif output_var == "b_x":
    for i in range(len(file_vars[0])):
        var.append(file_vars[0][i]+bx)
  elif output_var == "b_y":
    for i in range(len(file_vars[0])):
        var.append(file_vars[0][i]+by)
  elif output_var == "b_z":
    for i in range(len(file_vars[0])):
        var.append(file_vars[0][i]+bz)
  elif output_var == "joule_heating_temp":
    for i in range(len(file_vars[0])):
      result = file_vars[0][i]/(1.380649e-16*file_vars[1][i])
      # result = np.ma.masked_where(file_vars[2][i]<3.1e4,result)
      # result = np.divide(result,lambda_array(file_vars[2][i]),result,where=lambda_array(file_vars[2][i])>0.0)
      var.append(result)
  elif output_var == "thermal_conduction_temp":
    for i in range(len(file_vars[0])):
      result = file_vars[0][i]/(1.380649e-16*file_vars[1][i])
      # result = np.ma.masked_where(file_vars[2][i]<3.1e4,result)
      # result = np.divide(result,lambda_array(file_vars[2][i]),result,where=lambda_array(file_vars[2][i])>0.0)
      var.append(result)
  elif output_var == "lambda":
    for i in range(len(file_vars[0])):
      result = -file_vars[0][i]/(file_vars[1][i]**2)
      result = np.ma.masked_where(file_vars[2][i]<3.1e4,result)
      result = np.divide(result,lambda_array(file_vars[2][i]),result,where=lambda_array(file_vars[2][i])>0.0)
      var.append(result)
  elif output_var == "anomalous_factor":
    for i in range(len(file_vars[0])):
      mag_2d = np.sqrt((bx+file_vars[0][i])**2 + (by+file_vars[1][i])**2)
      grad_bx_mag = np.sqrt(grad_x(bx+file_vars[0][i],x,y)**2 + grad_y(bx+file_vars[0][i],x,y)**2)
      grad_by_mag = np.sqrt(grad_x(by+file_vars[0][1],x,y)**2 + grad_y(by+file_vars[1][i],x,y)**2)
      factor = np.abs(laplacian(np.sqrt(grad_bx_mag**2 + grad_by_mag**2)/(mag_2d),x,y))
      factor = factor**2
      factor = np.clip(1.0e50*factor,0.0,1.0)
      factor = scipy.ndimage.gaussian_filter(factor,sigma=3,mode='constant',radius=12)
      factor = np.clip(10.0*factor,0.0,1.0)
      factor = factor**1.5
      var.append(factor)
  elif output_var == "rad":
    for i in range(len(file_vars[0])):
      var.append(-file_vars[0][i])
  elif output_var == "losses_over_rad":
    for i in range(len(file_vars[0])):
      losses = -(file_vars[0][i]+file_vars[1][i])
      rad_losses = -file_vars[0][i]
      result = np.ma.masked_where(rad_losses <= 0.0, losses)
      var.append(np.divide(losses,rad_losses,where=np.greater(rad_losses,0.0),out=result))
  elif output_var == "heating":
    for i in range(len(file_vars[0])):
        var.append(-(file_vars[0][i] + file_vars[1][i]))
  elif output_var == "current_density_z":
    for i in range(len(file_vars[0])):
      curl_z = grad_x(by+file_vars[1][i],x,y) - grad_y(bx+file_vars[0][i],x,y)
      var.append(np.abs(3.0e10/(4.0*np.pi)*curl_z))
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
      n = file_vars[3][i]
      c_1 = 1.9e-22
      b_pow = 0.0
      n_pow = 2.0
      roc_pow = 0.0
      heat = np.ma.masked_invalid(c_1*(b_mag**b_pow)*(n**n_pow)/(roc**roc_pow))
      rad = np.ma.masked_invalid(file_vars[4][i])
      rad = np.ma.masked_where(rad==0.0,rad)
      cond = np.ma.masked_invalid(file_vars[5][i])
      cond = np.ma.masked_where(cond==0.0,cond)
      with np.errstate(divide='ignore',invalid='ignore'):
        ratio = np.ma.masked_invalid((heat - (rad - cond))/(rad - cond))
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
      mag_press = ((bx+file_vars[1][i])**2 + (by+file_vars[2][i])**2)/(8.0*np.pi)
      if np.all(mag_press == 0.0):
        var.append(np.ones_like(file_vars[0][i]))
      else:
        var.append(np.log10(file_vars[0][i]/np.ma.masked_where(mag_press == 0.0, mag_press)))
  elif output_var == "b_mag":
    for i in range(len(file_vars[0])):
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

def extract_tracer_particles_from_file(filename, display_interval, start_time, end_time, streamline_mode, label_list):
  input_file = open(filename)
  particle_sets = []
  particle_labels = []

  time = -1.0
  line = input_file.readline()
  output_number = 0
  t = []
  while True:
    #read through to next time step (or file end)
    while line and line[0:2] != "t=":
      line = input_file.readline()
    if not line:
      break
    time = float(line.split('=')[1])

    while time < start_time or not (output_number == 0 or display_interval == 0 or ((time - t[0])/display_interval >= output_number and not math.isinf(t[-1]))):
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

    curr_set = ([],[])
    curr_labels = []

    print(time)
    line = input_file.readline()
    while line and line[0:2] != "t=":
      comment_split = line.split('#')

      if len(comment_split) > 1:
        comment = comment_split[1].strip()
      else:
        comment = ""
      
      # print(comment)
      use_line = True
      if streamline_mode and not (comment == "" or comment[0] == 's'):
        use_line = False
      if len(label_list) != 0 and comment not in label_list:
        use_line = False
      # print(use_line)
      if use_line:
        # print("added")
        curr_labels.append(comment)
        curr_x,curr_y = comment_split[0].split(',')
        curr_set[0].append(float(''.join(curr_x.split())))
        curr_set[1].append(float(''.join(curr_y.split())))

      line = input_file.readline()
    output_number += 1
    t.append(time)
    
    particle_sets.append(curr_set)
    particle_labels.append(curr_labels)

    if not line:
      break

  # print(particle_sets)
  # print(particle_labels)
  input_file.close()
  return particle_sets, particle_labels

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
