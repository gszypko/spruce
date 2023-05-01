#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.image import NonUniformImage
from scipy.interpolate import LinearNDInterpolator
import argparse
from viewer_lib import *
import random as rnd

def compute_xia_heat(file_vars,c_1=343.4,b_pow=1.75,n_pow=0.125,roc_pow=0.75):
    global bx,by,bz
    i = 0
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
    return c_1*(b_mag**b_pow)*(n**n_pow)/(roc**roc_pow), b_mag, n, roc

# returns gradient with components corresponding to {c_1,b_pow,n_pow,roc_pow}
def compute_gradient(curr_state,xia_heat,rad_heat,b_mag,n,roc):
    result = np.array([0.0,0.0,0.0,0.0])
    arefinite = np.logical_and(np.isfinite(xia_heat),np.isfinite(rad_heat))
    param_coeff = [1.0/curr_state[0],np.log(b_mag,where=arefinite),np.log(n,where=arefinite),-np.log(roc,where=roc>0.0)]
    for i in range(1,4):
        result[i] = np.sum((param_coeff[i]*xia_heat)*(xia_heat - rad_heat),where=np.logical_and(arefinite,np.isfinite(param_coeff[i])))
    result *= 2.0/np.size(xia_heat)
    return result

def compute_loss(xia_heat,rad_heat):
    return np.sum(loss_map(xia_heat,rad_heat),where=np.logical_and(np.isfinite(xia_heat),np.isfinite(rad_heat)))/np.size(xia_heat)

def loss_map(xia_heat,rad_heat):
    result = np.zeros_like(rad_heat)
    np.log10(xia_heat/rad_heat,where=np.logical_and(rad_heat>0.0,xia_heat>0.0),out=result)
    return result**4

parser = argparse.ArgumentParser(description='View the output from mhdtoy.')
parser.add_argument('filename', help="the name of the file output from mhdtoy")
parser.add_argument('mode', help="the name of the file output from mhdtoy")
args = parser.parse_args()

if args.mode == "anneal":
    ANNEAL_MODE = True
elif args.mode == "gradient":
    ANNEAL_MODE = False
else:
    print("mode must be either anneal or gradient")
    exit()

start_time = 500
end_time = 700
num_ghost = 2
display_interval = 100
file_vec_name = None

xdim, ydim, X, Y, file_vars, vec_x, vec_y, bx, by, bz, t =\
  extract_data_from_file(args.filename, np.array(["bi_x","bi_y","bi_z","n","rad"]), display_interval, \
                         num_ghost, num_ghost, num_ghost, num_ghost, \
                            start_time, end_time, file_vec_name)

x = X[:,0]
y = Y[0,:]

# xia_heat, b_mag, n, roc = compute_xia_heat(file_vars)
rad_heat = file_vars[-1][0]
# state components: {c_1,b_pow,n_pow,roc_pow}
# state_init = np.array([343.4,1.75,0.125,0.75])
# state_init = np.array([299,1.75,0.125,0.75])
#9.78748685e-07 8.45000000e-01 1.30500000e+00 9.30000000e-01
state_init = np.array([3.1e-17,0.0,1.490000e+00,0.0])
# state_init = np.array([6.0e-7,1.5,1.0,0.75])
#6.62262243e-09 5.55000000e-01 1.56000000e+00 1.08500000e+00
# state_init = np.array([6.6e-9,0.5,1.5,1.0])
iter_limit = 20000

curr_state = state_init
print(f'initial state: {curr_state}')

if ANNEAL_MODE:
    temperature = 0.01
    temperature_factor = 0.9995
    min_temperature = 0.000001
    state_increments = [10.0**0.005,0.005,0.001,0.005]
    # state_increments = [10.0**0.02,0.05,0.05,0.05]
    # state_increments = [10.0**0.005,0.0,0.0,0.0]
    state_maxes = [np.inf,4.0,4.0,4.0]
    state_mins = [-np.inf,0.0,0.0,0.0]
    curr_xia = compute_xia_heat(file_vars,curr_state[0],curr_state[1],curr_state[2],curr_state[3])[0]
    curr_loss = compute_loss(curr_xia,rad_heat)
    best_loss = curr_loss
    best_state = curr_state
    loss_history = [curr_loss]
    for i in range(iter_limit):
        # attempt candidate modification
        cand_state = curr_state.copy()
        # cand_param = int(4.0*rnd.random())
        cand_param = 2*int(2.0*rnd.random())
        # cand_param = 0
        increase = rnd.random() >= 0.5
        if cand_param == 0:
            if increase:
                cand_state[0] *= state_increments[0]
            else:
                cand_state[0] /= state_increments[0]
        else:
            if increase:
                if cand_state[cand_param] < state_maxes[cand_param]:
                    cand_state[cand_param] += state_increments[cand_param]
                else:
                    continue
            else:
                if cand_state[cand_param] > state_mins[cand_param]:
                    cand_state[cand_param] -= state_increments[cand_param]
                else:
                    continue
        cand_xia = compute_xia_heat(file_vars,cand_state[0],cand_state[1],cand_state[2],cand_state[3])[0]
        cand_loss = compute_loss(cand_xia,rad_heat)

        verbose = i%100 == 0

        # keep candidate modification if an improvement (and sometimes keep it if not an improvement)
        thermal_chance = np.exp(-(cand_loss - curr_loss)/temperature)
        if verbose: print(f'Iteration {i}/{iter_limit}')
        if verbose: print(f'Candidate state: {cand_state}, loss: {cand_loss}')
        if verbose: print(f'Old state: {curr_state}, loss: {curr_loss}')
        if verbose: print(f'Temperature: {temperature}, Thermal jump chance: {thermal_chance}')
        if cand_loss < curr_loss or rnd.random() < thermal_chance:
            if cand_loss > curr_loss:
                if verbose: print(f'Modification accepted thermally')
            else:
                if verbose: print(f'Modification accepted normally')
            curr_state = cand_state
            curr_loss = cand_loss
            if curr_loss < min(loss_history):
                print(f'New global minimum found. Loss: {curr_loss}, State: {curr_state}')
                best_loss = curr_loss
                best_state = curr_state
            loss_history.append(curr_loss)
        else:
            if verbose: print(f'Modification rejected')
        if verbose: print()
        temperature *= temperature_factor
        if temperature < min_temperature:
            print("Minimum temperature reached.")
            break

    print(f'best state: {best_state}, loss: {best_loss}, initial state: {state_init}')

    xia_heat, b_mag, n, roc = compute_xia_heat(file_vars,best_state[0],best_state[1],best_state[2],best_state[3])

    plt.figure()
    plt.plot(loss_history)
    plt.yscale('log')
    plt.title("Loss History")
    plt.xlabel("Iteration")
    plt.ylabel("Loss")

else:
    learning_rate = np.array([0.0,1.0e-6,1.0e-6,1.0e-6])
    for i in range(iter_limit):
        xia_heat, b_mag, n, roc = compute_xia_heat(file_vars,curr_state[0],curr_state[1],curr_state[2],curr_state[3])
        loss = compute_loss(xia_heat,rad_heat)
        curr_grad = compute_gradient(curr_state,xia_heat,rad_heat,b_mag,n,roc)
        curr_state = curr_state - curr_grad*learning_rate
        if i%50 == 0:
            print(f'loss: {loss}, gradient: {curr_grad}, new_state: {curr_state}')

    print(f'final state: {curr_state}')

    xia_heat, b_mag, n, roc = compute_xia_heat(file_vars,curr_state[0],curr_state[1],curr_state[2],curr_state[3])

xia_max = np.nanmax(xia_heat)
xia_min = np.nanmin(xia_heat)
rad_max = np.nanmax(rad_heat)
rad_min = np.nanmin(rad_heat)
heat_max = max(xia_max,rad_max)
heat_min = min(xia_min,rad_min)


plt.figure()
plt.imshow(np.transpose(xia_heat), origin='lower', interpolation='nearest',norm=matplotlib.colors.SymLogNorm(linthresh=1e-4, base=10))
plt.title("xia_heat")
plt.clim(vmin=heat_min,vmax=heat_max)
plt.colorbar()
plt.figure()
plt.imshow(np.transpose(rad_heat), origin='lower', interpolation='nearest',norm=matplotlib.colors.SymLogNorm(linthresh=1e-4, base=10))
plt.title("rad_heat")
plt.clim(vmin=heat_min,vmax=heat_max)
plt.colorbar()
plt.figure()
plt.imshow(np.transpose(loss_map(xia_heat,rad_heat)), origin='lower', interpolation='nearest',norm=matplotlib.colors.LogNorm())
plt.title("loss")
plt.colorbar()
plt.show()

