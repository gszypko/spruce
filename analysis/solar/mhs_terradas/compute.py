#!/usr/bin/env python3
# mhs_terradas/compute.py
# Computes normalized magnetohydrostatic solution for
# analytically-defined (at boundary) coronal hole or
# active region magnetic field structure
# Output (normalized.txt) should be fed to mhs_terradas/denormalize.py
# to generate a .state file in physical units
# Implementation of J. Terradas et al. 2022, A&A, http://arxiv.org/abs/2202.06800

import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from scipy.interpolate import LinearNDInterpolator
import multiprocessing as mp
import sys
import os
import csv

if len(sys.argv) != 3:
    print("Need to specify output directory, ACTIVE_REGION")
    exit()
out_directory = sys.argv[1]
ACTIVE_REGION = sys.argv[2] == "T"
# GAUSSIAN = sys.argv[3] == "T"
GAUSSIAN = False

NOGRAV_MODE = False

if not os.path.exists(out_directory):
    os.makedirs(out_directory)

DOMAIN_WIDTH = 8.0
DOMAIN_HEIGHT = 12.0
X_DIM = 202
Z_DIM = 302
# N_GHOST = 1

Z_NONUNIFORM_START = 302
Z_NONUNIFORM_FACTOR = 1.0
Z_NONUNIFORM_FACTOR_FACTOR = 1.0

if ACTIVE_REGION:
    X_1 = 2.25 #location of center for CH; location of one magnetic pole for AR
    X_2 = 1.75 #location of other magnetic pole for AR; not used for CH
else: #CORONAL_HOLE
    X_1 = -2.0 #location of center for CH; location of one magnetic pole for AR
    X_2 = 0.0 #location of other magnetic pole for AR; not used for CH
B_1_MULT = -1.0
B_2_MULT = 1.0

# In units of P_C and T_C, respectively
PRESS_AR = 40.0
TEMP_AR = 4.0
PRESS_CH = 0.25
TEMP_CH = 0.8

CHI = 20.0 # pressure scale height (along field lines)
# BETA = 0.004

W_0 = 0.2 # for gaussian profile
X_0 = 0.2 # for parabolic profile


def a_boundary_parabolic(x_0, b_0):
    def parabolic(x):
        return (abs(x) < x_0) * (b_0*x - b_0/(x_0**2)*x**3/3.0) \
        + (x >= x_0) * (2.0/3.0*b_0*x_0) \
        + (x <= -x_0) * (-2.0/3.0*b_0*x_0)
    return parabolic

def a_boundary_gaussian(w_0, b_0):
    def gaussian(x):
        return b_0*w_0*np.sqrt(np.pi)*0.5*sp.erf(x/w_0)
    return gaussian

def a_boundary(x):
    if GAUSSIAN: func = a_boundary_gaussian(W_0,1)
    else: func = a_boundary_parabolic(X_0,1)
    if ACTIVE_REGION: return B_1_MULT*func(x-X_1) + B_2_MULT*func(x-X_2)
    else: return func(x-X_1)


def compute_b(a, x_pos, z_pos):
    result_x = np.zeros_like(a)
    result_z = np.zeros_like(a)
    xdim = result_x.shape[0]
    zdim = result_x.shape[1]
    for i in range(xdim):
        for j in range(zdim):
            if i != 0 and i != xdim - 1 and j != zdim - 1:
                if j == 0:
                    result_z[i,j] = (a[i+1,j] - a[i-1,j])/(x_pos[i+1]-x_pos[i-1])
                else:
                    result_x[i,j] = -(a[i,j+1] - a[i,j-1])/(z_pos[j+1]-z_pos[j-1])
                    result_z[i,j] = (a[i+1,j] - a[i-1,j])/(x_pos[i+1]-x_pos[i-1])
    return result_x,result_z

def compute_a0_from_boundary(x,z,num_bins=1000):
    if GAUSSIAN:
        dt = 1.0/num_bins
        t_vals = np.linspace(0.0,1.0,num_bins,endpoint=False) + 0.5*dt
        result = 0.0
        for t in t_vals:
            xi = (1-t)/t
            with np.errstate(divide='ignore'):
                result += (a_boundary(xi)/((x-xi)**2 + z**2) + a_boundary(-xi)/((x+xi)**2 + z**2))*dt/t**2
        result *= z/np.pi
    else:
        assert GAUSSIAN == False
        b0 = 1.0
        a0 = (2.0/3.0)*X_0*b0
        x_eff = x - X_1
        with np.errstate(divide='ignore'):
            result = a0/(4.0*np.pi*X_0**3) * ( \
                z*(-3.0*x_eff**2 + 3.0*X_0**2 + z**2)*np.log( ((x_eff - X_0)**2 + z**2)/((x_eff + X_0)**2 + z**2) ) \
                + (2.0*x_eff**3 - 6.0*x_eff*X_0**2 - 6.0*x_eff*z**2 + 4.0*X_0**3)*np.arctan((x_eff - X_0)/(z)) \
                + (-2.0*x_eff**3 + 6.0*x_eff*X_0**2 + 6.0*x_eff*z**2 + 4.0*X_0**3)*np.arctan((x_eff + X_0)/(z)) \
                - 8.0*x_eff*X_0*z \
            )
            if ACTIVE_REGION:
                result *= B_1_MULT
                x_eff = x - X_2
                result += B_2_MULT*a0/(4.0*np.pi*X_0**3) * ( \
                    z*(-3.0*x_eff**2 + 3.0*X_0**2 + z**2)*np.log( ((x_eff - X_0)**2 + z**2)/((x_eff + X_0)**2 + z**2) ) \
                    + (2.0*x_eff**3 - 6.0*x_eff*X_0**2 - 6.0*x_eff*z**2 + 4.0*X_0**3)*np.arctan((x_eff - X_0)/(z)) \
                    + (-2.0*x_eff**3 + 6.0*x_eff*X_0**2 + 6.0*x_eff*z**2 + 4.0*X_0**3)*np.arctan((x_eff + X_0)/(z)) \
                    - 8.0*x_eff*X_0*z \
                )
    return (z>0.0)*(result) + (z==0.0)*a_boundary(x)

def press_profile_ch(a, a_ref, p_ch):
    if NOGRAV_MODE: return (1.0/p_ch - 1.0)*(a_ref/a_ref)**2 + 1.0
    else: return (1.0/p_ch - 1.0)*(a/a_ref)**2 + 1.0

def temp_profile_ch(a, a_ref, t_ch):
    return (1.0 - t_ch)*(a/a_ref)**2 + t_ch
    # return (1.0 - t_ch)*(a_ref/a_ref)**2 + t_ch

def press_profile_ar(a, a_ref, p_ar):
    if NOGRAV_MODE: return (1.0 - 1.0/p_ar)*(a_ref/a_ref)**2 + (1.0/p_ar)
    else: return (1.0 - 1.0/p_ar)*(a/a_ref)**2 + (1.0/p_ar)

def temp_profile_ar(a, a_ref, t_ar):
    return (t_ar - 1.0)*(a/a_ref)**2 + 1.0
    # return (t_ar - 1.0)*(a_ref/a_ref)**2 + 1.0

def press_profile(a,a_ref):
    if ACTIVE_REGION: return press_profile_ar(a,a_ref,PRESS_AR)
    else: return press_profile_ch(a,a_ref,PRESS_CH)

def temp_profile(a,a_ref):
    if ACTIVE_REGION: return temp_profile_ar(a,a_ref,TEMP_AR)
    else: return temp_profile_ch(a,a_ref,TEMP_CH)

def temperature(a, z, a_ref):
    result = temp_profile(a,a_ref)
    if CHI == 0: return result
    else: return result*np.exp(z/CHI)

def pressure(a, z, a_ref):
    if NOGRAV_MODE:
        result = press_profile(a,a_ref)
        return result*np.exp(0.0*z)
    else:
        result = press_profile(a,a_ref)
        if CHI == 0.0: exponent = -z/temp_profile(a,a_ref)
        else: exponent = CHI/temp_profile(a,a_ref)*(np.exp(-z/CHI)-1.0) 
        return result*np.exp(exponent)

def press_profile_ch_deriv(a, a_ref, p_ch):
    if NOGRAV_MODE: return 0.0*a
    else: return 2.0*(1.0/p_ch - 1.0)*(a/a_ref)

def temp_profile_ch_deriv(a, a_ref, t_ch):
    return 2.0*(1.0 - t_ch)*(a/a_ref)
    # return 0.0*a

def press_profile_ar_deriv(a, a_ref, p_ar):
    if NOGRAV_MODE: return 0.0*a
    else: return 2.0*(1.0 - 1.0/p_ar)*(a/a_ref)

def temp_profile_ar_deriv(a, a_ref, t_ar):
    return 2.0*(t_ar - 1.0)*(a/a_ref)
    # return 0.0*a

def press_profile_deriv(a,a_ref):
    if ACTIVE_REGION: return press_profile_ar_deriv(a,a_ref,PRESS_AR)
    else: return press_profile_ch_deriv(a,a_ref,PRESS_CH)

def temp_profile_deriv(a,a_ref):
    if ACTIVE_REGION: return temp_profile_ar_deriv(a,a_ref,TEMP_AR)
    else: return temp_profile_ch_deriv(a,a_ref,TEMP_CH)

def f(a,a_ref,z_var):
    if CHI == 0.0: z = z_var
    else: z = -CHI*(np.exp(-z_var/CHI) - 1.0)
    if NOGRAV_MODE: return 0.0*np.exp(-z/temp_profile(a,a_ref))*(press_profile_deriv(a,a_ref) + z*temp_profile_deriv(a,a_ref)*press_profile(a,a_ref)/(temp_profile(a,a_ref)**2))
    else: return np.exp(-z/temp_profile(a,a_ref))*(press_profile_deriv(a,a_ref) + z*temp_profile_deriv(a,a_ref)*press_profile(a,a_ref)/(temp_profile(a,a_ref)**2))

def compute_a1(f1,f2,t_bins,s_bins,x,z):
    dt = 1.0/t_bins
    ds = 1.0/s_bins
    sum = 0.0
    for i in range(t_bins):
        t = (i + 0.5)*dt
        for j in range(s_bins):
            s = (j + 0.5)*ds
            t_arg = (1.0-t)/t
            s_arg = (1.0-s)/s
            sum += f1[i,j]*np.log( ((x-t_arg)**2 + (z-s_arg)**2)/((x-t_arg)**2 + (z+s_arg)**2) ) +\
                    f2[i,j]*np.log( ((x+t_arg)**2 + (z-s_arg)**2)/((x+t_arg)**2 + (z+s_arg)**2) )
    return -1.0/(4.0*np.pi)*sum

def precompute_f(t_bins,s_bins,a_ref,res=1000):
    f1 = np.zeros((t_bins,s_bins))
    f2 = np.zeros((t_bins,s_bins))
    dt = 1.0/t_bins
    ds = 1.0/s_bins
    for i in range(t_bins):
        t = (i + 0.5)*dt
        for j in range(s_bins):
            s = (j + 0.5)*ds
            t_arg = (1.0-t)/t
            s_arg = (1.0-s)/s
            a0_1 = compute_a0_from_boundary(t_arg,s_arg,res)
            a0_2 = compute_a0_from_boundary(-t_arg,s_arg,res)
            f1[i,j] = f(a0_1,a_ref,s_arg)*dt*ds/(t**2 * s**2)
            f2[i,j] = f(a0_2,a_ref,s_arg)*dt*ds/(t**2 * s**2)
    return f1,f2

glob_dict = {}

def init_worker(precomps0,precomps1,precomps_shape,counter,total):
    glob_dict['precomps0'] = precomps0
    glob_dict['precomps1'] = precomps1
    glob_dict['precomps_shape'] = precomps_shape
    glob_dict['counter'] = counter
    glob_dict['total'] = total

def mp_compute_a1(x,z,t_bins=200,s_bins=200):
    precomps0 = np.frombuffer(glob_dict['precomps0']).reshape(glob_dict['precomps_shape'])
    precomps1 = np.frombuffer(glob_dict['precomps1']).reshape(glob_dict['precomps_shape'])
    result = compute_a1(precomps0,precomps1,t_bins,s_bins,x,z)
    glob_dict['counter'].value += 1
    counter = glob_dict['counter'].value
    if counter%int(glob_dict['total']/100) == 0:
        percent = int(round(counter/glob_dict['total']*100))
        print(f"\r[{'█'*percent}{' '*(100-percent)}]{percent}%",end='\r')
    return result

if __name__ == '__main__':

    # x = np.linspace(-DOMAIN_WIDTH*(1.0/2.0 + N_GHOST/X_DIM),DOMAIN_WIDTH*(1.0/2.0 + N_GHOST/X_DIM),X_DIM + 2*N_GHOST)
    # z = np.linspace(-N_GHOST*DOMAIN_HEIGHT/Z_DIM,DOMAIN_HEIGHT*(1.0 + N_GHOST/Z_DIM),Z_DIM + 2*N_GHOST)
    x = np.linspace(-DOMAIN_WIDTH/2.0,DOMAIN_WIDTH/2.0,X_DIM)
    z = np.linspace(0.0,DOMAIN_HEIGHT,Z_DIM)

    dx_1d = (x[-1]-x[0])/(x.shape[0]-1)*np.ones_like(x)
    dz_1d = (z[-1]-z[0])/(z.shape[0]-1)*np.ones_like(z)

    base_dz = dz_1d[0]
    curr_dz = base_dz
    curr_factor = Z_NONUNIFORM_FACTOR
    for k in range(Z_NONUNIFORM_START,Z_DIM):
        curr_dz *= curr_factor
        dz_1d[k] = curr_dz
        z[k] = z[k-1] + 0.5*dz_1d[k-1] + 0.5*curr_dz
        curr_factor *= Z_NONUNIFORM_FACTOR_FACTOR

    z_2d = np.outer(np.ones_like(x),z)
    x_2d = np.outer(x,np.ones_like(z))
    dx_2d = np.outer(dx_1d,np.ones_like(z))
    dz_2d = np.outer(np.ones_like(x),dz_1d)
    zero = np.zeros_like(x_2d)
    one = np.ones_like(x_2d)

    # Computing full-domain A and B from boundary
    full_a0 = compute_a0_from_boundary(x_2d,z_2d)
    if ACTIVE_REGION:
        a_ref = np.max(full_a0)
    else:
        a_ref = np.min(full_a0)
    full_b0 = compute_b(full_a0,x,z)
    if ACTIVE_REGION: streamplot_pow = 1.5
    else: streamplot_pow = 2.0
    full_b0_mag = np.sqrt(full_b0[0][:,1]**2 + full_b0[1][:,1]**2)**streamplot_pow
    full_b0_mag_norm = full_b0_mag/np.trapz(full_b0_mag)
    cdf = [np.trapz(full_b0_mag_norm[:(i+1)]) for i in range(len(full_b0_mag_norm))]
    if ACTIVE_REGION: stream_points = np.interp(np.linspace(0.55,0.95,num=17),cdf,x)
    else: stream_points = np.interp(np.linspace(0.05,0.95,num=17),cdf,x)

    plot_extent = (x[0]-0.5*dx_1d[0],x[-1]+0.5*dx_1d[-1],z[0]-0.5*dz_1d[0],z[-1]+0.5*dz_1d[-1])

    x_uniform = np.arange(min(x)+dx_1d[0], max(x), dx_1d[0])
    z_uniform = np.arange(min(z)+dz_1d[0], max(z), dz_1d[0])
    x_uniform_grid, z_uniform_grid = np.meshgrid(x_uniform, z_uniform)  # 2D grid for interpolation
    x_grid, z_grid = np.meshgrid(x, z)  # 2D grid for interpolation

    interp_x = LinearNDInterpolator(list(zip(x_grid.flatten(), z_grid.flatten())), full_b0[0].transpose().flatten(),fill_value=0.0)
    interp_z = LinearNDInterpolator(list(zip(x_grid.flatten(), z_grid.flatten())), full_b0[1].transpose().flatten(),fill_value=0.0)
    full_b0_x_uniform = interp_x(x_uniform_grid, z_uniform_grid)
    full_b0_z_uniform = interp_z(x_uniform_grid, z_uniform_grid)

    # plt.subplots()
    # plt.imshow(np.transpose(full_a0),origin='lower',extent=plot_extent)
    fig, ax = plt.subplots()
    im = NonUniformImage(ax, animated=False, origin='lower', extent=plot_extent,\
        interpolation='nearest')
    im.set_clim(vmin=np.min(full_a0),vmax=np.max(full_a0))
    im.set_data(x,z,np.transpose(full_a0))
    # im.set_data(x,z,np.transpose(full_b0[0]))
    # im.set_data(x,z,np.transpose(full_b0[1]))
    ax.add_image(im)
    # ax.set_aspect('equal')
    ax.set(xlim=(x[0],x[-1]), ylim=(z[0],z[-1]))
    fig.colorbar(im,label=r'$A/B_0h$')
    # ax.streamplot(x,z,full_b0[0].transpose(),full_b0[1].transpose(),start_points=np.column_stack((stream_points,np.zeros_like(stream_points))),color=(0.0,0.0,0.0,0.2),density=100,linewidth=1.0,arrowstyle='->',maxlength=100.0)
    ax.streamplot(x_uniform,z_uniform,\
        full_b0_x_uniform,full_b0_z_uniform,\
            start_points=np.column_stack((stream_points,z[1]*np.ones_like(stream_points))),\
                color=(0.0,0.0,0.0,0.2),broken_streamlines=False,linewidth=1.0,arrowstyle='->')
    plt.title("Potential Field")
    plt.xlabel(r'$x/h$')
    plt.ylabel(r'$z/h$')
    plt.tight_layout()
    plt.savefig(out_directory+"/potential.png")

    full_a1 = np.zeros_like(full_a0)
    t_bins = 200
    s_bins = 200
    print("Precomputing pressure gradients...")
    precomps = precompute_f(t_bins,s_bins,a_ref)
    num_el = precomps[0].shape[0] * precomps[0].shape[1]

    pc0 = mp.RawArray('d',num_el)
    pc0_np = np.frombuffer(pc0).reshape(precomps[0].shape)
    np.copyto(pc0_np,precomps[0])
    pc1 = mp.RawArray('d',num_el)
    pc1_np = np.frombuffer(pc1).reshape(precomps[0].shape)
    np.copyto(pc1_np,precomps[1])
    counter = mp.Value('i',0)
    total = full_a1.shape[0]*full_a1.shape[1]

    args_list = []
    for i in range(full_a1.shape[0]):
        for j in range(full_a1.shape[1]):
            args_list.append((x[i],z[j]))
    print("Computing perturbation field...")
    print(f"[{' '*100}]0%",end='\r')
    with mp.Pool(processes=2,initializer=init_worker, initargs=(pc0,pc1,precomps[0].shape,counter,total)) as pool:
        full_a1 = np.reshape(pool.starmap(mp_compute_a1, args_list),(X_DIM,Z_DIM))
    print(f"[{'█'*100}]100%",end='\r')
    print('',end=None)

    full_b1 = compute_b(full_a1,x,z)

    # Computing pressure and temperature according to Terradas profile
    temp = temperature(full_a0,z_2d,a_ref)
    press = np.zeros_like(temp)
    for i in range(press.shape[0]):
        for j in range(press.shape[1]):
            press[i,j] = pressure(full_a0[i,j],z_2d[i,j],a_ref)


    plt.figure(1)
    plt.imshow(np.transpose(temp),origin='lower',extent=plot_extent)
    plt.title("Temperature")
    plt.xlabel(r'$x/h$')
    plt.ylabel(r'$z/h$')
    plt.colorbar(label=r'$T/T_C$')
    plt.tight_layout()
    plt.savefig(out_directory+"/temperature.png")

    plt.figure(2)
    plt.imshow(np.transpose(press),origin='lower',extent=plot_extent)
    plt.title("Pressure")
    plt.xlabel(r'$x/h$')
    plt.ylabel(r'$z/h$')
    plt.colorbar(label=r'$p/p_{00}$')
    plt.tight_layout()
    plt.savefig(out_directory+"/pressure.png")

    with open(out_directory+"/normalized.txt", 'w', newline='') as f:
        writer = csv.writer(f)
        labels = ["#","CHI"]
        values = ["#",CHI]
        if ACTIVE_REGION:
            writer.writerow(["# ACTIVE REGION"])
            labels.extend(["X_1","X_2","B1_FACTOR","B2_FACTOR","PRESS_AR","TEMP_AR"])
            values.extend([X_1,X_2,B_1_MULT,B_2_MULT,PRESS_AR,TEMP_AR])
        else:
            writer.writerow(["# CORONAL HOLE"])
            labels.extend(["PRESS_CH","TEMP_CH","X_1"])
            values.extend([PRESS_CH,TEMP_CH,X_1])
        if GAUSSIAN:
            writer.writerow(["# GAUSSIAN PROFILE"])
            labels.extend(["W_0"])
            values.extend([W_0])
        else:
            writer.writerow(["# PARABOLIC PROFILE"])
            labels.extend(["X_0"])
            values.extend([X_0])
        writer.writerow(labels)
        writer.writerow(values)

        writer.writerow([str(X_DIM-2),str(Z_DIM-2)])
        names = ["d_x","d_y","pos_x","pos_y","press","temp","b0_x","b0_y","b1_x","b1_y"]
        grids = [dx_2d,dz_2d,x_2d,z_2d,press,temp,full_b0[0],full_b0[1],full_b1[0],full_b1[1]]
        for i in range(len(names)):
            writer.writerow([names[i]])
            writer.writerows(grids[i][1:-1,1:-1])
            # writer.writerows(grids[i][N_GHOST:-N_GHOST,N_GHOST:-N_GHOST])