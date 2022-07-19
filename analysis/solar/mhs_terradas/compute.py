#!/usr/bin/env python3
# mhs_terradas/compute.py
# Computes normalized magnetohydrostatic solution for
# analytically-defined (at boundary) coronal hole or
# active region magnetic field structure
# Output (normalized.txt) should be fed to mhs_terradas/denormalize.py
# to generate a .state file in physical units
# Implementation of J., T., R., S., R., O., et al. 2022, A&A, http://arxiv.org/abs/2202.06800

import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import csv

if len(sys.argv) != 4:
    print("Need to specify output directory, ACTIVE_REGION, GAUSSIAN")
    exit()
out_directory = sys.argv[1]
ACTIVE_REGION = sys.argv[2] == "T"
GAUSSIAN = sys.argv[3] == "T"

if not os.path.exists(out_directory):
    os.makedirs(out_directory)

DOMAIN_WIDTH = 4.0
DOMAIN_HEIGHT = 4.0
X_DIM = 102
Z_DIM = 102
# N_GHOST = 1

X_1 = -0.25
X_2 = 0.25
B_1_MULT = 1.0
B_2_MULT = -1.0

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
    else: return func(x)


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
    dt = 1.0/num_bins
    t_vals = np.linspace(0.0,1.0,num_bins,endpoint=False) + 0.5*dt
    result = 0.0
    for t in t_vals:
        xi = (1-t)/t
        result += (a_boundary(xi)/((x-xi)**2 + z**2) + a_boundary(-xi)/((x+xi)**2 + z**2))*dt/t**2
    result *= z/np.pi
    return (z>0.0)*(result) + (z==0.0)*a_boundary(x)

def press_profile_ch(a, a_ref, p_ch):
    return (1.0/p_ch - 1.0)*(a/a_ref)**2 + 1.0

def temp_profile_ch(a, a_ref, t_ch):
    return (1.0 - t_ch)*(a/a_ref)**2 + t_ch

def press_profile_ar(a, a_ref, p_ar):
    return (1.0 - 1.0/p_ar)*(a/a_ref)**2 + (1.0/p_ar)

def temp_profile_ar(a, a_ref, t_ar):
    return (t_ar - 1.0)*(a/a_ref)**2 + 1.0

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
    result = press_profile(a,a_ref)
    if CHI == 0.0: exponent = -z/temp_profile(a,a_ref)
    else: exponent = CHI/temp_profile(a,a_ref)*(np.exp(-z/CHI)-1.0) 
    return result*np.exp(exponent)


def press_profile_ch_deriv(a, a_ref, p_ch):
    return 2.0*(1.0/p_ch - 1.0)*(a/a_ref)

def temp_profile_ch_deriv(a, a_ref, t_ch):
    return 2.0*(1.0 - t_ch)*(a/a_ref)

def press_profile_ar_deriv(a, a_ref, p_ar):
    return 2.0*(1.0 - 1.0/p_ar)*(a/a_ref)

def temp_profile_ar_deriv(a, a_ref, t_ar):
    return 2.0*(t_ar - 1.0)*(a/a_ref)

def press_profile_deriv(a,a_ref):
    if ACTIVE_REGION: return press_profile_ar_deriv(a,a_ref,PRESS_AR)
    else: return press_profile_ch_deriv(a,a_ref,PRESS_CH)

def temp_profile_deriv(a,a_ref):
    if ACTIVE_REGION: return temp_profile_ar_deriv(a,a_ref,TEMP_AR)
    else: return temp_profile_ch_deriv(a,a_ref,TEMP_CH)

def f(a,a_ref,z_var):
    if CHI == 0.0: z = z_var
    else: z = -CHI*(np.exp(-z_var/CHI) - 1.0)
    return np.exp(-z/temp_profile(a,a_ref))*(press_profile_deriv(a,a_ref) + z*temp_profile_deriv(a,a_ref)*press_profile(a,a_ref)/(temp_profile(a,a_ref)**2))

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
    print(str(-1.0/(4.0*np.pi)*sum))
    return -1.0/(4.0*np.pi)*sum

def precompute_f(t_bins,s_bins,a_ref):
    f1 = np.zeros((t_bins,s_bins))
    f2 = np.zeros((t_bins,s_bins))
    dt = 1.0/t_bins
    ds = 1.0/s_bins
    for i in range(t_bins):
        print(str(i)+"/"+str(t_bins))
        t = (i + 0.5)*dt
        for j in range(s_bins):
            s = (j + 0.5)*ds
            t_arg = (1.0-t)/t
            s_arg = (1.0-s)/s
            a0_1 = compute_a0_from_boundary(t_arg,s_arg,1000)
            a0_2 = compute_a0_from_boundary(-t_arg,s_arg,1000)
            f1[i,j] = f(a0_1,a_ref,s_arg)*dt*ds/(t**2 * s**2)
            f2[i,j] = f(a0_2,a_ref,s_arg)*dt*ds/(t**2 * s**2)
    return f1,f2

# x = np.linspace(-DOMAIN_WIDTH*(1.0/2.0 + N_GHOST/X_DIM),DOMAIN_WIDTH*(1.0/2.0 + N_GHOST/X_DIM),X_DIM + 2*N_GHOST)
# z = np.linspace(-N_GHOST*DOMAIN_HEIGHT/Z_DIM,DOMAIN_HEIGHT*(1.0 + N_GHOST/Z_DIM),Z_DIM + 2*N_GHOST)
x = np.linspace(-DOMAIN_WIDTH/2.0,DOMAIN_WIDTH/2.0,X_DIM)
z = np.linspace(0.0,DOMAIN_HEIGHT,Z_DIM)
dx = (x[-1]-x[0])/(x.shape[0]-1)
dz = (z[-1]-z[0])/(z.shape[0]-1)
x_2d = np.outer(x,np.ones_like(z))
z_2d = np.outer(np.ones_like(x),z)

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

plt.figure(3)
plt.imshow(np.transpose(full_a0),origin='lower',extent=(x[0]-0.5*dx,x[-1]+0.5*dx,z[0]-0.5*dz,z[-1]+0.5*dz))
plt.colorbar(label=r'$A/B_0h$')
plt.streamplot(x,z,full_b0[0].transpose(),full_b0[1].transpose(),start_points=np.column_stack((stream_points,np.zeros_like(stream_points))),color=(0.0,0.0,0.0,0.2),density=100,linewidth=1.0,arrowstyle='->',maxlength=100.0)
plt.title("Potential Field")
plt.xlabel(r'$x/h$')
plt.ylabel(r'$z/h$')
plt.tight_layout()
plt.savefig(out_directory+"/potential.png")

full_a1 = np.zeros_like(full_a0)
t_bins = 200
s_bins = 200
precomps = precompute_f(t_bins,s_bins,a_ref)
for i in range(full_a1.shape[0]):
    for j in range(full_a1.shape[1]):
        print(str(i) + " " + str(j))
        full_a1[i,j] = compute_a1(precomps[0],precomps[1],t_bins,s_bins,x[i],z[j])
full_b1 = compute_b(full_a1,x,z)

# Computing pressure and temperature according to Terradas profile
plt.figure(1)
temp = temperature(full_a0,z_2d,a_ref)
press = np.zeros_like(temp)
for i in range(press.shape[0]):
    for j in range(press.shape[1]):
        press[i,j] = pressure(full_a0[i,j],z_2d[i,j],a_ref)
plt.imshow(np.transpose(temp),origin='lower',extent=(x[0]-0.5*dx,x[-1]+0.5*dx,z[0]-0.5*dz,z[-1]+0.5*dz))
plt.title("Temperature")
plt.xlabel(r'$x/h$')
plt.ylabel(r'$z/h$')
plt.colorbar(label=r'$T/T_C$')
plt.tight_layout()
plt.savefig(out_directory+"/temperature.png")

plt.figure(2)
plt.imshow(np.transpose(press),origin='lower',extent=(x[0]-0.5*dx,x[-1]+0.5*dx,z[0]-0.5*dz,z[-1]+0.5*dz))
plt.title("Pressure")
plt.xlabel(r'$x/h$')
plt.ylabel(r'$z/h$')
plt.colorbar(label=r'$p/p_{00}$')
plt.tight_layout()
plt.savefig(out_directory+"/pressure.png")

zero = np.zeros_like(x_2d)
one = np.ones_like(x_2d)
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
        labels.extend(["PRESS_CH","TEMP_CH"])
        values.extend([PRESS_CH,TEMP_CH])
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
    grids = [dx*one,dz*one,x_2d,z_2d,press,temp,full_b0[0],full_b0[1],full_b1[0],full_b1[1]]
    for i in range(len(names)):
        writer.writerow([names[i]])
        writer.writerows(grids[i][1:-1,1:-1])
        # writer.writerows(grids[i][N_GHOST:-N_GHOST,N_GHOST:-N_GHOST])