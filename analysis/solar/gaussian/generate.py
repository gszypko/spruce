#!/usr/bin/env python3
# gaussian/generate.py
# Generates a field-free gaussian density distribution

import scipy.special as sp
import numpy as np
import sys
import os
import csv

if len(sys.argv) != 3:
    print("Need to specify output directory and filename")
    exit()
out_directory = sys.argv[1]
filename = sys.argv[2]

if not os.path.exists(out_directory):
    os.makedirs(out_directory)

DOMAIN_WIDTH = 2.0e10
DOMAIN_HEIGHT = 2.0e10
X_DIM = 100
Z_DIM = 100

GAMMA = 5.0/3.0
ION_MASS = 1.6726e-24 #g

MAX_TEMP = 1.0e6
MAX_RHO = 1.0e-14
MIN_RHO = 1.0e-18

BE_X = 100.0/np.sqrt(2.0)
BE_Y = 100.0/np.sqrt(2.0)
BI_X = 0.0
BI_Y = 0.0

G = 0
#G = -5.0e4

x = np.linspace(-DOMAIN_WIDTH/2.0,DOMAIN_WIDTH/2.0,X_DIM)
z = np.linspace(-DOMAIN_HEIGHT/2.0,DOMAIN_HEIGHT/2.0,Z_DIM)
dx = (x[-1]-x[0])/(x.shape[0]-1)
dz = (z[-1]-z[0])/(z.shape[0]-1)
x_2d = np.outer(x,np.ones_like(z))
z_2d = np.outer(np.ones_like(x),z)

sigma = 0.1*DOMAIN_WIDTH
gaussian_1d = np.exp(-(x**2)/(2.0*sigma**2))
rho = (MAX_RHO-MIN_RHO)*np.outer(gaussian_1d,gaussian_1d) + MIN_RHO

zero = np.zeros_like(x_2d)
one = np.ones_like(x_2d)
with open(out_directory+"/"+filename, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["#SIMPLE ISOTHERMAL GAUSSIAN DENSITY"])
    writer.writerow(["#TEMP",str(MAX_TEMP)])
    writer.writerow(["#MAX_RHO",str(MAX_RHO)])
    writer.writerow(["#WIDTH",str(DOMAIN_WIDTH)])
    writer.writerow(["#HEIGHT",str(DOMAIN_HEIGHT)])
    writer.writerow(["xdim","ydim"])
    writer.writerow([str(X_DIM),str(Z_DIM)])
    writer.writerow(["ion_mass"])
    writer.writerow([str(ION_MASS)])
    writer.writerow(["adiabatic_index"])
    writer.writerow([str(GAMMA)])
    writer.writerow(["t=0"])
    names = ["d_x","d_y","pos_x","pos_y","rho","temp","mom_x","mom_y","be_x","be_y","bi_x","bi_y","grav_x","grav_y"]
    vars = [dx*one,dz*one,x_2d,z_2d,rho,MAX_TEMP*one,zero,zero,BE_X*one,BE_Y*one,BI_X*one,BI_Y*one,zero,G*one]
    for i in range(len(names)):
        writer.writerow([names[i]])
        writer.writerows(vars[i])
