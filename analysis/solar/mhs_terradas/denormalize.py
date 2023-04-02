#!/usr/bin/env python3
# mhs_terradas/denormalize.py
# Undoes unit normalization in output from mhs_terradas/compute.py
# and generates a .state file compatible with mhdtoy

import csv
import numpy as np
import os
import sys

BETA = 0.004
B_0 = -125.0*30.0 #G
PRESS_00 = BETA*B_0**2/(8.0*np.pi) #dyn cm^-2
T_C = 1.0e6 #K
G = 2.748e4 #cm s^-2
MU = 0.5*1.6726e-24 #g
K_B = 1.3807e-16 #erg K^-1
H = K_B*T_C/MU/G
GAMMA = 5.0/3.0
ION_MASS = 1.6726e-24 #g


if len(sys.argv) < 3:
    print("Need to specify input directory, output folder, and optionally output filename")
    exit()

in_directory = sys.argv[1]
if not os.path.exists(in_directory):
    print("Input directory does not exist")
    exit()

in_filename = in_directory + "/normalized.txt"
if not os.path.exists(in_filename):
    print("Input text file does not exist (must be named normalized.txt)")
    exit()

out_path = sys.argv[2]
if not os.path.exists(out_path):
    os.makedirs(out_path)

if len(sys.argv) == 4:
    out_filename = sys.argv[3]
else:
    out_filename = "mhs.state"

vars = {}
comments = []
xdim = 0
ydim = 0

with open(in_filename, newline='') as f:
    reader = csv.reader(f)
    line = next(reader)
    while line[0][0] == "#":
        comments.append(line)
        line = next(reader)
    dims = np.asarray(line).astype(int)
    xdim = dims[0]
    ydim = dims[1]
    for row in reader:
        curr_name = row[0]
        curr_grid = []
        for i in range(xdim):
            curr_grid.append(next(reader))
            assert(len(curr_grid[-1]) == ydim)
        vars[curr_name] = np.asarray(curr_grid).astype(np.float64)

zero = np.zeros((xdim,ydim))
one = np.ones((xdim,ydim))
vars["rho"] = vars["press"]/vars["temp"]
for name in ["grav_x","mom_x","mom_y","mom_z","be_x","be_y","be_z"]:
    vars[name] = zero
for name in ["grav_y"]:
    vars[name] = -1.0*one
for name in ["bi_z"]:
    vars[name] = 0.01*one
vars["bi_x"] = BETA/2*vars["b1_x"]
vars["bi_y"] = BETA/2*vars["b1_y"]
vars["be_x"] = vars["b0_x"]
vars["be_y"] = vars["b0_y"]
for name in ["press","b0_x","b0_y","b1_x","b1_y"]:
    vars.pop(name)

with open(out_path+"/"+out_filename, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(comments)
    writer.writerow(["#","DENORMALIZATION FACTORS:"])
    writer.writerow(["#","BETA","B_0","T_C","G","MU","K_B","GAMMA","ION_MASS"])
    writer.writerow(["#",str(BETA),str(B_0),str(T_C),str(G),str(MU),str(K_B),str(GAMMA),str(ION_MASS)])
    writer.writerow(["xdim","ydim"])
    writer.writerow([str(xdim),str(ydim)])
    writer.writerow(["ion_mass"])
    writer.writerow([str(ION_MASS)])
    writer.writerow(["adiabatic_index"])
    writer.writerow([str(GAMMA)])
    writer.writerow(["t=0"])
    names = ["d_x","d_y","pos_x","pos_y","rho","temp","mom_x","mom_y","mom_z","be_x","be_y","be_z","bi_x","bi_y","bi_z","grav_x","grav_y"]
    mults = [H,H,H,H,PRESS_00*MU/K_B/T_C,T_C,zero,zero,zero,B_0,B_0,B_0,B_0,B_0,B_0,G,G]
    for i in range(len(names)):
        name = names[i]
        mult = mults[i]
        writer.writerow([name])
        writer.writerows(mult*vars[name])
