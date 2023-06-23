#!/usr/bin/env python3
# mhs_terradas/denormalize.py
# Undoes unit normalization in output from mhs_terradas/compute.py
# and generates a .state file compatible with mhdtoy

import csv
import numpy as np
import os
import sys

TRIXI_MODE = False
NOGRAV_MODE = False

BETA = 0.004
b_0 = [-125.0*30.0*0.1,20.0*30.0*0.1] #G
press_00 = [BETA*b_0[0]**2/(8.0*np.pi),BETA*b_0[1]**2/(8.0*np.pi)]
T_C = 1.0e6 #K
MU = 0.5*1.6726e-24 #g
K_B = 1.3807e-16 #erg K^-1
if NOGRAV_MODE:
    G = 0.0
    H = 6007869626.05265
else:
    G = 2.748e4 #cm s^-2
    H = K_B*T_C/MU/G
GAMMA = 5.0/3.0
ION_MASS = 1.6726e-24 #g

# OVERLAP = 20

if len(sys.argv) < 4:
    print("Need to specify AR directory, CH directory, output folder, and optionally output filename")
    exit()

in_directory_ar = sys.argv[1]
if not os.path.exists(in_directory_ar):
    print("AR Input directory does not exist")
    exit()
in_directory_ch = sys.argv[2]
if not os.path.exists(in_directory_ar):
    print("CH Input directory does not exist")
    exit()

in_filename_ar = in_directory_ar + "/normalized.txt"
if not os.path.exists(in_filename_ar):
    print("AR Input text file does not exist (must be named normalized.txt)")
    exit()
in_filename_ch = in_directory_ch + "/normalized.txt"
if not os.path.exists(in_filename_ch):
    print("CH Input text file does not exist (must be named normalized.txt)")
    exit()

out_path = sys.argv[3]
if not os.path.exists(out_path):
    os.makedirs(out_path)

if len(sys.argv) == 5:
    out_filename = sys.argv[4]
else:
    out_filename = "mhs.state"

vars = []
comments = []
xdim = []
ydim = []

for in_filename in [in_filename_ar,in_filename_ch]:
    this_vars = {}
    this_comments = []
    this_xdim = 0
    this_ydim = 0
    with open(in_filename, newline='') as f:
        reader = csv.reader(f)
        line = next(reader)
        while line[0][0] == "#":
            this_comments.append(line)
            line = next(reader)
        dims = np.asarray(line).astype(int)
        this_xdim = dims[0]
        this_ydim = dims[1]
        for row in reader:
            curr_name = row[0]
            curr_grid = []
            for i in range(this_xdim):
                curr_grid.append(next(reader))
                assert(len(curr_grid[-1]) == this_ydim)
            this_vars[curr_name] = np.asarray(curr_grid).astype(np.float64)
    vars.append(this_vars)
    comments.append(this_comments)
    xdim.append(this_xdim)
    ydim.append(this_ydim)

assert ydim[0] == ydim[1]


counter = 0
for this_var in vars:
    zero = np.zeros((xdim[counter],ydim[counter]))
    one = np.ones((xdim[counter],ydim[counter]))
    this_var["rho"] = this_var["press"]/this_var["temp"]
    for name in ["grav_x","mom_x","mom_y","mom_z","be_x","be_y","be_z"]:
        this_var[name] = zero
    for name in ["grav_y"]:
        this_var[name] = -1.0*one
    for name in ["bi_z"]:
        this_var[name] = 0.01*one
    this_var["bi_x"] = BETA/2*this_var["b1_x"]
    this_var["bi_y"] = BETA/2*this_var["b1_y"]
    this_var["be_x"] = this_var["b0_x"]
    this_var["be_y"] = this_var["b0_y"]
    for name in ["press","b0_x","b0_y","b1_x","b1_y"]:
        this_var.pop(name)
    counter += 1

assert xdim[0] == xdim[1] and ydim[0] == ydim[1]

with open(out_path+"/"+out_filename, 'w', newline='') as f:
    writer = csv.writer(f)
    if not TRIXI_MODE:
        for comment in comments: writer.writerows(comment)
        writer.writerow(["#","JOINED INTO INTERCHANGE RECONNECTION REGION"])
        writer.writerow(["#","DENORMALIZATION FACTORS:"])
        writer.writerow(["#","BETA","B_0_AR","B_0_CH","T_C","G","MU","K_B","GAMMA","ION_MASS"])
        writer.writerow(["#",str(BETA),str(b_0[0]),str(b_0[1]),str(T_C),str(G),str(MU),str(K_B),str(GAMMA),str(ION_MASS)])
    writer.writerow(["xdim","ydim"])
    writer.writerow([str(xdim[0]),str(ydim[0])])
    if not TRIXI_MODE:
        writer.writerow(["ion_mass"])
        writer.writerow([str(ION_MASS)])
    writer.writerow(["adiabatic_index"])
    writer.writerow([str(GAMMA)])
    if not TRIXI_MODE:
        writer.writerow(["t=0"])
    names = ["d_x","d_y","pos_x","pos_y","rho","temp","mom_x","mom_y","mom_z","be_x","be_y","be_z","bi_x","bi_y","bi_z","grav_x","grav_y"]
    mults = [[H,H,H,H,press_00[0]*MU/K_B/T_C,T_C,zero,zero,zero,b_0[0],b_0[0],b_0[0],b_0[0],b_0[0],b_0[0],G,G],\
                [H,H,H,H,press_00[1]*MU/K_B/T_C,T_C,zero,zero,zero,b_0[1],b_0[1],b_0[1],b_0[1],b_0[1],b_0[1],G,G]]
    for i in range(len(names)):
        name = names[i]
        ar_var = mults[0][i]*(vars[0][name])
        ch_var = mults[1][i]*(vars[1][name])
        if name in ["d_x","d_y","pos_x","pos_y","grav_x","grav_y"]:
            combined_var = ar_var
        elif name == "temp":
            ar_rho = mults[0][names.index("rho")]*(vars[0]["rho"])
            ch_rho = mults[1][names.index("rho")]*(vars[1]["rho"])
            combined_var = (ar_var*ar_rho + ch_var*ch_rho)/(ar_rho+ch_rho)
        else:
            combined_var =  ar_var + ch_var
        if TRIXI_MODE:            
            # writer.writerow([name])
            # writer.writerows(combined_var)
            # grids[n] = (grids[rho]/m_pd.m_ion_mass).max(m_pd.density_min);
            # grids[press] = 2*grids[n]*K_B*grids[temp].max(m_pd.temp_min);
            # grids[thermal_energy] = grids[press]/(m_pd.m_adiabatic_index - 1.0);
            if name in ["d_x","d_y","pos_x","pos_y","grav_x","grav_y","rho","mom_x","mom_y","mom_z"]:
                writer.writerow([name])
                writer.writerows(combined_var)
            elif name in ["temp"]:
                number_density = (mults[0][names.index("rho")]*(vars[0]["rho"]) + mults[1][names.index("rho")]*(vars[1]["rho"]))/ION_MASS
                pressure = 2.0*number_density*K_B*combined_var
                thermal_energy = pressure/(GAMMA - 1.0)
                writer.writerow(["thermal_energy"])
                writer.writerows(thermal_energy)
            elif name in ["be_x","be_y","be_z"]:
                external_name = name
                internal_name = external_name.replace("e","i")
                combined_var = combined_var + mults[0][names.index(internal_name)]*(vars[0][internal_name]) + mults[1][names.index(internal_name)]*(vars[1][internal_name])
                writer.writerow([name.replace("e","")])
                writer.writerows(combined_var/np.sqrt(4.0*np.pi))
        else:
            writer.writerow([name])
            writer.writerows(combined_var)

