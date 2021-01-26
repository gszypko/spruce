#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

vec_interval = 4 #interval between B field vectors to display

input_file = open(sys.argv[1])
display_interval = float(sys.argv[2])
output_var = sys.argv[3] #must be one of rho, temp, press

#determine grid size
dim = input_file.readline().split(',')
xdim = int(dim[0])
ydim = int(dim[1])

#read in static magnetic field
bx = np.zeros([xdim,ydim], dtype=float)
by = np.zeros([xdim,ydim], dtype=float)
rows_x = input_file.readline().split(';')
rows_y = input_file.readline().split(';')
for i in range(xdim):
  bx[i] = np.asarray(rows_x[i].split(','))
  by[i] = np.asarray(rows_y[i].split(','))

var = []
t = []
output_number = 0
line = ""

while True:
  #read through to next time step (or file end)
  line = input_file.readline()
  while line and line[0:2] != "t=":
    line = input_file.readline()
  if not line:
    break
  time = float(line.split('=')[1])

  #read through to correct variable
  line = input_file.readline()
  while line and line.rstrip() != output_var:
    if line[0:2] == "t=":
      sys.exit("Specified output variable not found in file")
    input_file.readline()
    line = input_file.readline()
  if not line:
    break

  if output_number == 0 or time/display_interval >= output_number:
    output_number += 1
    t.append(time)
    this_var = np.zeros([xdim,ydim], dtype=float)
    rows = input_file.readline().split(';')
    for i in range(xdim):
      this_var[i] = np.asarray(rows[i].split(','))
    var.append(this_var)
  
input_file.close()

# fig = plt.figure()
fig, ax = plt.subplots()

frame = 0
# im = plt.imshow(np.transpose(var[frame]), animated=True, origin='lower', interpolation='bilinear')
# fig.colorbar(im)
im = ax.imshow(np.transpose(var[frame]), animated=True, origin='lower', interpolation='nearest', vmin=0.0, vmax=(np.max(var[0])+1.0))
ax.set(xlabel="x",ylabel="y",title=output_var+", t="+str(t[frame]))
fig.colorbar(im)
X, Y = np.mgrid[0:xdim, 0:ydim]
ax.quiver(X[::vec_interval, ::vec_interval], \
  Y[::vec_interval, ::vec_interval], \
    bx[::vec_interval, ::vec_interval], \
      by[::vec_interval, ::vec_interval])


def updatefig(*args):
    global frame
    # global im
    global ax
    frame = (frame + 1)%len(var)
    im.set_array(np.transpose(var[frame]))
    ax.set(xlabel="x",ylabel="y",title=output_var+", t="+str(t[frame]))
    return im, ax

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=False)
plt.show()
