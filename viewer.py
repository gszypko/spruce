#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

input_file = open(sys.argv[1])
display_interval = float(sys.argv[2])
output_var = sys.argv[3] #must be one of rho, temp, press

#determine grid size
dim = input_file.readline().split(',')
xdim = int(dim[0])
ydim = int(dim[1])

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
