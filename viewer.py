#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

input_file = open(sys.argv[1])
display_interval = float(sys.argv[2])

#determine grid size
dim = input_file.readline().split(',')
xdim = int(dim[0])
ydim = int(dim[1])

rho = []
t = []
output_number = 0

while True:
  line = input_file.readline()
  if not line:
    break
  time = float(line.split('=')[1])
  if output_number == 0 or time/display_interval >= output_number:
    output_number += 1
    t.append(time)
    this_rho = np.zeros([xdim,ydim], dtype=float)
    rows = input_file.readline().split(';')
    for i in range(xdim):
      this_rho[i] = np.asarray(rows[i].split(','))
    rho.append(this_rho)
  else:
    input_file.readline()
  
input_file.close()

# fig = plt.figure()
fig, ax = plt.subplots()

frame = 0
# im = plt.imshow(np.transpose(rho[frame]), animated=True, origin='lower', interpolation='bilinear')
# fig.colorbar(im)
im = ax.imshow(np.transpose(rho[frame]), animated=True, origin='lower', interpolation='nearest', vmin=0.0, vmax=(np.max(rho[0])+1.0))
ax.set(xlabel="x",ylabel="y",title="frame="+str(frame))
fig.colorbar(im)

def updatefig(*args):
    global frame
    # global im
    global ax
    frame = (frame + 1)%len(rho)
    im.set_array(np.transpose(rho[frame]))
    ax.set(xlabel="x",ylabel="y",title="t="+str(t[frame]))
    return im, ax

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=False)
plt.show()
