#!/usr/bin/env python2

import csv
import numpy as np
from matplotlib import pyplot as plt
from sys import argv

# centimeters
circle_r = 0.5
container_r = 5.0

t = []
x = []
y = []
a = []
N = -1

sought_time = 2.0

with open(argv[1], 'rb') as csvfile:
    csvreader = csv.reader(csvfile, delimiter='\t')
    for row in csvreader:
        if not row[0].startswith("#"):
            time = float(row[0])
            if time >= sought_time:
                t = time
                N = (len(row) - 1) / 6
                x = [ float(row[1 + 6*i]) for i in range(N) ]
                y = [ float(row[2 + 6*i]) for i in range(N) ]
                a = [ float(row[5 + 6*i]) for i in range(N) ]
                break;

# A circle
circle_t = np.linspace(0, 2*np.pi, 32)
circle_x = circle_r * np.cos(circle_t)
circle_y = circle_r * np.sin(circle_t)

circle_hd_t = np.linspace(0, 2*np.pi, 128)
circle_hd_x = circle_r * np.cos(circle_hd_t)
circle_hd_y = circle_r * np.sin(circle_hd_t)

t = np.array(t)
x = 100*np.array(x)
y = 100*np.array(y)
a = 100*np.array(a)

view_limit = 1.1*container_r

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure(figsize=[3,3])

ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
        xlim=(-view_limit, view_limit), ylim=(-view_limit, view_limit))
ax.xaxis.grid(False)
ax.yaxis.grid(False)

if len(argv) >= 4:
    ax.set_frame_on(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

ax.plot(circle_hd_x*container_r/circle_r,
         circle_hd_y*container_r/circle_r, color="black")
for xi, yi, ai in zip(x, y, a):
    ax.plot(circle_x + xi, circle_y + yi, color="black")

ax.set_xlabel("x [cm]")
ax.set_ylabel("y [cm]")

if len(argv) >= 3:
    plt.savefig(argv[2], bbox_inches="tight")
else:
    plt.show()
