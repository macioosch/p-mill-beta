#!/usr/bin/env python3
"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""

import csv
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from sys import argv

step = 5

t = []
x = []
y = []
a = []

with open(argv[1], newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter='\t')
    rows = 0
    for row in csvreader:
        if len(row) == 7:
            if rows % step == 0:
                t.append(float(row[0]))
                x.append(float(row[1]))
                y.append(float(row[2]))
                a.append(float(row[5]))
            rows += 1

# A circle
circle_r = 0.005
circle_t = np.linspace(0, 2*np.pi, 32)
circle_x = circle_r * np.cos(circle_t)
circle_y = circle_r * np.sin(circle_t)

min_x = min(x) - circle_r
min_y = min(y) - circle_r
view_width = max(x) - min_x + circle_r
view_height = max(y) - min_y + circle_r
view_limit = max(view_width, view_height)

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(min_x, min_x + view_limit), ylim=(min_y, min_y + view_limit))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    ax = np.append(circle_x + x[i], [None, x[i], x[i] - 0.04*view_limit*np.sin(a[i])])
    ay = np.append(circle_y + y[i] + circle_r,
            [None, y[i] + circle_r, y[i] + circle_r + 0.04*view_limit*np.cos(a[i])])
    line.set_data(ax, ay)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
        frames=len(t), interval=1, blit=True)

plt.show()
