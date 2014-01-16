#!/usr/bin/env python2
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
import progressbar as pb
from matplotlib import pyplot as plt
from matplotlib import animation
from sys import argv

step = 30
delta_a = -np.pi/2

t = []
x = []
y = []
a = []
N = -1

with open(argv[1], 'rb') as csvfile:
    csvreader = csv.reader(csvfile, delimiter='\t')
    rows = 0
    for row in csvreader:
        if not row[0].startswith("#"):
            if -1 == N:
                N = (len(row) - 1) / 6
                x = [ [] for i in range(N) ]
                y = [ [] for i in range(N) ]
                a = [ [] for i in range(N) ]
            if rows % step == 0:
                t.append(float(row[0]))
                for i in range(N):
                    x[i].append(float(row[1 + 6*i]))
                    y[i].append(float(row[2 + 6*i]))
                    a[i].append(float(row[5 + 6*i]))
            rows += 1

# A circle
circle_r = 0.005
circle_t = np.linspace(0, 2*np.pi, 32)
circle_x = circle_r * np.cos(circle_t)
circle_y = circle_r * np.sin(circle_t)

t = np.array(t)
x = [ np.array(i) for i in x ]
y = [ np.array(i) for i in y ]
a = [ np.array(i) for i in a ]

container_r = 0.04
view_limit = 1.1*container_r

def plot_track(x, y, a):
    ax.plot(x + circle_r*np.cos(a + delta_a), y + circle_r*np.sin(a + delta_a), "#aaaaaa")
    ax.plot(x, y, "#666666", lw=2)

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
        xlim=(-view_limit, view_limit), ylim=(-view_limit, view_limit))
        #xlim=(min_x, min_x + view_limit), ylim=(min_y, min_y + view_limit))
line = []
for xi, yi, ai in zip(x, y, a):
    plot_track(xi, yi, ai)
    line_i, = ax.plot([], [], lw=2)
    line.append(line_i)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

plt.plot(circle_x*container_r/circle_r, circle_y*container_r/circle_r)

ax.set_xlabel("x [cm]")
ax.set_ylabel("y [cm]")
indicator_length = max(circle_r, 0.04*view_limit)

# initialization function: plot the background of each frame
def init():
    for l in line:
        l.set_data([], [])
    time_text.set_text('')
    return line + [time_text]

# animation function.  This is called sequentially
def animate(i):
    def draw_circle(line, x, y, a):
        data_x = np.append(circle_x + x[i],
                [None, x[i], x[i] + indicator_length*np.cos(a[i] + delta_a)])
        data_y = np.append(circle_y + y[i],
                [None, y[i], y[i] + indicator_length*np.sin(a[i] + delta_a)])
        line.set_data(data_x, data_y)
    for l, xi, yi, ai in zip(line, x, y, a):
        draw_circle(l, xi, yi, ai)
    time_text.set_text('time = {:.3f} s'.format(t[i]))
    pbar.update(i+1)
    return line + [time_text]

# call the animator.  blit=True means only re-draw the parts that have changed.
widgets = ['Animating: ', pb.Percentage(), ' ', pb.Bar(), ' ', pb.ETA(), ' ']
pbar = pb.ProgressBar(widgets=widgets, maxval=len(t)).start()
anim = animation.FuncAnimation(fig, animate, init_func=init,
        frames=len(t), interval=1, blit=True)
pbar.finish()

# Save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
if len(argv) >= 3:
    anim.save(argv[2], fps=30, extra_args=['-vcodec', 'libx264'])
else:
    plt.show()
