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

step = 3
delta_a = -np.pi/2

t = []
x1 = []
y1 = []
a1 = []
x2 = []
y2 = []
a2 = []

with open(argv[1], 'rb') as csvfile:
    csvreader = csv.reader(csvfile, delimiter='\t')
    rows = 0
    for row in csvreader:
        if len(row) == 13 and row[0][0] != "#":
            if rows % step == 0:
                t.append(float(row[0]))
                x1.append(float(row[1]))
                y1.append(float(row[2]))
                a1.append(float(row[5]))
                x2.append(float(row[7]))
                y2.append(float(row[8]))
                a2.append(float(row[11]))
            rows += 1

# A circle
circle_r = 0.01
circle_t = np.linspace(0, 2*np.pi, 32)
circle_x = circle_r * np.cos(circle_t)
circle_y = circle_r * np.sin(circle_t)

t = np.array(t)
x1 = np.array(x1)
y1 = np.array(y1)
a1 = np.array(a1)
x2 = np.array(x2)
y2 = np.array(y2)
a2 = np.array(a2)

#min_x = min(x) - 2*circle_r
#min_y = min(y) - 2*circle_r
#view_width = max(x) - min_x + 2*circle_r
#view_height = max(y) - min_y + 2*circle_r
#view_limit = max(view_width, view_height)
view_limit = 0.05

def plot_track(x, y, a):
    ax.plot(x + circle_r*np.cos(a + delta_a), y + circle_r*np.sin(a + delta_a), "#aaaaaa")
    ax.plot(x, y, "#666666", lw=2)

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
        xlim=(-view_limit, view_limit), ylim=(-view_limit, view_limit))
        #xlim=(min_x, min_x + view_limit), ylim=(min_y, min_y + view_limit))
plot_track(x1, y1, a1)
plot_track(x2, y2, a2)
line1, = ax.plot([], [], lw=2)
line2, = ax.plot([], [], lw=2)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

ax.set_xlabel("x [cm]")
ax.set_ylabel("y [cm]")
indicator_length = max(circle_r, 0.04*view_limit)

# initialization function: plot the background of each frame
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    return line1, line2, time_text

# animation function.  This is called sequentially
def animate(i):
    def draw_circle(line, x, y, a):
        data_x = np.append(circle_x + x[i],
                [None, x[i], x[i] + indicator_length*np.cos(a[i] + delta_a)])
        data_y = np.append(circle_y + y[i],
                [None, y[i], y[i] + indicator_length*np.sin(a[i] + delta_a)])
        line.set_data(data_x, data_y)
    draw_circle(line1, x1, y1, a1)
    draw_circle(line2, x2, y2, a2)
    time_text.set_text('time = {:.3f} s'.format(t[i]))
    pbar.update(i+1)
    return line1, line2, time_text

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
