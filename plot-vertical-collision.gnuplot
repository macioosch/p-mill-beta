#!/usr/bin/env gnuplot
set terminal pdf enhanced size 6,6
set decimalsign ','
file = "output/steel-1-vertical.csv"
peaks = "output/steel-1-vertical-peaks.csv"
set output "plots/steel-1-vertical-collision.pdf"

g = 9.80665
h0 = 0.05
radius = 0.005
rho = 7.8e3
E = 2e11
G = 8e10
nu = 0.5*E/G - 1.0
Er = E / (2*(1-nu**2))
m = (4/3.) * pi * radius**3 * rho
energy(h, v) = 0.5*m*v**2 + m*g*h + (h > 0 ? 0 : (8.0/15)*Er*sqrt(radius)*(-h)**2.5)
h_static = -(9/16.)**(1/3.0) * ( g*m/(Er*sqrt(radius)) )**(2/3.0)
E0 = energy(h0, 0)
restitution = 0.7

set multiplot layout 3,1
set format xy "%.0f"
unset key
set grid

dt = 0.101
scale = 1e6
set xlabel sprintf("t - %.0f {/Symbol m}s", scale*dt)

set xrange [(0.10096-dt)*scale : (0.10104-dt)*scale]

set ylabel "y [{/Symbol m}m]"
plot file u (($1-dt)*scale):($3*1e6) w l

set ylabel "v_y [m s^{-1}]"
set format y "%.1f"
set yrange [-1.2:1.2]
plot file u (($1-dt)*scale):5 w l

set ylabel "E/E_0 - 1"
#set yrange [0.4:1.1]
set format y "%.0l · 10^{%L}"
set autoscale y
plot file u (($1-dt)*scale):(energy($3, $5)/E0 - 1) w l
