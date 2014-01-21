#!/usr/bin/env gnuplot
set terminal pdf enhanced size 6,8
file = "output/steel-1-vertical.csv"
set output "plots/steel-1-vertical-no-damping.pdf"

g = 9.80665
h0 = 0.01
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
restitution = 1.0

set multiplot layout 4,1
set format xy "%.2f"
unset key
set grid
set xlabel "t [s]"

set ylabel "y [m]"
set yrange [-.01:.06]
plot file u 1:3 w l

set ylabel "v_y [m s^{-1}]"
set yrange [-1.2:1.2]
plot file u 1:5 w l

set ylabel "|E/E_0 - 1|"
#set yrange [1e-5:3e0]
set yrange [1e-16:1e-5]
set logscale y
set format y " 10^{%L}"
unset mytics
plot file u 1:(abs(energy($3, $5)/E0 - 1)) w l

set ylabel "|E_i / E_{t,i} - 1|"
set format x "%.0f"
set xlabel "i"
set xtics 1
set xrange [0.001:7.999]
set yrange [2e-5:2e-1]
set autoscale y
plot peaks u 1:(abs($2 / (E0*restitution**(2*$1)) - 1))
#plot peaks u 1:(E0*restitution**(2*$1)) w l t "theory",\
#     peaks u 1:2 w l t "simulation"
