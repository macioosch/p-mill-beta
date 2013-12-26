#!/usr/bin/env gnuplot
set terminal pdf enhanced size 8,5
file = "output/steel-1-vertical.csv"
peaks = "output/steel-1-vertical-peaks.csv"
set output "plots/steel-1-vertical.pdf"

g = 9.80665
h0 = 0.02
radius = 0.01
rho = 7.75e3
E = 2e8
G = 7.93e7
nu = 0.5*E/G - 1.0
Er = E / (2*(1-nu**2))
m = (4/3.) * pi * radius**3 * rho
energy(h, v) = 0.5*m*v**2 + m*g*h + (h > 0 ? 0 : (8.0/15)*Er*sqrt(radius)*(-h)**2.5)
h_static = -(9/16.)**(1/3.0) * ( g*m/(Er*sqrt(radius)) )**(2/3.0)
E0 = energy(h0 - h_static, 0)
restitution = 0.8

set multiplot layout 4,1
set format xy "%.2e"
unset key
#set xrange [0.063:0.066]

set ylabel "y [m]"
plot file u 1:3 w l
set ylabel "y - h_s [m]"
plot file u 1:5 w l
set ylabel "E/E_0 - 1"
plot file u 1:(energy($3, $5)/energy(h0, 0) - 1) w l

set ylabel "E_i [J]"
set logscale y
#set yrange [-1:1]
#set xrange [0:10]
#plot peaks u 1:(-$2/(E0*restitution**(2*$1)) + 1) w p
set xrange [0:20]
#set yrange [1e-9:1e-2]
set key bottom left
plot peaks u 1:(E0*restitution**(2*$1)) w l t "theory",\
     peaks u 1:2 w l t "simulation"
