#!/usr/bin/env gnuplot
set terminal pdf enhanced size 6,2
file = "output/steel-1-static.csv"
set output "plots/steel-1-static.pdf"

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
E0 = energy(h0 - h_static, 0)
restitution = 0.7

unset key
set grid
set xlabel "t [ms]"
set ylabel "{/Symbol d}_n [nm]"
set format x "%.1f"
set xrange [0:3]
plot file u ($1*1e3):(-$3*1e9) w l
