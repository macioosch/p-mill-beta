#!/usr/bin/env gnuplot
set terminal pdf enhanced size 5,5
file = "output/steel-test.csv"
set output "plots/steel-test.pdf"

m = 0.03246312408709453
g = 9.80665
h0 = 0.02
Er = 1.098901098901099e8
Rr = 0.005
energy(h, v) = 0.5*m*v**2 + m*g*h + (h > 0 ? 0 : (8.0/15)*Er*sqrt(Rr)*(-h)**2.5)

set multiplot layout 3,1
set format xy "%.2e"
unset key
#set xrange [0.063:0.066]

set ylabel "y [m]"
plot file u 1:2 w l t "$y$"
set ylabel "v_y [m s^{-1}]"
plot file u 1:3 w l t "$v_y$"
set ylabel "E/E_0 - 1"
plot file u 1:(energy($2, $3)/energy(h0, 0) - 1) w l t "energy"
