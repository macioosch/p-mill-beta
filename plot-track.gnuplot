#!/usr/bin/env gnuplot
# expects the "basename" variable as an option
set terminal pdf enhanced size 8,5
basename = "steel-1-horizontal"
filename = "output/".basename.".csv"
set output "plots/".basename.".pdf"

set multiplot layout 1,1
unset key

set xlabel "x [m]"
set ylabel "y [m]"
plot filename u 2:3 w l

