#!/usr/bin/env gnuplot
set terminal pdf enhanced size 6,9
basename = "steel-1-horizontal"
filename = "output/".basename.".csv"
set output "plots/".basename.".pdf"

g = 9.80665
h0 = 0.0
ball_rho = 7.8e3
ball_r = 0.01
ball_m = (4/3.) * pi * ball_r**3 * ball_rho
ball_I = 2/5.*ball_m*ball_r**2

energy(y, w, vx, vy) = (ball_m*(vx**2+vy**2) + ball_I*w**2)/2 + ball_m*g*(y-h0)

set multiplot layout 6,1
unset key
#set format x "%.2f"
#set format y "%.3f"

set xlabel "x [m]"
set ylabel "y [m]"
plot filename u 2:3 w l
set ylabel "v_x [m s^{-1}]"
plot filename u 2:4 w l
set ylabel "|v_y| [m s^{-1}]"
plot filename u 2:(abs($5)) w l
set ylabel "{/Symbol w} [s^{-1}]"
plot filename u 2:7 w l
set ylabel "(mv^2 + I{/Symbol w}^2)/2 + mg(y-h_0) [J]"
plot filename u 2:(energy($3, $7, $4, $5)) w l
set ylabel "|{/Symbol w}r + v_x| [m s^{-1}]"
#set logscale y
plot filename u 2:(abs($7*ball_r + $4)) w l
