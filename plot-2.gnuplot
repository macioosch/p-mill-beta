#!/usr/bin/env gnuplot
set terminal pdf enhanced size 6,7
basename = "steel-2"
filename = "output/".basename.".csv"
set output "plots/".basename.".pdf"

ball_rho = 7.8e3
ball_r = 0.01
ball_m = (4/3.) * pi * ball_r**3 * ball_rho
ball_I = 2/5.*ball_m*ball_r**2

energy(w, vx, vy) = (ball_m*(vx**2+vy**2) + ball_I*w**2)/2
e2(w1, vx1, vy1, w2, vx2, vy2) = energy(w1, vx1, vy1) + energy(w2, vx2, vy2)

v0 = 0.1
w0 = 1e2
e0 = energy(w0, v0, 0)
eCM0 = e2(w0, 0.5*v0, 0, 0, -0.5*v0, 0)

xCM0 = 0.005
yCM0 = 0.005

set multiplot layout 6,1
unset key
#set format x "%.2f"
#set format y "%.3f"

set xlabel "t [s]"
set ylabel "|1 - E/E_0|"
plot filename u 1:(abs(1 - e2($7, $4, $5, $13, $10, $11)/e0)) w l
set ylabel "|1 - E_{CM}/E_{CM,0}|"
plot filename u 1:(abs(1 - e2($7, $4-v0/2, $5, $13, $10-v0/2, $11)/eCM0)) w l
set ylabel "|1 - p_x/p_{x,0}|"
plot filename u 1:(abs(1 - ($4 + $10)/v0)) w l
set ylabel "|p_y/p_{x,0}|"
plot filename u 1:(abs(($5 + $11)/v0)) w l
set ylabel "x^'_{CM}"
plot filename u 1:(1 - (($2 + $8)/2 - $1*v0/2)/xCM0) w l
set ylabel "y^'_{CM}"
plot filename u 1:(1 - ($3 + $9)/2/yCM0) w l
