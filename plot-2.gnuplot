#!/usr/bin/env gnuplot
set terminal pdf enhanced size 6,5
basename = "steel-2"
filename = "output/".basename.".csv"
set output "plots/".basename.".pdf"

ball_rho = 7.8e3
ball_r = 0.005
ball_m = (4/3.) * pi * ball_r**3 * ball_rho
ball_I = 2/5.*ball_m*ball_r**2
E = 2e11
G = 8e10
nu = 0.5*E/G - 1.0
Er = E / (2*(1-nu**2))

dn(x1, y1, x2, y2) = 2*ball_r - sqrt((x2-x1)**2 + (y2-y1)**2)
eb(dn) = dn < 0 ? 0 : 8.0/15 * Er * sqrt(ball_r) * dn**2.5

energy(w, vx, vy) = (ball_m*(vx**2+vy**2) + ball_I*w**2)/2

e2(w1, vx1, vy1, x1, y1, w2, vx2, vy2, x2, y2) = \
    energy(w1, vx1, vy1) + energy(w2, vx2, vy2) + eb(dn(x1, y1, x2, y2))
L(x, vx, y, vy, w) = ball_m*(x*vy - y*vx) + w*ball_I

x0 = -3e-6
y0 = 0.000
v0 = 0.1
w0 = 0e2
e0 = energy(w0, v0, 0)
eCM0 = e2(w0, 0.5*v0, 0, x0, y0, 0, -0.5*v0, 0, 2*ball_r, 0)
L0 = L(x0, v0, y0, 0, w0)

xCM0 = (0.01+x0)/2
yCM0 = -0.0005

set multiplot layout 3,1
unset key
#set format x "%.2f"
#set format y "%.3f"

set xrange [0:100]
set format y "%.1l · 10^{%L}"
set decimalsign ','
set xlabel "t [{/Symbol m}s]"

set ylabel "|1 - E/E_0|"
set yrange [-0.5e-8:3.5e-8]
set ytics 1e-8
plot filename u (1e6*$1):(abs(1 - e2($7, $4, $5, $2, $3, $13, $10, $11, $8, $9)/e0)) w l
#set ylabel "|1 - E_{CM}/E_{CM,0}|"
#plot filename u 1:(abs(1 - e2($7, $4-v0/2, $5, $2, $3, $13, $10-v0/2, $11, $8, $9)/eCM0)) w l
set ylabel "|1 - p_x/p_{x,0}|"
set yrange [-0.5e-14:2.0e-14]
set ytics 5e-15
plot filename u (1e6*$1):(abs(1 - ($4 + $10)/v0)) w l
#set ylabel "|p_y/p_{x,0}|"
#plot filename u 1:(abs(($5 + $11)/v0)) w l
set ylabel "x^'_{CM} [m]"
set yrange [-1e-13:1e-13]
set ytics 5e-14
plot filename u (1e6*$1):(1 - (($2 + $8)/2 - $1*v0/2)/xCM0) w l
#set ylabel "y^'_{CM}"
#plot filename u 1:(1 - ($3 + $9)/2/yCM0) w l
#set ylabel "L"
#plot filename u 1:(abs(1 - (L($2,$4,$3,$5,$7) + L($8,$10,$9,$11,$13)) / L0)) w l
