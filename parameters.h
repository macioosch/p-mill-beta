#ifndef MAIN_H
#define MAIN_H

struct parameters {
    /* type of the simulation, supported right now:
     * 1: one ball and a half-plane (with optional gravity),
     * TODO: 2: 2 balls interacting only with each other,
     * TODO: 3: N balls in a stationary or moving container.
     */
    int type;

    // for type 1 & 2: explicit setting of the sole ball's kinetic params
    double y0, vx0, vy0, w0, g;

    // additionally for type 2, params of ball 2:
    double w1, x0;

    // parameters of the balls set by the user
    double E, e, G, mu_r, mu_s, nu, rho, rb;
    // other parameters
    double tmax, dt, output_lines;

    // parameters computed by the program
    double Er, Gr, m, mr, rr, beta;
};

#endif // MAIN_H
