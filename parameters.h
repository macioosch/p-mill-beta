#ifndef MAIN_H
#define MAIN_H

struct parameters {
    /* type of the simulation, supported right now:
     * 1: one ball and a half-plane (with optional gravity),
     * TODO: 2: 2 balls interacting only with each other,
     * TODO: 3: N balls in a stationary or moving container.
     */
    int type;

    // for type 1: explicit setting of the sole ball's kinetic params
    double y0, vx0, vy0, w0;

    // material parameters of the balls
    double E, e, G, mu, nu, rho;
    // geometric parameters of the balls
    double rb;
};

#endif // MAIN_H
