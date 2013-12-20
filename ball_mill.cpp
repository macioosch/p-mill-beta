#include <iostream>
#include <iomanip>
#include "ball_mill.h"

void type_1_init(parameters &params, std::vector<ball> &b)
{
    b.resize(1);
    b[0].x.resize(2);
    b[0].v.resize(2);

    b[0].a = 0.0;
    b[0].w = 0.0;
    b[0].x[0] = 0.0;
    b[0].x[1] = params.y0;
    b[0].v[0] = params.vx0;
    b[0].v[1] = params.vy0;
}

blas::vector<double> a_total(const parameters &params, const ball &b,
                             const double &t)
{
    static blas::vector<double> force (2);
    force[0] = 0.0;
    force[1] = -params.m * params.g;
    if (b.x[1] < 0.0)
        force[1] += (4/3.) * params.Er * sqrt(params.rr) * pow(-b.x[1], 1.5)
                + 2*params.beta * sqrt((5/3.)*params.mr*params.Er)
                * pow(params.rr*(-b.x[1]), 0.25) * b.v[1];
    return force / params.m;
}

void simulate(const parameters &params, std::vector<ball> &b)
{
    const int steps = std::ceil(params.tmax / params.dt);
    const int output_interval = steps / params.output_lines;
    blas::vector<double> x1(2), v1(2), a1(2),
                         x2(2), v2(2), a2(2),
                         x3(2), v3(2), a3(2),
                         x4(2), v4(2), a4(2);

    for (int i=0; i<steps; ++i) {
        x1 = b[0].x;
        v1 = b[0].v;
        a1 = a_total(params, b[0], params.dt*i);

        x2 = b[0].x = x1 + v1*0.5*params.dt;
        v2 = b[0].v = v1 + a1*0.5*params.dt;
        a2 = a_total(params, b[0], params.dt*(i+0.5));

        x3 = b[0].x = x1 + v2*0.5*params.dt;
        v3 = b[0].v = v1 + a2*0.5*params.dt;
        a3 = a_total(params, b[0], params.dt*(i+0.5));

        x4 = b[0].x = x1 + v3*params.dt;
        v4 = b[0].v = v1 + a3*params.dt;
        a4 = a_total(params, b[0], params.dt*(i+1));

        b[0].x = x1 + (params.dt/6.0)*(v1 + 2*v2 + 2*v3 + v4);
        b[0].v = v1 + (params.dt/6.0)*(a1 + 2*a2 + 2*a3 + a4);

        if (i % output_interval == 0)
            std::cout << std::setprecision(14)
                      << params.dt*(i+1)
                      << "\t" << b[0].x[1]
                      << "\t" << b[0].v[1]
                      << std::endl;
    }
}
