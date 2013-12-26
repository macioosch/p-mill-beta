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

    b[0].wall_collision = NULL;
}

blas::vector<double> a_total_1(const parameters &params, const ball &b)
{
    static blas::vector<double> force (2);
    static double delta_t;
    force[0] = 0.0;
    force[1] = -params.m * params.g;
    if (b.x[1] < 0.0) {
        // the ball is colliding with the wall
        if (NULL != b.wall_collision)
            delta_t = b.x[0] - b.wall_collision->x0[0];
        else
            delta_t = 0.0;
        force[0] += -8.0 * params.Gr * sqrt(params.rb * (-b.x[1])) * delta_t;
        force[1] += (4/3.) * params.Er * sqrt(params.rb) * pow(-b.x[1], 1.5)
                + 2*params.beta * sqrt((5/3.)*params.m*params.Er)
                * pow(params.rb*(-b.x[1]), 0.25) * b.v[1];
    }
    return force / params.m;
}

void simulate_1(const parameters &params, std::vector<ball> &b)
{
    const int steps = std::ceil(params.tmax / params.dt);
    const int output_interval = steps / params.output_lines;
    blas::vector<double> x1(2), v1(2), a1(2),
                         x2(2), v2(2), a2(2),
                         x3(2), v3(2), a3(2),
                         x4(2), v4(2), a4(2);

    for (int i=0; i<steps; ++i) {
        // RK4 INTEGRATION
        x1 = b[0].x;
        v1 = b[0].v;
        a1 = a_total_1(params, b[0]);//, params.dt*i);

        x2 = b[0].x = x1 + v1*0.5*params.dt;
        v2 = b[0].v = v1 + a1*0.5*params.dt;
        a2 = a_total_1(params, b[0]);//, params.dt*(i+0.5));

        x3 = b[0].x = x1 + v2*0.5*params.dt;
        v3 = b[0].v = v1 + a2*0.5*params.dt;
        a3 = a_total_1(params, b[0]);//, params.dt*(i+0.5));

        x4 = b[0].x = x1 + v3*params.dt;
        v4 = b[0].v = v1 + a3*params.dt;
        a4 = a_total_1(params, b[0]);//, params.dt*(i+1));

        b[0].x = x1 + (params.dt/6.0)*(v1 + 2*v2 + 2*v3 + v4);
        b[0].v = v1 + (params.dt/6.0)*(a1 + 2*a2 + 2*a3 + a4);

        // checking for COLLISIONS
        if (b[0].x[1] <= 0.0 && NULL == b[0].wall_collision) {
            // a collision with wall just started
            b[0].wall_collision = new collision;
            b[0].wall_collision->n.resize(2);
            b[0].wall_collision->n[0] = 0.0;
            b[0].wall_collision->n[1] = 1.0;
            b[0].wall_collision->x0 = b[0].x;
        } else if (b[0].x[1] > 0.0 && NULL != b[0].wall_collision) {
            // a collision with wall just ended
            delete b[0].wall_collision;
            b[0].wall_collision = NULL;
        }

        // occasional OUTPUT
        if (i % output_interval == 0)
            std::cout << std::setprecision(14)
                      << params.dt*(i+1)
                      << "\t" << b[0].x[0]
                      << "\t" << b[0].x[1]
                      << "\t" << b[0].v[0]
                      << "\t" << b[0].v[1]
                      << "\t" << b[0].a
                      << "\t" << b[0].w
                      << std::endl;
    }
}
