#include <iostream>
#include <iomanip>
#include <math.h>
#include "ball_mill.h"

void type_1_init(parameters &params, std::vector<ball> &b)
{
    b.resize(1);
    b[0].x.resize(2);
    b[0].v.resize(2);

    b[0].a = 0.0;
    b[0].w = params.w0;
    b[0].x[0] = 0.0;
    b[0].x[1] = params.y0;
    b[0].v[0] = params.vx0;
    b[0].v[1] = params.vy0;

    b[0].wall_collision = NULL;
}

blas::vector<double> a_total_1(const parameters &params, const ball &b)
{
    static blas::vector<double> force (2);
    static double delta_n, delta_t, fmax;
    force[0] = 0.0;
    force[1] = -params.m * params.g;
    if (b.x[1] < 0.0) {
        // the ball is colliding with the wall
        delta_n = -b.x[1];
        if (NULL != b.wall_collision)
            delta_t = b.x[0] + params.rb*std::sin(b.a-b.wall_collision->a0)
                    - b.wall_collision->point[0];
        else
            delta_t = 0.0;
        force[0] += 4*params.beta * sqrt(5/3.*params.m*params.Gr)
                * pow(params.rb*delta_n, 0.25) * (b.v[0] + b.w*(params.rb-delta_n))
                - 8.0 * params.Gr * sqrt(params.rb * delta_n) * delta_t;
        force[1] += (4/3.) * params.Er * sqrt(params.rb) * pow(delta_n, 1.5)
                + 2*params.beta * sqrt(5/3.*params.m*params.Er)
                * pow(params.rb*delta_n, 0.25) * b.v[1];
        fmax = fabs(force[1] * params.mu_s);
        if (fabs(force[0]) > fmax)
            force[0] = fmax * force[0] / fabs(force[0]);
    }
    return force / params.m;
}

void simulate_1(const parameters &params, std::vector<ball> &b)
{
    const int steps = std::ceil(params.tmax / params.dt);
    const int output_interval = steps / params.output_lines;
    blas::vector<double> x1(2), dx1(2), d2x1(2),
                         x2(2), dx2(2), d2x2(2),
                         x3(2), dx3(2), d2x3(2),
                         x4(2), dx4(2), d2x4(2);
    double a1, da1, d2a1,
           a2, da2, d2a2,
           a3, da3, d2a3,
           a4, da4, d2a4;

    for (int i=0; i<steps; ++i) {
        // RK4 INTEGRATION
        x1 = b[0].x;
        a1 = b[0].a;
        dx1 = b[0].v;
        da1 = b[0].w;
        d2x1 = a_total_1(params, b[0]);//, params.dt*i);
        if (x1[1] >= 0.0)
            d2a1 = 0.0;
        else
            d2a1 = (params.m*d2x1[0]*(params.rb+x1[1])
                    - fabs((4/3.)*params.Er*sqrt(params.rb)*pow(-x1[1], 1.5))*params.mu_r*params.rb*da1)
                    / (2/5. * params.m * pow(params.rb, 2));

        x2 = b[0].x = x1 + dx1*0.5*params.dt;
        a2 = b[0].a = a1 + da1*0.5*params.dt;
        dx2 = b[0].v = dx1 + d2x1*0.5*params.dt;
        da2 = b[0].w = da1 + d2a1*0.5*params.dt;
        d2x2 = a_total_1(params, b[0]);//, params.dt*(i+0.5));
        if (x2[1] >= 0.0)
            d2a2 = 0.0;
        else
            d2a2 = (params.m*d2x2[0]*(params.rb+x2[1])
                    - fabs((4/3.)*params.Er*sqrt(params.rb)*pow(-x2[1], 1.5))*params.mu_r*params.rb*da2)
                    / (2/5. * params.m * pow(params.rb, 2));

        x3 = b[0].x = x1 + dx2*0.5*params.dt;
        a3 = b[0].a = a1 + da2*0.5*params.dt;
        dx3 = b[0].v = dx1 + d2x2*0.5*params.dt;
        da3 = b[0].w = da1 + d2a2*0.5*params.dt;
        d2x3 = a_total_1(params, b[0]);//, params.dt*(i+0.5));
        if (x3[1] >= 0.0)
            d2a3 = 0.0;
        else
            d2a3 = (params.m*d2x3[0]*(params.rb+x3[1])
                    - fabs((4/3.)*params.Er*sqrt(params.rb)*pow(-x3[1], 1.5))*params.mu_r*params.rb*da3)
                    / (2/5. * params.m * pow(params.rb, 2));

        x4 = b[0].x = x1 + dx3*params.dt;
        a4 = b[0].a = a1 + da3*params.dt;
        dx4 = b[0].v = dx1 + d2x3*params.dt;
        da4 = b[0].w = da1 + d2a3*params.dt;
        d2x4 = a_total_1(params, b[0]);//, params.dt*(i+1));
        if (x4[1] >= 0.0)
            d2a4 = 0.0;
        else
            d2a4 = (params.m*d2x4[0]*(params.rb+x4[1])
                    - fabs((4/3.)*params.Er*sqrt(params.rb)*pow(-x4[1], 1.5))*params.mu_r*params.rb*da4)
                    / (2/5. * params.m * pow(params.rb, 2));

        b[0].x = x1 + (params.dt/6.0)*(dx1 + 2*dx2 + 2*dx3 + dx4);
        b[0].a = a1 + (params.dt/6.0)*(da1 + 2*da2 + 2*da3 + da4);
        b[0].v = dx1 + (params.dt/6.0)*(d2x1 + 2*d2x2 + 2*d2x3 + d2x4);
        b[0].w = da1 + (params.dt/6.0)*(d2a1 + 2*d2a2 + 2*d2a3 + d2a4);

        // checking for COLLISIONS
        if (b[0].x[1] <= 0.0 && NULL == b[0].wall_collision) {
            // a collision with wall just started
            b[0].wall_collision = new collision_with_wall;
            b[0].wall_collision->point = b[0].x;
            b[0].wall_collision->point[1] = 0.0;
            b[0].wall_collision->r0 =
                    b[0].wall_collision->point - b[0].x;
            b[0].wall_collision->a0 = b[0].a;
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
