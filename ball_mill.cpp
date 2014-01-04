#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>
#include "ball_mill.h"

void type_1_init(parameters &params, std::vector<ball> &b)
{
    b.resize(1);

    b[0].a = 0.0;
    b[0].w = params.w0;
    b[0].x[0] = 0.0;
    b[0].x[1] = params.y0;
    b[0].v[0] = params.vx0;
    b[0].v[1] = params.vy0;

    b[0].wall_delta_t = NULL;
}

Eigen::Vector2d a_total_1(const parameters &params, const ball &b)
{
    static Eigen::Vector2d force (2);
    static double delta_n, delta_t, fmax;
    static const double
            fnc_const = 4/3.*params.Er * sqrt(params.rb),
            fnd_const = 2*params.beta * sqrt(5/3.*params.m*params.Er) * pow(params.rb, 1/4.),
            ftc_const = -8*params.Gr * sqrt(params.rb),
            ftd_const = 4*params.beta * sqrt(5/3.*params.m*params.Gr) * pow(params.rb, 1/4.);
    force[0] = 0.0;
    force[1] = 0.0;
    if (b.x[1] < 0.0) {
        // the ball is colliding with the wall
        delta_n = -b.x[1];
        delta_t = 0.0;
        if (NULL != b.wall_delta_t)
            delta_t = b.wall_delta_t[0];
        force[0] += ftd_const * pow(delta_n, 1/4.) * (b.v[0] + b.w*params.rb)
                + ftc_const * sqrt(delta_n) * delta_t;
        force[1] += fnc_const * pow(delta_n, 3/2.)
                + fnd_const * pow(delta_n, 1/4.) * b.v[1];
        fmax = fabs(force[1] * params.mu_s);
        if (fabs(force[0]) > fmax)
            force[0] = fmax * force[0] / fabs(force[0]);
    }
    force[1] += -params.m * params.g;
    return force / params.m;
}

double fsign(const double &x)
{
    if (x > 0.0)
        return 1.0;
    else if (x < 0.0)
        return -1.0;
    return 0.0;
}

double angular_acceleration(const parameters &params, const double &a_tangential, const double &w,
                            const double &delta_n)
{
    static const double torque_const = 4/3.*params.Er * sqrt(params.rb);
    static const double moment_inertia = 2/5.*params.m * pow(params.rb, 2);
    return (params.m * a_tangential - fabs(torque_const*pow(delta_n, 3/2.)) * params.mu_r * fsign(w))
            * (params.rb-delta_n) / moment_inertia;
}

void simulate_1(const parameters &params, std::vector<ball> &b)
{
    const int steps = std::ceil(params.tmax / params.dt);
    int output_interval;
    Eigen::Vector2d x1(2), dx1(2), d2x1(2),
                    x2(2), dx2(2), d2x2(2),
                    x3(2), dx3(2), d2x3(2),
                    x4(2), dx4(2), d2x4(2), ac;
    double a1, da1, d2a1,
               da2, d2a2,
               da3, d2a3,
               da4, delta_t_0;
    bool reset_delta_t;

    if (steps > params.output_lines)
        output_interval = steps / params.output_lines;
    else
        output_interval = 1;

    for (int i=0; i<steps; ++i) {
        // RK4 INTEGRATION
        x1 = b[0].x;
        a1 = b[0].a;
        dx1 = b[0].v;
        da1 = b[0].w;
        d2x1 = a_total_1(params, b[0]);//, params.dt*i);
        delta_t_0 = 0.0;
        if (NULL != b[0].wall_delta_t && x2[1] < 0.0)
            delta_t_0 = b[0].wall_delta_t[0];
        d2a1 = 0.0;
        if (x1[1] < 0.0) d2a1 = angular_acceleration(params, d2x1[0], da1, -x1[1]);

        x2 = b[0].x = x1 + dx1*0.5*params.dt;
        b[0].a = a1 + da1*0.5*params.dt;
        dx2 = b[0].v = dx1 + d2x1*0.5*params.dt;
        da2 = b[0].w = da1 + d2a1*0.5*params.dt;
        if (NULL != b[0].wall_delta_t && x2[1] < 0.0)
            b[0].wall_delta_t[0] = delta_t_0 + dx1[0]*0.5*params.dt
                    + tan(da1*0.5*params.dt) * (params.rb+x2[1]);
        d2x2 = a_total_1(params, b[0]);//, params.dt*(i+0.5));
        d2a2 = 0.0;
        if (x2[1] < 0.0) d2a2 = angular_acceleration(params, d2x2[0], da2, -x2[1]);

        x3 = b[0].x = x1 + dx2*0.5*params.dt;
        b[0].a = a1 + da2*0.5*params.dt;
        dx3 = b[0].v = dx1 + d2x2*0.5*params.dt;
        da3 = b[0].w = da1 + d2a2*0.5*params.dt;
        if (NULL != b[0].wall_delta_t && x3[1] < 0.0)
            b[0].wall_delta_t[0] = delta_t_0 + dx2[0]*0.5*params.dt
                    + tan(da2*0.5*params.dt) * (params.rb+x3[1]);
        d2x3 = a_total_1(params, b[0]);//, params.dt*(i+0.5));
        d2a3 = 0.0;
        if (x3[1] < 0.0) d2a3 = angular_acceleration(params, d2x3[0], da3, -x3[1]);

        x4 = b[0].x = x1 + dx3*params.dt;
        b[0].a = a1 + da3*params.dt;
        dx4 = b[0].v = dx1 + d2x3*params.dt;
        da4 = b[0].w = da1 + d2a3*params.dt;
        if (NULL != b[0].wall_delta_t && x4[1] < 0.0)
            b[0].wall_delta_t[0] = delta_t_0 + dx3[0]*params.dt
                    + tan(da3*params.dt) * (params.rb+x4[1]);
        d2x4 = a_total_1(params, b[0]);//, params.dt*(i+1));
        //d2a4 = 0.0;
        //if (x4[1] < 0.0) d2a4 = angular_acceleration(params, d2x4[0], da4, -x4[1]);

        ac = d2x1 + 2*d2x2 + 2*d2x3 + d2x4;

        // checking for COLLISIONS
        reset_delta_t = false;
        if (b[0].x[1] <= 0.0 && NULL == b[0].wall_delta_t) {
            // a collision with wall just started
            b[0].wall_delta_t = new double(0.0);
        } else if (b[0].x[1] > 0.0 && NULL != b[0].wall_delta_t) {
            // a collision with wall just ended
            delete b[0].wall_delta_t;
            b[0].wall_delta_t = NULL;
        } else if (b[0].x[1] <= 0.0 && NULL != b[0].wall_delta_t) {
            // a collision continues: tangential force may need truncation
            if (fabs(ac[0]) > params.mu_s*fabs(ac[1]+params.g)) {
                ac[0] = params.mu_s*fabs(ac[1]+params.g) * ac[0]/fabs(ac[0]);
                reset_delta_t = true;
            }
        }

        b[0].x = x1 + params.dt/6.0*(dx1 + 2*dx2 + 2*dx3 + dx4);
        b[0].a = a1 + params.dt/6.0*(da1 + 2*da2 + 2*da3 + da4);
        b[0].v = dx1 + params.dt/6.0*ac;
        b[0].w = da1;
        if (b[0].x[1] < 0.0) b[0].w += params.dt*angular_acceleration(params, ac[0], (da1 + 2*da2 + 2*da3 + da4)/6.0, -b[0].x[1]);
        //b[0].w = da1 + params.dt/6.0*(d2a1 + 2*d2a2 + 2*d2a3 + d2a4);

        if (reset_delta_t)
            b[0].wall_delta_t[0] = 0.0;
        else if (NULL != b[0].wall_delta_t && b[0].x[1] < 0.0)
            b[0].wall_delta_t[0] = delta_t_0 + params.dt/6.0*(dx1 + 2*dx2 + 2*dx3 + dx4)[0]
                    + tan(params.dt/6.0*(da1 + 2*da2 + 2*da3 + da4)) * (params.rb+b[0].x[1]);

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
