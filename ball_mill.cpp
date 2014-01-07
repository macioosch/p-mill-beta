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

void type_2_init(parameters &params, std::vector<ball> &b)
{
    b.resize(2);

    b[0].a = 0.0;
    b[0].w = params.w0;
    b[0].x[0] = -params.x0;
    b[0].x[1] = params.y0;
    b[0].v[0] = params.vx0;
    b[0].v[1] = 0.0;
    b[0].wall_delta_t = NULL;

    b[1].a = 0.0;
    b[1].w = params.w1;
    b[1].x[0] = 2*params.rb;
    b[1].x[1] = 0.0;
    b[1].v[0] = 0.0;
    b[1].v[1] = 0.0;
    b[1].wall_delta_t = NULL;
}

Eigen::Vector2d a_total_1(const parameters &params, const ball &b)
{
    static Eigen::Vector2d force;
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

double angular_acceleration(const parameters &params, const double &a_tangential,
                            const double &w, const double &delta_n)
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
    Eigen::Vector2d x1, dx1, d2x1,
                    x2, dx2, d2x2,
                    x3, dx3, d2x3,
                    x4, dx4, d2x4, ac;
    double a1, da1, d2a1,
               da2, d2a2,
               da3, d2a3,
               da4, delta_t_0;
    bool reset_delta_t;

    if (steps > params.output_lines)
        output_interval = steps / params.output_lines;
    else
        output_interval = 1;

    for (int i=0; i<=steps; ++i) {
        // occasional OUTPUT
        if (i % output_interval == 0)
            std::cout << std::setprecision(14)
                      << params.dt*(i)
                      << "\t" << b[0].x[0]
                      << "\t" << b[0].x[1]
                      << "\t" << b[0].v[0]
                      << "\t" << b[0].v[1]
                      << "\t" << b[0].a
                      << "\t" << b[0].w
                      << std::endl;

        // RK4 INTEGRATION
        x1 = b[0].x;
        a1 = b[0].a;
        dx1 = b[0].v;
        da1 = b[0].w;
        d2x1 = a_total_1(params, b[0]);
        delta_t_0 = 0.0;
        if (NULL != b[0].wall_delta_t && x1[1] < 0.0)
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
        d2x2 = a_total_1(params, b[0]);
        d2a2 = 0.0;
        if (x2[1] < 0.0) d2a2 = angular_acceleration(params, d2x2[0], da2, -x2[1]);

        x3 = b[0].x = x1 + dx2*0.5*params.dt;
        b[0].a = a1 + da2*0.5*params.dt;
        dx3 = b[0].v = dx1 + d2x2*0.5*params.dt;
        da3 = b[0].w = da1 + d2a2*0.5*params.dt;
        if (NULL != b[0].wall_delta_t && x3[1] < 0.0)
            b[0].wall_delta_t[0] = delta_t_0 + dx2[0]*0.5*params.dt
                    + tan(da2*0.5*params.dt) * (params.rb+x3[1]);
        d2x3 = a_total_1(params, b[0]);
        d2a3 = 0.0;
        if (x3[1] < 0.0) d2a3 = angular_acceleration(params, d2x3[0], da3, -x3[1]);

        x4 = b[0].x = x1 + dx3*params.dt;
        b[0].a = a1 + da3*params.dt;
        dx4 = b[0].v = dx1 + d2x3*params.dt;
        da4 = b[0].w = da1 + d2a3*params.dt;
        if (NULL != b[0].wall_delta_t && x4[1] < 0.0)
            b[0].wall_delta_t[0] = delta_t_0 + dx3[0]*params.dt
                    + tan(da3*params.dt) * (params.rb+x4[1]);
        d2x4 = a_total_1(params, b[0]);
        //d2a4 = 0.0;
        //if (x4[1] < 0.0) d2a4 = angular_acceleration(params, d2x4[0], da4, -x4[1]);

        b[0].x = x1 + params.dt/6.0*(dx1 + 2*dx2 + 2*dx3 + dx4);

        ac = (d2x1 + 2*d2x2 + 2*d2x3 + d2x4)/6.0;

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

        b[0].a = a1 + params.dt/6.0*(da1 + 2*da2 + 2*da3 + da4);
        b[0].v = dx1 + params.dt*ac;
        b[0].w = da1;
        if (b[0].x[1] < 0.0) b[0].w += params.dt*angular_acceleration(params, ac[0], (da1 + 2*da2 + 2*da3 + da4)/6.0, -b[0].x[1]);
        //b[0].w = da1 + params.dt/6.0*(d2a1 + 2*d2a2 + 2*d2a3 + d2a4);

        if (reset_delta_t)
            b[0].wall_delta_t[0] = 0.0;
        else if (NULL != b[0].wall_delta_t && b[0].x[1] < 0.0)
            b[0].wall_delta_t[0] = delta_t_0 + params.dt/6.0*(dx1 + 2*dx2 + 2*dx3 + dx4)[0]
                    + tan(params.dt/6.0*(da1 + 2*da2 + 2*da3 + da4)) * (params.rb+b[0].x[1]);
    }
}

double cross2d(const Eigen::Vector2d &v1, const Eigen::Vector2d &v2)
{
    return v1[0]*v2[1] - v1[1]*v2[0];
}

Eigen::Matrix2d a_total_2(const parameters &params, std::vector<ball> &b,
                          Eigen::Vector2d *delta_t_given)
{
    static Eigen::Matrix2d acceleration_matrix, cross_matrix;
    static Eigen::Vector2d force, force_n, force_t, r01, v01, acceleration,
            n_dir, t_dir, delta_t;
    static double delta_n, fmax;
    static const double
            fnc_const = 4/3.*params.Er * sqrt(params.rb),
            fnd_const = 2*params.beta * sqrt(5/3.*params.m*params.Er) * pow(params.rb, 1/4.),
            ftc_const = -8*params.Gr * sqrt(params.rb),
            ftd_const = 4*params.beta * sqrt(5/3.*params.m*params.Gr) * pow(params.rb, 1/4.);
    cross_matrix << 0, 1, -1, 0;

    force[0] = force[1] = 0.0;
    r01 = b[1].x - b[0].x;
    if (r01.norm() < 2*params.rb) {
        // the balls are colliding
        delta_n = 2*params.rb - r01.norm();
        delta_t << 0, 0;
        if (NULL != delta_t_given)
            delta_t = delta_t_given[0];
        n_dir = r01.normalized();
        v01 = b[1].v - b[0].v + (b[0].w + b[1].w)*params.rb*cross_matrix*n_dir;
        if (n_dir[0]*v01[1] - n_dir[1]*v01[0] != 0.0)
            t_dir << n_dir[1] * (n_dir[1]*v01[0] - n_dir[0]*v01[1]),
                     n_dir[0] * (n_dir[0]*v01[1] - n_dir[1]*v01[0]);
        else if (n_dir[0]*delta_t[1] - n_dir[1]*delta_t[0] != 0.0)
            t_dir << n_dir[1] * (n_dir[1]*delta_t[0] - n_dir[0]*delta_t[1]),
                     n_dir[0] * (n_dir[0]*delta_t[1] - n_dir[1]*delta_t[0]);
        else
            t_dir << -n_dir[1], n_dir[0];
        t_dir.normalize();

        force_n = (fnc_const * pow(delta_n, 3/2.)
                + fnd_const * pow(delta_n, 1/4.) * v01.dot(n_dir)) * n_dir;
        force_t = (ftc_const * sqrt(delta_n) * delta_t.norm()
                + ftd_const * pow(delta_n, 1/4.) * v01.dot(t_dir)) * t_dir;

        fmax = force_n.norm() * params.mu_s;
        if (force_t.norm() <= fmax)
            force = force_n + force_t;
        else
            force = force_n + fmax * force_t.normalized();
    }
    acceleration = force / params.m;
    acceleration_matrix.col(0) = -acceleration;
    acceleration_matrix.col(1) = acceleration;
    return acceleration_matrix;
}

Eigen::Vector2d angular_acceleration_2(
        const parameters &params, const ball &b1, const ball &b2,
        const Eigen::Vector2d &acceleration)
{
    static const double moment_inertia = 2/5.*params.m * pow(params.rb, 2);
    static Eigen::Vector2d aa, n_dir, r1;
    static double torque_var;
    n_dir = (b2.x - b1.x).normalized();
    r1 = (b2.x - b1.x)/2.0;
    torque_var = -fabs(acceleration.dot(n_dir)) * params.m * params.mu_r * r1.norm();
    aa << cross2d(r1, params.m * acceleration) + torque_var * fsign(b1.w),
          cross2d(r1, params.m * acceleration) + torque_var * fsign(b2.w);
    aa /= moment_inertia;
    return aa;
}

void simulate_2(const parameters &params, std::vector<ball> &b)
{
    const int steps = std::ceil(params.tmax / params.dt);
    int output_interval;

    std::vector<ball> b1, b2, b3, b4;
    Eigen::Matrix2d d2x1, d2x2, d2x3, d2x4, cross_matrix, ac;
    Eigen::Vector2d d2a1, d2a2, d2a3, d2a4, *delta_t = NULL, delta_t_0,
            r01, v01, n_dir, ac_n, ac_t;
    bool colliding, reset_delta_t;
    cross_matrix << 0, -1, 1, 0;

    if (steps > params.output_lines)
        output_interval = steps / params.output_lines;
    else
        output_interval = 1;

    for (int i=0; i<=steps; ++i) {
        // occasional OUTPUT
        if (i % output_interval == 0) {
            std::cout << std::setprecision(14) << params.dt*i;
            for (int j=0; j<2; ++j)
                 std::cout << "\t" << b[j].x[0]
                           << "\t" << b[j].x[1]
                           << "\t" << b[j].v[0]
                           << "\t" << b[j].v[1]
                           << "\t" << b[j].a
                           << "\t" << b[j].w;
            std::cout << std::endl;
        }

        // RK4 INTEGRATION
        b1 = b;
        delta_t_0 << 0, 0;
        colliding = (b1[1].x - b1[0].x).norm() < 2*params.rb;
        if (NULL != delta_t && colliding)
            delta_t_0 = delta_t[0];
        d2x1 = a_total_2(params, b1, delta_t);
        d2a1 << 0, 0;
        if (colliding)
            d2a1 = angular_acceleration_2(params, b1[0], b1[1], d2x1.col(0));

        b2 = b;
        for (int j=0; j<2; ++j) {
            b2[j].x = b1[j].x + b1[j].v*0.5*params.dt;
            b2[j].a = b1[j].a + b1[j].w*0.5*params.dt;
            b2[j].v = b1[j].v + d2x1.col(j)*0.5*params.dt;
            b2[j].w = b1[j].w + d2a1(j)*0.5*params.dt;
        }
        r01 = b2[1].x - b2[0].x;
        colliding = r01.norm() < 2*params.rb;
        if (NULL != delta_t && colliding) {
            n_dir = r01.normalized();
            v01 = b1[1].v - b1[0].v + (b1[0].w + b1[1].w)*params.rb*cross_matrix*n_dir;
            delta_t[0] = delta_t_0 + v01*0.5*params.dt;
        }
        d2x2 = a_total_2(params, b2, delta_t);
        d2a2 << 0, 0;
        if (colliding)
            d2a2 = angular_acceleration_2(params, b2[0], b2[1], d2x2.col(0));

        b3 = b;
        for (int j=0; j<2; ++j) {
            b3[j].x = b1[j].x + b2[j].v*0.5*params.dt;
            b3[j].a = b1[j].a + b2[j].w*0.5*params.dt;
            b3[j].v = b1[j].v + d2x2.col(j)*0.5*params.dt;
            b3[j].w = b1[j].w + d2a2(j)*0.5*params.dt;
        }
        r01 = b3[1].x - b3[0].x;
        colliding = r01.norm() < 2*params.rb;
        if (NULL != delta_t && colliding) {
            n_dir = r01.normalized();
            v01 = b2[1].v - b2[0].v + (b2[0].w + b2[1].w)*params.rb*cross_matrix*n_dir;
            delta_t[0] = delta_t_0 + v01*0.5*params.dt;
        }
        d2x3 = a_total_2(params, b3, delta_t);
        d2a3 << 0, 0;
        if (colliding)
            d2a3 = angular_acceleration_2(params, b3[0], b3[1], d2x3.col(0));

        b4 = b;
        for (int j=0; j<2; ++j) {
            b4[j].x = b1[j].x + b3[j].v*params.dt;
            b4[j].a = b1[j].a + b3[j].w*params.dt;
            b4[j].v = b1[j].v + d2x3.col(j)*params.dt;
            b4[j].w = b1[j].w + d2a3(j)*params.dt;
        }
        r01 = b4[1].x - b4[0].x;
        colliding = r01.norm() < 2*params.rb;
        if (NULL != delta_t && colliding) {
            n_dir = r01.normalized();
            v01 = b3[1].v - b3[0].v + (b3[0].w + b3[1].w)*params.rb*cross_matrix*n_dir;
            delta_t[0] = delta_t_0 + v01*params.dt;
        }
        d2x4 = a_total_2(params, b4, delta_t);
        d2a4 << 0, 0;
        if (colliding)
            d2a4 = angular_acceleration_2(params, b4[0], b4[1], d2x4.col(0));

        for (int j=0; j<2; ++j)
            b[j].x = b1[j].x + params.dt/6.0*(b1[j].v + 2*b2[j].v + 2*b3[j].v + b4[j].v);

        ac = (d2x1 + 2*d2x2 + 2*d2x3 + d2x4)/6.0;

        // checking for COLLISIONS
        reset_delta_t = false;
        r01 = b[1].x - b[0].x;
        colliding = r01.norm() < 2*params.rb;
        if (colliding && NULL == delta_t) {
            // a collision with wall just started
            delta_t = new Eigen::Vector2d;
            delta_t[0] << 0, 0;
        } else if (!colliding && NULL != delta_t) {
            // a collision with wall just ended
            delete delta_t;
            delta_t = NULL;
        } else if (colliding && NULL != delta_t) {
            // a collision continues: tangential force may need truncation
            n_dir = r01.normalized();
            ac_n = ac.col(0).dot(n_dir) * n_dir;
            ac_t = ac.col(0) - ac_n;
            if (ac_t.norm() > params.mu_s*ac_n.norm()) {
                ac.col(0) = ac_n + ac_t.normalized()*params.mu_s*ac_n.norm();
                ac.col(1) = -ac.col(0);
                reset_delta_t = true;
            }
        }

        for (int j=0; j<2; ++j) {
            b[j].a = b1[j].a + params.dt/6.0*(b1[j].w + 2*b2[j].w + 2*b3[j].w + b4[j].w);
            b[j].v = b1[j].v + params.dt * ac.col(j);
            b[j].w = b1[j].w + params.dt/6.0*(d2a1(j) + 2*d2a2(j) + 2*d2a3(j) + d2a4(j));
        }

        if (reset_delta_t)
            delta_t[0] << 0, 0;
        else if (NULL != delta_t && colliding) {
            n_dir = r01.normalized();
            v01 = b[1].v - b[0].v + (b[0].w + b[1].w)*params.rb*cross_matrix*n_dir;
            delta_t[0] = delta_t_0 + v01*params.dt;
        }
    }
}
