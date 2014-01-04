#ifndef BALL_MILL_H
#define BALL_MILL_H

#include <eigen3/Eigen/Core>
#include <vector>
#include "parse_cli_opts.h"

struct ball {
    Eigen::Vector2d x, v;
    double a, w;
    double *wall_delta_t;
};

void type_1_init(parameters &params, std::vector<ball> &b);
void simulate_1(const parameters &params, std::vector<ball> &b);

#endif // BALL_MILL_H
