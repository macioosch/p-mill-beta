#ifndef BALL_MILL_H
#define BALL_MILL_H

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "parse_cli_opts.h"

namespace blas = boost::numeric::ublas;

struct collision_with_wall {
    // collision in "point", r0 = point - b.x, a0 = b.a
    blas::vector<double> r0, point;
    double a0;
};

struct ball {
    blas::vector<double> x, v;
    double a, w;
    collision_with_wall* wall_collision;
};

void type_1_init(parameters &params, std::vector<ball> &b);
void simulate_1(const parameters &params, std::vector<ball> &b);

#endif // BALL_MILL_H
