#ifndef BALL_MILL_H
#define BALL_MILL_H

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "parse_cli_opts.h"

namespace blas = boost::numeric::ublas;

struct ball {
    blas::vector<double> x, v;
    double a, w;
};

void type_1_init(parameters &params, std::vector<ball> &b);
void simulate(parameters &params, std::vector<ball> &b);

#endif // BALL_MILL_H
