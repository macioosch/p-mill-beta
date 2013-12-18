#ifndef BALL_MILL_H
#define BALL_MILL_H

#include <boost/fusion/container/vector.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "parse_cli_opts.h"

using namespace std;
namespace blas = boost::numeric::ublas;
namespace fu = boost::fusion;

struct ball {
    blas::vector<double> x;
    blas::vector<double> v;
    double a, w;
};

int type_1_init(parameters &params, vector<ball> &b);

#endif // BALL_MILL_H
