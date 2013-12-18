#ifndef BALL_MILL_H
#define BALL_MILL_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace blas = boost::numeric::ublas;

struct ball {
    blas::vector<double> x;
    blas::vector<double> v;
    double a, w;
};

#endif // BALL_MILL_H
