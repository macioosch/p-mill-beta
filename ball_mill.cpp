#include "ball_mill.h"

int type_1_init(parameters &params, vector<ball> &b)
{
    b.resize(1);

    b[0].a = 0.0;
    b[0].w = 0.0;
    b[0].x.resize(2);
    b[0].x(0) = 0.0;
    b[0].x(1) = params.y0;
    b[0].v.resize(2);
    b[0].v(0) = params.vx0;
    b[0].v(1) = params.vy0;

    return 0;
}
