#include <iostream>
#include "parse_cli_opts.h"
#include "ball_mill.h"

using namespace std;

int main(int argc, char *argv[])
{
    bool return_status = 0;
    parameters params;
    return_status += parse_cli_opts(argc, argv, params);

    ball b;

    b.a = 0.0;
    b.w = 0.0;
    b.x.resize(2);
    b.x(0) = 0.0;
    b.x(1) = params.y0;
    b.v.resize(2);
    b.v(0) = params.vx0;
    b.v(1) = params.vy0;

    return return_status;
}
