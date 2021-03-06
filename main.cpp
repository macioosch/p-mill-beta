#include <iostream>
#include <vector>
#include "ball_mill.h"
#include "parse_cli_opts.h"

int main(int argc, char *argv[])
{
    bool return_status = 0;
    parameters params;
    return_status += parse_cli_opts(argc, argv, params);

    std::vector<ball> b;

    if (1 == params.type) {
        type_1_init(params, b);
        simulate_1(params, b);
    } else if (2 == params.type) {
        type_2_init(params, b);
        simulate_2(params, b);
    } else if (3 == params.type) {
        type_3_init(params, b);
        simulate_3(params, b);
    } else if (4 == params.type) {
        // (sic!) - type 2 init, type 3 sim
        type_2_init(params, b);
        simulate_3(params, b);
    }

    return return_status;
}
