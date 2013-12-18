#include <iostream>
#include "parse_cli_opts.h"
#include "ball_mill.h"

using namespace std;

int main(int argc, char *argv[])
{
    bool return_status = 0;
    parameters params;
    return_status += parse_cli_opts(argc, argv, params);

    vector<ball> b;

    if (1 == params.type) {
        type_1_init(params, b);
        cout << "x: " << b[0].x << ", v: " << b[0].v << endl;
    } else if (2 == params.type) {
        // type_2_init(params, b);
    } else if (3 == params.type) {
        // type_3_init(params, b);
    }

    return return_status;
}
