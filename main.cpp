#include <iostream>
#include "parse_cli_opts.h"

using namespace std;

int main(int argc, char *argv[])
{
    bool return_status = 0;
    parameters params;
    return_status += parse_cli_opts(argc, argv, params);

    cout << "Accessing params from main loop: " << params.type << endl;

    return return_status;
}
