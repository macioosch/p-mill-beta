#include <iostream>
#include <math.h>
#include <boost/program_options.hpp>
#include "parse_cli_opts.h"

using namespace std;
namespace po = boost::program_options;

int parse_cli_opts(int argc, char *argv[], parameters &params) {

    po::options_description opt_general("General options");
    opt_general.add_options()
        ("dt", po::value<double>(&(params.dt)), "time step [s]")
        ("help", "produce this help message")
        ("output_lines", po::value<double>(&(params.output_lines)),
            "number of output lines, every tmax/output_lines s")
        ("tmax", po::value<double>(&(params.tmax)), "total simulation time [s]")
        ("type", po::value<int>(&(params.type)), "simulation type")
    ;

    po::options_description opt_balls("Balls' options");
    opt_balls.add_options()
        ("e", po::value<double>(&(params.e)), "balls' restitution coefficient")
        ("E", po::value<double>(&(params.E)), "balls' Young's modulus [Pa]")
        ("G", po::value<double>(&(params.G)), "balls' shear modulus [Pa]")
        ("mu_r", po::value<double>(&(params.mu_r)), "balls' rolling friction coefficient")
        ("mu_s", po::value<double>(&(params.mu_s)), "balls' static friction coefficient")
        ("rho", po::value<double>(&(params.rho)), "balls' density [kg m^-3]")
        ("rb", po::value<double>(&(params.rb)), "balls' radius [m]")
    ;

    po::options_description opt_type1(
                "Options for type 1 simulation, a ball over a plane");
    opt_type1.add_options()
        ("g", po::value<double>(&(params.g)),
            "gravitational acceleration [m s^-2]")
        ("vx0", po::value<double>(&(params.vx0)),
            "ball's initial horizontal velocity [m s^-1]")
        ("vy0", po::value<double>(&(params.vy0)),
            "ball's initial vertical velocity [m s^-1]")
        ("w0", po::value<double>(&(params.w0)),
            "ball's initial angular velocity [s^-1]")
        ("y0", po::value<double>(&(params.y0)),
            "ball's initial distance from the plane [m]")
    ;

    po::options_description opt_all("Allowed options");
    opt_all.add(opt_general).add(opt_balls).add(opt_type1);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opt_all), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << opt_all << "\n";
        return 1;
    }

    if (vm.count("type")) {
        cout << "# Simulation type set to " << params.type << ".\n";
    } else {
        cerr << "Please set simulation type!\n";
        return 1;
    }

    // secondary parameters
    params.nu = 0.5 * params.E / params.G - 1.0;
    params.Er = params.E / (2.0*(1.0-pow(params.nu,2)));
    params.Gr = params.G / (2.0*(2.0-params.nu));
    params.m = params.rho * (4.0/3.0) * M_PI * pow(params.rb, 3);
    params.mr = params.m / 2.0;
    params.rr = params.rb / 2.0;
    params.beta = -pow(1.0 + pow(M_PI/log(params.e), 2), -0.5);

    return 0;
}
