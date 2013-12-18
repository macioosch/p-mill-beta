/* A placeholder.
 * TODO:
 * - function should return a struct of meaningful options.
 */
#include <iostream>
#include <boost/program_options.hpp>
#include "parse_cli_opts.h"

using namespace std;
namespace po = boost::program_options;

int parse_cli_opts(int argc, char *argv[], parameters &params) {

    po::options_description opt_general("General options");
    opt_general.add_options()
        ("help", "produce this help message")
        ("type", po::value<int>(&(params.type)), "simulation type")
    ;

    po::options_description opt_balls("Balls' options");
    opt_balls.add_options()
        ("e", po::value<double>(&(params.e)), "balls' restitution coefficient")
        ("E", po::value<double>(&(params.E)), "balls' Young's modulus [Pa]")
        ("G", po::value<double>(&(params.G)), "balls' shear modulus [Pa]")
        ("mu", po::value<double>(&(params.mu)),
            "balls' rolling friction coefficient")
        ("nu", po::value<double>(&(params.nu)), "balls' Poisson's ratio")
        ("rho", po::value<double>(&(params.rho)), "balls' density [kg m^-3]")
        ("rb", po::value<double>(&(params.rb)), "balls' radius [m]")
    ;

    po::options_description opt_type1(
                "Options for type 1 simulation, a ball over a plane");
    opt_type1.add_options()
        ("g", po::value<double>(&(params.g)),
            "gravitational acceleration [m s^-2]")
        ("vx0", po::value<double>(&(params.vx0)),
            "ball's initial horizontal velocity")
        ("vy0", po::value<double>(&(params.vy0)),
            "ball's initial vertical velocity")
        ("w0", po::value<double>(&(params.w0)),
            "ball's initial angular velocity")
        ("y0", po::value<double>(&(params.y0)),
            "ball's initial distance from the plane")
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
        cout << "Type was set to " << params.type << ".\n";
    } else {
        cerr << "Please set simulation type!\n";
        return 1;
    }

    return 0;
}
