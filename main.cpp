#include <iostream>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include "ShallowWater.h"

namespace po = boost::program_options;

int main(int argc, char** argv)
{
    // Boost program options
    po::options_description opts("Allowed options");
    opts.add_options()
        ("help", "produce help message")
        ("dt", po::value<double>()->default_value(0.1), "Time-step to use.")
        ("T", po::value<double>()->default_value(5.0), "Total integration time.")
        ("Nx", po::value<int>()->default_value(100), "Number of grid points in x.")
        ("Ny", po::value<int>()->default_value(100), "Number of grid points in y.")
        ("ic", po::value<double>()->default_value(1), "Index of the initial condition to use (1-4).");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << opts << "\n";
        return 1;
    }

    // Loading parameters into memory
    const double dt     = vm["dt"].as<double>();
    const double T      = vm["T"].as<double>();
    const int Nx        = vm["Nx"].as<int>();
    const int Ny        = vm["Ny"].as<int>();
    const int ic         = vm["ic"].as<double>();

    int analysis = 2; // [1] - BLAS analysis, [2] - for based analysis
    // Fixed parameters
    double dx = 1.;
    double dy = 1.; 
    
    // Testing class ShallowWater
    ShallowWater sol1(dt, T, Nx, Ny, ic, dx, dy, analysis);
//    ShallowWater sol1;
    sol1.sayHello();
    
    std::cout << "\ndt = " << sol1.getTimeStep() << std::endl;
    std::cout << "T = " << sol1.getIntegrationTime() << std::endl;
    std::cout << "Nx = " << sol1.getNx() << std::endl;
    std::cout << "Ny = " << sol1.getNy() << std::endl;
    std::cout << "Initial condition index = " << sol1.getIc() << std::endl;
    std::cout << "dx = " << sol1.getdx() << std::endl;
    std::cout << "dy = " << sol1.getdy() << std::endl;
    std::cout << std::endl;
    
    sol1.SetInitialCondition(); 

    sol1.TimeIntegrateForLoop();

    return 0;
}