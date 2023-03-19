#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include "ShallowWater.h"

namespace po = boost::program_options;

int main(int argc, char* argv[])
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


    
//    const double dt     = 0.1;
//    const double T      = 5;
//    const int Nx        = 11;
//    const int Ny        = 11;
//    const int ic         = 1;
    int analysis = 2; // [1] - BLAS analysis, [2] - for based analysis
    // Fixed parameters
    double dx = 1.;
    double dy = 1.; 
    
    // Testing class ShallowWater
    ShallowWater sol1(dt, T, Nx, Ny, ic, dx, dy, analysis);
//    ShallowWater sol1;
    sol1.sayHello();
    
    std::cout << "\nTime-step: " << sol1.getTimeStep() << std::endl;
    std::cout << "Total integration time: " << sol1.getIntegrationTime() << std::endl;
    std::cout << "Number of grid points in x = " << sol1.getNx() << std::endl;
    std::cout << "Number of grid points in y = " << sol1.getNy() << std::endl;
    std::cout << "Initial condition index = " << sol1.getIc() << std::endl;
    std::cout << std::endl;
    
    sol1.SetInitialCondition(); 
    
    if (analysis == 1){
        std::cout << "BLAS IMPLEMENTATION SELECTED" << std::endl;
        sol1.TimeIntegrate();    
    }
    else if (analysis == 2){
        std::cout << "FOR LOOP IMPLEMENTATION SELECTED" << std::endl;
        sol1.TimeIntegrateForLoop();
    }
    
    sol1.WriteFile();
    
    
    int y1 = 88;
    int x1 = 26;
    int y2 = 26;
    int x2 = 21;
    std::cout << std::setprecision (16) << std::fixed << "h[" << y1 << "," << x1 << "] = " << *(sol1.geth() + y1 +Ny*(x1)) << std::endl;
    std::cout << std::setprecision (16) << std::fixed << "h[" << y2 << "," << x2 << "] = " << *(sol1.geth() + y2 +Ny*x2) << std::endl;
    return 0;
}