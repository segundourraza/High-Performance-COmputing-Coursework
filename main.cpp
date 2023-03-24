#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
//#include <boost/timer/timer.hpp>
#include <chrono>

#include "ShallowWater.h"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
    // Boost program options
    po::options_description opts("Allowed options");
    opts.add_options()
        ("help", "produce help message")
        ("dt", po::value<double>()->default_value(0.1), "Time-step to use.")
        ("T", po::value<double>()->default_value(80), "Total integration time.")
        ("Nx", po::value<int>()->default_value(100), "Number of grid points in x.")
        ("Ny", po::value<int>()->default_value(100), "Number of grid points in y.")
        ("ic", po::value<double>()->default_value(1), "Index of the initial condition to use (1-4).")
        ("mode", po::value<int>()->default_value(2), "Analysis mode. [1] - BLAS analysis, [2] - for based analysis");
        
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
    const int ic        = vm["ic"].as<double>();
    int analysis        = vm["mode"].as<int>();; // [1] - BLAS analysis, [2] - for based analysis
    
    // Fixed parameters
    double dx = 1.;
    double dy = 1.; 
    
    // Testing class ShallowWater
    ShallowWater sol1(dt, T, Nx, Ny, ic, dx, dy, analysis);
    std::cout << "\nSIMUALTION PARAMETERS:" << std::endl;
    std::cout << "\t" << "Time-step:\t" << "\t" <<  "\t" <<  sol1.getTimeStep() << std::endl;
    std::cout << "\t" << "Total integration time:\t" << "\t" << sol1.getIntegrationTime() << std::endl;
    std::cout << "\t" << "Number of grid points in x:\t" << sol1.getNx() << std::endl;
    std::cout << "\t" << "Number of grid points in y:\t" << sol1.getNy() << std::endl;
    std::cout << "\t" << "Spatial step in x:\t" << "\t" <<  dx << std::endl;
    std::cout << "\t" << "Spatial step in y:\t" << "\t" <<  dy << std::endl;
    std::cout << "\t" << "Initial condition index:\t" << sol1.getIc() << std::endl;
    
    sol1.SetInitialCondition(); 
    
    if (analysis == 1){
        std::cout << "\t" << "Implemenatation mode:\t\tBLAS\n" << std::endl;
        sol1.TimeIntegrateBLAS();    
    }
    else if (analysis == 2){
        std::cout << "\t" << "Implemenatation mode:\t\t" << "FOR LOOP" << std::endl;
//        sol1.TimeIntegrateForLoop();
        sol1.TimeIntegrate();
    }
    
    sol1.WriteFile();
    
    int y1 = 88;
    int x1 = 26;
    int y2 = 26;
    int x2 = 21;
    
    std::cout << "\nSIMUALTION RESULTS:" << std::endl;
    std::cout << "\t" << std::setprecision (16) << std::fixed << "h[" << y1 << "," << x1 << "] = " << "\t" <<*(sol1.geth() + y1 +Ny*(x1)) << std::endl;
    std::cout << "\t" << std::setprecision (16) << std::fixed << "h[" << y2 << "," << x2 << "] = " << "\t" <<*(sol1.geth() + y2 +Ny*x2) << std::endl;
    
    return 0;
}