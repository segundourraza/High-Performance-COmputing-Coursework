#include <iostream>
#include <ShallowWater.h>
#include <iomanip>

int main(int argc, char** argv)
{
    std::cout << std::setprecision(5) << std::fixed;
    std::cout << "You have entered " << argc
         << " arguments:" << "\n";
    
//    double dt = atof(argv[1]);
//    double T = atof(argv[2]);
//    int nx = atoi(argv[3]);
//    int ny = atoi(argv[4]);
//    int ic = atoi(argv[5]);

    // Hard code parameters initially
    double dt = 0.1;
    double T = 5;
    int nx = 100;
    int ny = 100;
    int ic = 3;
   
    // Fixed parameters
    double dx = 1.;
    double dy = 1.; 
    
    // Testing class ShallowWater
    ShallowWater sol1(dt, T, nx, ny, ic, dx, dy);
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
//    sol1.PrintMatrix(sol1.getNx(), sol1.geth(), sol1.getNy());

    std::cout << std::endl;
    sol1.TimeIntegrate();
    
    // Testing Limits structure
    Limits lim1(0., dx*nx, 0., dy*ny);
    std::cout << "\nx lower = " << lim1.xl << std::endl;
    std::cout << "x upper = " << lim1.xu << std::endl;
    std::cout << "y lower = " << lim1.yl << std::endl;
    std::cout << "y upper = " << lim1.yu << std::endl;
    

    return 0;
}