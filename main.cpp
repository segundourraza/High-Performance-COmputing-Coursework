#include <iostream>
#include <ShallowWater.h>

int main(int argc, char** argv)
{
    std::cout << "You have entered " << argc
         << " arguments:" << "\n";
    
//    double dt = atof(argv[1]);
//    double T = atof(argv[2]);
//    int nx = atoi(argv[3]);
//    int ny = atoi(argv[4]);
//    int ic = atoi(argv[5]);

    // Hard code parameters initially
    double dt = 0.1;
    double T = 25.1;
    int nx = 10;
    int ny = 12;
    int ic = 1;
   
    // Fixed parameters
    double dx = 1.;
    double dy = 1.; 
    
    std::cout << "dt = " << dt << std::endl;
    std::cout << "T = " << T << std::endl;
    std::cout << "Nx = " << nx << std::endl;
    std::cout << "Ny = " << ny << std::endl;
    std::cout << "Initial condition index = " << ic << std::endl;
    
    return 0;
}