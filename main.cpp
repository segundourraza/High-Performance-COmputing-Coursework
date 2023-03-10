#include <iostream>
#include <ShallowWater.h>

int main(int argc, char** argv)
{
    std::cout << "You have entered " << argc
         << " arguments:" << "\n";
    
    double dt = atof(argv[1]);
    double T = atof(argv[2]);
    int nx = atoi(argv[3]);
    int ny = atoi(argv[4]);
    int ic = atoi(argv[5]);
    
    std::cout << "dt = " << dt << std::endl;
    std::cout << "T = " << T << std::endl;
    std::cout << "Nx = " << nx << std::endl;
    std::cout << "Ny = " << ny << std::endl;
    std::cout << "Initial condition index = " << ic << std::endl;
    
    return 0;
}