#include <iostream>
#include <cmath>

struct Limits{
    double xl, xu, yl, yu;
public:
    Limits(double xa, double xb, double ya, double yb); // Constructor declaration
        
    ~Limits(); // Default destructor declaration
};

#ifndef SHALLOWWATER_H
#define SHALLOWWATER_H

class ShallowWater
{
    // Default initialisation
    double dt = 0.1;
    double T = 25.1;
    int Nx = 10;
    int Ny = 12;
    int ic = 1;
    double dx = 1.;
    double dy = 1.; 
    
public:
    // Constructors
    ShallowWater(); // Default Constructorn declaration
    ShallowWater(double dtt, double Tt, int Nxx, int Nyy, int icc, double dxx, double dyy); // Constructor declaration
    
    // Methods
    double** SetInitialCondition(){// ARRAY OF POINTER TO POINTER WILL BE GENERATED. 
    // Output[i] will contain the ith column of the initial condition, corresponding
    // to the points [x0,y0], [x0,y1], ... [x0,yn].
    
//    int Nx = ShallowWater.Nx;
//    int Ny = ShallowWater.Ny;
//    double dx = ShallowWater.dx;
//    double dy = ShallowWater.dy;
//    int ic = ShallowWater.ic;
//    
    // Create array of pointer containing data from each column dynamically
    double** g = new double*[Nx];
    
    // Dynamically associate each column pointer to an array
    for (int i = 0; i<Nx; i++){
        g[i] = new double[Ny];
    }
    
    // Populate 2 dimensional array depending on Index of initial condition
    switch (ic){
        case 1:
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    g[i][j] = (double) (std::exp(-(i*dx-50)*(i*dx-50)/25));
                }
            }
            break;
        case 2:
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    g[i][j] = (double) ( std::exp(-(j*dy-50)*(j*dy-50)/25));
                }
            }
            break;
        case 3: 
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    g[i][j] = (double) (std::exp(-((i*dx-50)*(i*dx-50) + (j*dy-50)*(j*dy-50))/25));
                }
            }
            break;
        case 4: 
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    g[i][j] = (double) (std::exp(-((i*dx-25)*(i*dx-25) + (j*dy-25)*(j*dy-25))/25) + std::exp(-((i*dx-75)*(i*dx-75) + (j*dy-75)*(j*dy-75))/25));
                }
            }
            break;
    }
    // Print 2 dimensional array
    for (int i = 0; i<Ny; i++){
        for (int j = 0; j<Nx; j++) { 
            std::cout << g[j][i] << ",\t";
        }
        std::cout << std::endl;
    }
    
    return g;
}

    void sayHello();
    
    // 'Getter' functions
    double getTimeStep();
    double getIntegrationTime();
    int getNx();
    int getNy();
    int getIc();
    double getdx();
    double getdy();
    
    ~ShallowWater(); // Destructor
    
};

#endif