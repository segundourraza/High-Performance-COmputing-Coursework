#include <ShallowWater.h>
#include <iostream>

// constructor definition
ShallowWater::ShallowWater(){
}   // Default Constructor

ShallowWater::ShallowWater(double dtt, double Tt, int Nxx, int Nyy, int icc, double dxx, double dyy) : dt(dtt), T(Tt), Nx(Nxx), Ny(Nyy), ic(icc), dx(dxx), dy(dyy){
}   // Constructor using initialization list to avoid calling default class constructor and then over writting


ShallowWater::~ShallowWater(){
std::cout << "Class destroyed" << std::endl;

delete[] g;
}   // Custom destructor definition

// Method definition
void ShallowWater::sayHello(){
std::cout << "Hello from class 'ShallowWater'!" << std::endl;
}

void ShallowWater::SetInitialCondition(){// ARRAY OF POINTER TO POINTER WILL BE GENERATED. 
    // Output[i] will contain the ith column of the initial condition, corresponding
    // to the points [x0,y0], [x0,y1], ... [x0,yn].

    g = new double[Nx*Ny];
    
    // Populate 2 dimensional array depending on Index of initial condition
    switch (ic){
        case 1:
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
//                    g[i][j] = (double) (std::exp(-(i*dx-50)*(i*dx-50)/25));
//                    g[i][j] = i+1 + (j+1)*10;
                    g[i*Nx + j] = (double) (std::exp(-(i*dx-Nx/2)*(i*dx-Nx/2)/(Nx/4)));
                }
            }
            break;
        case 2:
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    g[i*Nx + j] = (double) ( std::exp(-(j*dy-50)*(j*dy-50)/25));
                }
            }
            break;
        case 3: 
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    g[i*Nx + j] = (double) (std::exp(-((i*dx-50)*(i*dx-50) + (j*dy-50)*(j*dy-50))/25));
                }
            }
            break;
        case 4: 
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    g[i*Nx + j] = (double) (std::exp(-((i*dx-25)*(i*dx-25) + (j*dy-25)*(j*dy-25))/25) + std::exp(-((i*dx-75)*(i*dx-75) + (j*dy-75)*(j*dy-75))/25));
                }
            }
            break;
    }
}

void ShallowWater::PrintMatrix(){ 
    for (int i = 0; i<Ny; i++){
        for (int j = 0; j<Nx; j++) { 
        std::cout << g[j*Nx + i] << ",\t";
        }
    std::cout << std::endl;
    }
}


// 'Getter' function definition
double ShallowWater::getTimeStep(){ return dt;}
double ShallowWater::getIntegrationTime(){ return T;}
int ShallowWater::getNx(){return Nx;}
int ShallowWater::getNy(){return Ny;}
int ShallowWater::getIc(){return ic;}
double ShallowWater::getdx(){ return dx;}
double ShallowWater::getdy(){return dy;}
double* ShallowWater::getg(){return g;}


// Structure 'limits'

 Limits::Limits(double xa = 0, double xb = 0, double ya = 0, double yb = 0): xl{xa}, xu{xb}, yl{ya}, yu{yb} {
    if(xa>xb) std::swap(xl, xu);
    if(ya>yb) std::swap(yl, yu);
}   // Default Constructor of structure 'Limit'
    
Limits::~Limits(){
    std::cout << "Structure destroyed" << std::endl;
}   // Default destructor