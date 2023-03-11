#include <ShallowWater.h>
#include <iostream>
#include <cblas.h>

// constructor definition
ShallowWater::ShallowWater(){
}   // Default Constructor

ShallowWater::ShallowWater(double dtt, double Tt, int Nxx, int Nyy, int icc, double dxx, double dyy) : dt(dtt), T(Tt), Nx(Nxx), Ny(Nyy), ic(icc), dx(dxx), dy(dyy){
}   // Constructor using initialization list to avoid calling default class constructor and then over writting


ShallowWater::~ShallowWater(){
std::cout << "Class destroyed" << std::endl;

delete[] h;
delete[] u;
delete[] v;
}   // Custom destructor definition

// Method definition
void ShallowWater::sayHello(){
std::cout << "Hello from class 'ShallowWater'!" << std::endl;
}

void ShallowWater::SetInitialCondition(){// ARRAY OF POINTER TO POINTER WILL BE GENERATED. 
    // Output[i] will contain the ith column of the initial condition, corresponding
    // to the points [x0,y0], [x0,y1], ... [x0,yn].

    h = new double[Nx*Ny];
    u = new double[Nx*Ny];
    v = new double[Nx*Ny];
    // Populate 2 dimensional array depending on Index of initial condition
    switch (ic){
        case 1:
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
//                    g[i][j] = (double) (std::exp(-(i*dx-50)*(i*dx-50)/25));
//                    g[i][j] = i+1 + (j+1)*10;
                    h[i*Nx + j] = (double) (std::exp(-(i*dx-Nx/2)*(i*dx-Nx/2)/(Nx/4)));
                    u[i*Nx + j] = 0;
                    v[i*Nx + j] = 0;
                }
            }
            break;
        case 2:
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
//                    h[i*Nx + j] = (double) ( std::exp(-(j*dy-50)*(j*dy-50)/25));
                    h[i*Nx + j] = (double) (std::exp(-(j*dy-Nx/2)*(j*dy-Nx/2)/(Nx/4)));
                    u[i*Nx + j] = 0;
                    v[i*Nx + j] = 0;
                }
            }
            break;
        case 3: 
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    h[i*Nx + j] = (double) (std::exp(-((i*dx-50)*(i*dx-50) + (j*dy-50)*(j*dy-50))/25));
                    u[i*Nx + j] = 0;
                    v[i*Nx + j] = 0;
                }
            }
            break;
        case 4: 
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    h[i*Nx + j] = (double) (std::exp(-((i*dx-25)*(i*dx-25) + (j*dy-25)*(j*dy-25))/25) + std::exp(-((i*dx-75)*(i*dx-75) + (j*dy-75)*(j*dy-75))/25));
                    u[i*Nx + j] = 0;
                    v[i*Nx + j] = 0;
                }
            }
            break;
    }
}

void ShallowWater::PrintMatrix(double* A){ 
    for (int i = 0; i<Ny; i++){
        for (int j = 0; j<Nx; j++) { 
        std::cout << A[j*Nx + i] << ",\t";
        }
    std::cout << std::endl;
    }
}


void ShallowWater::GetDerivatives(const char& dir, double* f, double* df){
    // Performs differentiation of array f in direction 'dir'. If dir == 'x', column  wise
    // differentiation is performed. If dir == 'y', row wise differentiation is performed.
    // On output, f is rewritten with the derivatives    
    double coeffs[6] = {-0.016667, 0.15, -0.75, 0.75, -0.15, 0.016667};
    if (dir == 'x'){
        double step = dx;
        for (int iy = 0 ; iy < Ny; iy++){
            for (int ix = 0; ix < Nx; ix++){
                if (ix == 0){
                    double vect[6] = {f[iy + (Nx-3)*Ny] , f[iy + (Nx-2)*Ny], f[iy + (Nx-1)*Ny], f[iy + (1)*Ny], f[iy + (2)*Ny], f[iy + (3)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                }else if (ix == 1){
                    double vect[6] = {f[iy + (Nx-2)*Ny] , f[iy + (Nx-1)*Ny], f[iy + 0*Ny], f[iy + (2)*Ny], f[iy + (3)*Ny], f[iy + (4)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                } else if (ix == 2){
                    double vect[6] = {f[iy + (Nx-1)*Ny] , f[iy + (0)*Ny], f[iy + (1)*Ny], f[iy + (3)*Ny], f[iy + (4)*Ny], f[iy + (5)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                }else if (ix == Nx-3){
                    double vect[6] = {f[iy + (ix-3)*Ny] , f[iy + (ix-2)*Ny], f[iy + (ix-1)*Ny], f[iy + (ix+1)*Ny], f[iy + (ix+2)*Ny], f[iy + (0)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                } else if (ix == Nx-2){
                    double vect[6] = {f[iy + (ix-3)*Ny] , f[iy + (ix-2)*Ny], f[iy + (ix-1)*Ny], f[iy + (ix+1)*Ny], f[iy + (0)*Ny], f[iy + (1)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;   
                } else if (ix == Nx-1){
                    double vect[6] = {f[iy + (ix-3)*Ny] , f[iy + (ix-2)*Ny], f[iy + (ix-1)*Ny], f[iy + (0)*Ny], f[iy + (1)*Ny], f[iy + (2)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                }else {
                    df[iy +ix*Ny] = cblas_ddot(6, coeffs, 1, &f[iy + (ix-3)*Ny], Ny)/step;
                }
            }
        }
    }else if (dir =='y'){
        double step = dy;
        for (int ix = 0 ; ix < Nx; ix++){
            for (int iy = 0; iy < Ny; iy++){
                if (iy == 0){
                    double vect[6] = {f[Ny-3 + ix*Ny] , f[Ny-2 + ix*Ny], f[Ny-1 + ix*Ny], f[iy+1 + (ix)*Ny], f[iy+2 + (ix)*Ny], f[iy+3 + (ix)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                }else if (iy == 1){
                    double vect[6] = {f[Ny-2 + ix*Ny] , f[Ny-1 + ix*Ny], f[iy-1 + ix*Ny], f[iy+1 + (ix)*Ny], f[iy+2 + (ix)*Ny], f[iy+3 + (ix)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                } else if (iy == 2){
                    double vect[6] = {f[Ny-1 + ix*Ny] , f[iy-2 + ix*Ny], f[iy-1 + ix*Ny], f[iy+1 + (ix)*Ny], f[iy+2 + (ix)*Ny], f[iy+3 + (ix)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                }else if (iy == Ny-3){
                    double vect[6] = {f[iy -3 + ix*Ny] , f[iy -2 + ix*Ny], f[iy -1 + ix*Ny], f[iy+1 + (ix)*Ny], f[iy+2 + (ix)*Ny], f[0 + (ix)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                } else if (iy == Ny-2){
                    double vect[6] = {f[iy -3 + ix*Ny] , f[iy -2 + ix*Ny], f[iy -1 + ix*Ny], f[iy+1 + (ix)*Ny], f[0 + (ix)*Ny], f[1 + (ix)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                } else if (iy == Ny-1){
                    double vect[6] = {f[iy -3 + ix*Ny] , f[iy -2 + ix*Ny], f[iy -1 + ix*Ny], f[0 + (ix)*Ny], f[1+ (ix)*Ny], f[2+ (ix)*Ny]};
                    df[iy + ix*Ny] = cblas_ddot(6, coeffs, 1, vect, 1)/step;
                }else {
                    df[iy +ix*Ny] = cblas_ddot(6, coeffs, 1, &f[iy-3 + ix*Ny], 1)/step;
                }
            }
        }
    }else {
        throw std::invalid_argument("ERROR: Invalid differentiation direction. Direction must be 'x' or 'y'");
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
double* ShallowWater::geth(){return h;}
double* ShallowWater::getu(){return u;}
double* ShallowWater::getv(){return v;}

// Structure 'limits'

 Limits::Limits(double xa = 0, double xb = 0, double ya = 0, double yb = 0): xl{xa}, xu{xb}, yl{ya}, yu{yb} {
    if(xa>xb) std::swap(xl, xu);
    if(ya>yb) std::swap(yl, yu);
}   // Default Constructor of structure 'Limit'
    
Limits::~Limits(){
    std::cout << "Structure destroyed" << std::endl;
}   // Default destructor