#include <ShallowWater.h>
#include <iostream>
#include <cblas.h>


#define g 9.81

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
                    u[i*Nx + j] =  99;
                    v[i*Nx + j] = -99;
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

void ShallowWater::PrintMatrix(const int& N, double* A, const int& lda){ 
    for (int i = 0; i<lda; i++){
        for (int j = 0; j<N; j++) { 
        std::cout << A[j*lda + i] << ",\t";
        }
    std::cout << std::endl;
    }
}


void ShallowWater::EvaluateFuncBLAS(double* uu, double* vv, double* hh, double* f){
    // Declering and defining variables
    int m = 3*Ny*Nx;
    int n = 3*Ny*Nx;
    int kl = 1;
    int ku = 1;
    int lda = 1 + kl + ku;
    
    
    // Declering and defining derivatives array
    double* dhdx = new double[Nx*Ny];
    double* dudx = new double[Nx*Ny];
    double* dvdx = new double[Nx*Ny];
    for (int i = 0; i< Ny; i++){
        
    }
    
    double* dhdy = new double[Nx*Ny];
    double* dudy = new double[Nx*Ny];
    double* dvdy = new double[Nx*Ny];

    

    double* A = new double[3*n];
    double* B = new double[3*n];
    double* C = new double[n];
    double* dXdx = new double[n];
    double* dXdy = new double[n];
    
    // Populate matrices
    for (int i = 0; i< n; i+=3){
        // Construct A
        if (i == 0)
        A[i*lda + ku] = A[(i+1)*lda + ku] =  A[(i+2)*lda + ku] =  uu[i/3];
        A[(i+1)*lda] = hh[i/3];
        A[i*lda + kl + ku] = g;
        // Construct B
        B[i*lda + ku] = B[(i+1)*lda + ku] =  B[(i+2)*lda + ku] =  vv[i/3];
        B[(i+2)*lda] = hh[i/3];
        B[(i+1)*lda + kl + ku] = g;
        // Construct C
        C[i] = uu[i/3]*dhdx[i/3] + vv[i/3]*dhdy[i/3];
//        C[i] = u[i/3]*h[i/3] + v[i/3]*h[i/3];
        // Construct dXdx and dXdy
        dXdx[i] = dhdx[i/3];
        dXdx[i+1] = dudx[i/3];
        dXdx[i+2] = dvdx[i/3];
        dXdy[i] = dhdy[i/3];
        dXdy[i+1] = dudy[i/3];
        dXdy[i+2] = dvdy[i/3];
        
    }
//    for(int i = 0; i< n; i++){
//        std::cout << dXdx[i] << "\t" << "\t" << dXdy[i]<< std::endl;
//    }
//    ShallowWater::PrintMatrix(1, C, n);
    f = new double[n];
    cblas_dgbmv(CblasColMajor, CblasNoTrans, m, n, kl, ku, -1.0, A, lda, dXdx, 1, 0.0,  f, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, m, n, kl, ku, -1.0, B, lda, dXdy, 1, 1.0, f, 1);
    cblas_daxpy( n, 1.0, C, 1, f, 1);
    
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] dXdx;
    delete[] dXdy;
    delete[] dhdx;
    delete[] dhdy;
    delete[] dudx;
    delete[] dudy;
    delete[] dvdx;
    delete[] dvdy;
}

int ShallowWater::PBCI(const int& index, const int& N ){
    if (index < 0){
        return (N-index);
    } else {
        return index;
    }
}

void ShallowWater::TimeIntegrate(){
    
    int m = Ny; 
    int n = Nx;
    int lda = m; // Column major
    
    double* dhdx = new double[Nx*Ny];
    double* dudx = new double[Nx*Ny];
    double* dvdx = new double[Nx*Ny];

    double* dhdy = new double[Nx*Ny];
    double* dudy = new double[Nx*Ny];
    double* dvdy = new double[Nx*Ny];
    
    double coeffs[6] = {-0.016667, 0.15, -0.75, 0.75, -0.15, 0.016667};
    // Serial Impelemntation
    
    for (double t = dt; t < T; t+=dt){
    
        // Calculate derivatives in direction x
        for (int iy = 0; iy < Ny; iy+=){
            for (int ix = 0; ix<Nx; ix+=){
                if 
                double* vect = {u[iy + (ix)*lda
            }
        }
        
        
        
        // Calculate derivatives in direction y
        
        
        
        
    }
    
    delete[] dhdx;
    delete[] dhdy;
    delete[] dudx;
    delete[] dudy;
    delete[] dvdx;
    delete[] dvdy;
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