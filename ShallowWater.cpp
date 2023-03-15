#include <ShallowWater.h>
#include <iostream>
#include <fstream> 
#include <cblas.h>
#include <cstdlib>

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
    
    for (int i = 0; i<Nx*Ny; i++){
        u[i] =  0;
        v[i] =  0;
    }
    // Populate 2 dimensional array depending on Index of initial condition
    
    switch (ic){
        case 1:
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
//                    g[i][j] = (double) (std::exp(-(i*dx-50)*(i*dx-50)/25));
//                    g[i][j] = i+1 + (j+1)*10;
                    h[i*Nx + j] = (double) (10 + std::exp(-(i*dx-Nx/2.)*(i*dx-Nx/2.)/(Nx/4.)));
                }
            }
            break;
        case 2:
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
//                    h[i*Nx + j] = (double) ( std::exp(-(j*dy-50)*(j*dy-50)/25));
                    h[i*Nx + j] = (double) (10+ std::exp(-(j*dy-Nx/2)*(j*dy-Nx/2)/(Nx/4)));
                }
            }
            break;
        case 3: 
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    h[i*Nx + j] = (double) (10 + std::exp(-((i*dx-50)*(i*dx-50) + (j*dy-50)*(j*dy-50))/25));
                }
            }
            break;
        case 4: 
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    h[i*Nx + j] = (double) (10 + std::exp(-((i*dx-25)*(i*dx-25) + (j*dy-25)*(j*dy-25))/25) + std::exp(-((i*dx-75)*(i*dx-75) + (j*dy-75)*(j*dy-75))/25));
                }
            }
            break;
    }
}

void ShallowWater::PrintVector(const int& N, const double* x){
    for (int i = 0; i<N; i++){
        std::cout << x[i] << std::endl;
    }
}

void ShallowWater::PrintMatrix(const int& N, const double* A, const int& lday, const int& inc){ 
    for (int i = 0; i<lday; i+=inc){
        for (int j = 0; j<N; j++) {
            if(std::abs(A[j*lday + (i)]) < 1e-13){
                std::cout << 0 << ",\t";
            }
            else{
                std::cout << A[j*lday + (i)] << ",\t";
            }
        }
    std::cout << "\n" << std::endl;
    }
}

//void ShallowWater::TimeIntegrate(){
//    
//    int m = Ny; 
//    int n = Nx;
//    int lday = m; // Column major
//    
//    double* dhdx = new double[Nx*Ny];
//    double* dudx = new double[Nx*Ny];
//    double* dvdx = new double[Nx*Ny];
//
//    double* dhdy = new double[Nx*Ny];
//    double* dudy = new double[Nx*Ny];
//    double* dvdy = new double[Nx*Ny];
//
//    double coeffs[6] = {-0.016667, 0.15, -0.75, 0.75, -0.15, 0.016667};
//    double RK4[4] = {dt/6, dt/3, dt/3, dt/6};
//    double Kcoeffs[4] = {0, dt/2, dt/2, dt};
//    
//    // Serial Impelemntation
//
////    for (double t = dt; t < T; t+=dt){
//    for (int t = 0; t<1; t++){
//        
//        for (int k = 0; k < 1; k++){
//            
//            // Calculate derivatives in direction x and y (ASSUME SQUARE)
//            
//            // Boundary points for ix <3 and ix > Nx-3
//            for (int iy = 0; iy < Ny; iy++){
//                
//                // x derivatives
//                
//                // Left most points
//                dhdx[iy] = coeffs[0]*h[iy + (Nx-3)*lday] + coeffs[1]*h[iy + (Nx-2)*lday] + coeffs[2]*h[iy + (Nx-1)*lday] + coeffs[3]*h[iy + (1)*lday] + coeffs[4]*h[iy + (2)*lday] + coeffs[5]*h[iy + (3)*lday];
//                dhdx[iy + 1*lday] = coeffs[0]*h[iy + (Nx-2)*lday] + coeffs[1]*h[iy + (Nx-1)*lday] + coeffs[2]*h[iy + (0)*lday] + coeffs[3]*h[iy + (2)*lday] + coeffs[4]*h[iy + (3)*lday] + coeffs[5]*h[iy + (4)*lday];
//                dhdx[iy + 2*lday] = coeffs[0]*h[iy +  (Nx-1)*lday] + coeffs[1]*h[iy + (0)*lday] + coeffs[2]*h[iy + (1)*lday] + coeffs[3]*h[iy + (3)*lday] + coeffs[4]*h[iy + (4)*lday] + coeffs[5]*h[iy + (5)*lday];
//                
//                dudx[iy] = coeffs[0]*u[iy + (Nx-3)*lday] + coeffs[1]*u[iy + (Nx-2)*lday] + coeffs[2]*u[iy + (Nx-1)*lday] + coeffs[3]*u[iy + (1)*lday] + coeffs[4]*u[iy + (2)*lday] + coeffs[5]*u[iy + (3)*lday];
//                dudx[iy + 1*lday] = coeffs[0]*u[iy + (Nx-2)*lday] + coeffs[1]*u[iy + (Nx-1)*lday] + coeffs[2]*u[iy + (0)*lday] + coeffs[3]*u[iy + (2)*lday] + coeffs[4]*u[iy + (3)*lday] + coeffs[5]*u[iy + (4)*lday];
//                dudx[iy + 2*lday] = coeffs[0]*u[iy + (Nx-1)*lday] + coeffs[1]*u[iy + (0)*lday] + coeffs[2]*u[iy + (1)*lday] + coeffs[3]*u[iy + (3)*lday] + coeffs[4]*u[iy + (4)*lday] + coeffs[5]*u[iy + (5)*lday];
//                
//                dvdx[iy] = coeffs[0]*v[iy + (Nx-3)*lday] + coeffs[1]*v[iy + (Nx-2)*lday] + coeffs[2]*v[iy + (Nx-1)*lday] + coeffs[3]*v[iy + (1)*lday] + coeffs[4]*v[iy + (2)*lday] + coeffs[5]*v[iy + (3)*lday];
//                dvdx[iy + 1*lday] = coeffs[0]*v[iy + (Nx-2)*lday] + coeffs[1]*v[iy + (Nx-1)*lday] + coeffs[2]*v[iy + (0)*lday] + coeffs[3]*v[iy + (2)*lday] + coeffs[4]*v[iy + (3)*lday] + coeffs[5]*v[iy + (4)*lday];
//                dvdx[iy + 2*lday] = coeffs[0]*v[iy + (Nx-1)*lday] + coeffs[1]*v[iy + (0)*lday] + coeffs[2]*v[iy + (1)*lday] + coeffs[3]*v[iy + (3)*lday] + coeffs[4]*v[iy + (4)*lday] + coeffs[5]*v[iy + (5)*lday];
//                
//                // Right most points
//                dhdx[iy+(Nx-1)*lday] = coeffs[0]*h[iy + (Nx-4)*lday] + coeffs[1]*h[iy + (Nx-3)*lday] + coeffs[2]*h[iy + (Nx-2)*lday] + coeffs[3]*h[iy + (0)*lday] + coeffs[4]*h[iy + (1)*lday] + coeffs[5]*h[iy + (2)*lday];
//                dhdx[iy+(Nx-2)*lday] = coeffs[0]*h[iy + (Nx-5)*lday] + coeffs[1]*h[iy + (Nx-4)*lday] + coeffs[2]*h[iy + (Nx-3)*lday] + coeffs[3]*h[iy + (Nx-1)*lday] + coeffs[4]*h[iy + (0)*lday] + coeffs[5]*h[iy + (1)*lday];
//                dhdx[iy+(Nx-3)*lday] = coeffs[0]*h[iy + (Nx-6)*lday] + coeffs[1]*h[iy + (Nx-5)*lday] + coeffs[2]*h[iy + (Nx-4)*lday] + coeffs[3]*h[iy + (Nx-2)*lday] + coeffs[4]*h[iy + (Nx-1)*lday] + coeffs[5]*h[iy + (0)*lday];
//                
//                dudx[iy+(Nx-1)*lday] = coeffs[0]*u[iy + (Nx-4)*lday] + coeffs[1]*u[iy + (Nx-3)*lday] + coeffs[2]*u[iy + (Nx-2)*lday] + coeffs[3]*u[iy + (0)*lday] + coeffs[4]*u[iy + (1)*lday] + coeffs[5]*u[iy + (2)*lday];
//                dudx[iy+(Nx-2)*lday] = coeffs[0]*u[iy + (Nx-5)*lday] + coeffs[1]*u[iy + (Nx-4)*lday] + coeffs[2]*u[iy + (Nx-3)*lday] + coeffs[3]*u[iy + (Nx-1)*lday] + coeffs[4]*u[iy + (0)*lday] + coeffs[5]*u[iy + (1)*lday];
//                dudx[iy+(Nx-3)*lday] = coeffs[0]*u[iy + (Nx-6)*lday] + coeffs[1]*u[iy + (Nx-5)*lday] + coeffs[2]*u[iy + (Nx-4)*lday] + coeffs[3]*u[iy + (Nx-2)*lday] + coeffs[4]*u[iy + (Nx-1)*lday] + coeffs[5]*u[iy + (0)*lday];
//                
//                dudx[iy+(Nx-1)*lday] = coeffs[0]*v[iy + (Nx-4)*lday] + coeffs[1]*v[iy + (Nx-3)*lday] + coeffs[2]*v[iy + (Nx-2)*lday] + coeffs[3]*v[iy + (0)*lday] + coeffs[4]*v[iy + (1)*lday] + coeffs[5]*v[iy + (2)*lday];
//                dudx[iy+(Nx-2)*lday] = coeffs[0]*v[iy + (Nx-5)*lday] + coeffs[1]*v[iy + (Nx-4)*lday] + coeffs[2]*v[iy + (Nx-3)*lday] + coeffs[3]*v[iy + (Nx-1)*lday] + coeffs[4]*v[iy + (0)*lday] + coeffs[5]*v[iy + (1)*lday];
//                dudx[iy+(Nx-3)*lday] = coeffs[0]*v[iy + (Nx-6)*lday] + coeffs[1]*v[iy + (Nx-5)*lday] + coeffs[2]*v[iy + (Nx-4)*lday] + coeffs[3]*v[iy + (Nx-2)*lday] + coeffs[4]*v[iy + (Nx-1)*lday] + coeffs[5]*v[iy + (0)*lday];
//                
//                // y derivatives
//                
//            }
//            
//            // Inner points
//            
//            for (int iy = 0; iy < Ny; iy++){
//                for (int ix = 3; ix<Nx-3; ix++){
//                    dhdx[iy+lday*ix] = coeffs[0]*h[iy + (ix-3)*lday] + coeffs[1]*h[iy + (ix-2)*lday] + coeffs[2]*h[iy + (ix-1)*lday] + coeffs[3]*h[iy + (ix+1)*lday] + coeffs[4]*h[iy + (ix+2)*lday] + coeffs[5]*h[iy + (ix+3)*lday];
//                    
//                    dudx[iy+lday*ix] = coeffs[0]*u[iy + (ix-3)*lday] + coeffs[1]*u[iy + (ix-2)*lday] + coeffs[2]*u[iy + (ix-1)*lday] + coeffs[3]*u[iy + (ix+1)*lday] + coeffs[4]*u[iy + (ix+2)*lday] + coeffs[5]*u[iy + (ix+3)*lday];
//                    
//                    dvdx[iy+lday*ix] = coeffs[0]*v[iy + (ix-3)*lday] + coeffs[1]*v[iy + (ix-2)*lday] + coeffs[2]*v[iy + (ix-1)*lday] + coeffs[3]*v[iy + (ix+1)*lday] + coeffs[4]*v[iy + (ix+2)*lday] + coeffs[5]*v[iy + (ix+3)*lday];
//                }
//            }
//            
//            std::cout << std::endl;
//            ShallowWater::PrintMatrix(n, dhdx, lday);
//            
//             
//            // Calculate derivatives in direction y
//            
//            
//            
//            
//        }
//    }
//    delete[] dhdx;
//    delete[] dhdy;
//    delete[] dudx;
//    delete[] dudy;
//    delete[] dvdx;
//    delete[] dvdy;
//}



void ShallowWater::TimeIntegrate(){ 
    
    int ldsy = 3*Ny;
    int dimS = ldsy*Nx;
    int kl = 3; 
    int ku = 3;
    int lday = 1+ kl + ku;
    
    // Initialize variables
    double* A = new double[lday*Nx];
    double* S = new double[dimS];
    double* Snew = new double[dimS];
    double* k2 = new double[dimS];
    double* k1 = new double[dimS];
    
    double coeffs[6] = {-0.016667, 0.15, -0.75, 0.75, -0.15, 0.016667};
    double RK4coeffs[4] = {dt/6, dt/3, dt/3, dt/6};
    double kcoeffs[3] = {dt/2, dt/2, dt};
    
    // Construct state vector S = [u, v, h]^T. [u,v,h] for all nodes stack on yop of each other in a column major way
    ShallowWater::ConstructSVector(S);
    
    // Construct derivative matrix A
    ShallowWater::PopulateA(Nx, A, lday, coeffs);
    
    // Start integration loop 
    double t = 0;
    while (t < T){
        t += dt; 
        
        // Calculate k1 and propagate Snew
        EvaluateFuncBlasV2(kl, ku, A, lday, S, ldsy, coeffs, k1);
        for (int i = 0; i<dimS; i++){
            Snew[i] = S[i] + RK4coeffs[0]*k1[i];
            S[i] += kcoeffs[0]*k1[i];
        }
        
        // Calculate k2 and propagate Snew
        EvaluateFuncBlasV2(kl, ku, A, lday, S, ldsy, coeffs, k2);
        for (int i = 0; i<dimS; i++){
            Snew[i] += RK4coeffs[1]*k2[i];
            S[i] += kcoeffs[1]*k2[i] -kcoeffs[0]*k1[i];
        }
        
        // Calculate k3 and propagate Snew
        EvaluateFuncBlasV2(kl, ku, A, lday, S, ldsy, coeffs, k1);
        for (int i = 0; i<dimS; i++){
            Snew[i] += RK4coeffs[2]*k1[i];
            S[i] +=  kcoeffs[2]*k1[i]- kcoeffs[1]*k2[i];
        }
        
        // Calculate k4 and update S for next iteration
        EvaluateFuncBlasV2(kl, ku, A, lday, S, ldsy, coeffs, k2);
        for (int i = 0; i<dimS; i++){
            S[i] = Snew[i] + RK4coeffs[3]*k2[i];
        }
    
        std::cout << "\rTime step: " << t<< " seconds"<< std::endl;
    }  
    std::cout << "Writting file!" << std::endl;
    WriteFile(S);
}

void ShallowWater::EvaluateFuncBlasV2(const int& kla, const int& kua, const double* A, const int& lday, double* S, const int& ldsy, const double* coeffs, double* k){
    // k = F(S) = - B*d(S)/dx - C*d(S)/dy
    int dimS = ldsy*Nx;
    
    double* dSdx = new double[dimS];
    double* dSdy = new double[dimS];
    
    // Step 1: Evaluate derivatives of State S
    getDerivatives(kla, kua, A, lday, S, ldsy, dSdx, dSdy, coeffs);
    
    // Step 2: Construct banded matrix B
    int kl = 2;
    int ku = 2;
    int ldy = 1 + kl + ku;
    double* B = new double[ldy*dimS];
    
    for (int i = 0; i < dimS; i+=3){ 
        if (i != 0) {B[(i-1)*ldy] = g;}
        B[i*ldy+ku] = B[(i+1)*ldy+ku] = B[(i+2)*ldy+ku] = S[i];
        if (i!=dimS-1){B[(i+1)*ldy-1] = S[i+2];}
    }
    B[(dimS-1)*ldy] = g;
     
    // Step 3: Evaluate b*d(S)/dx
    cblas_dgbmv(CblasColMajor, CblasNoTrans, dimS, dimS, kl, ku, -1.0, B, ldy, dSdx, 1, 0, k, 1);
        
    // Step 4: Construct banded matrix C
    kl = 1;
    ku = 1;
    ldy = 1 + kl + ku;
    double* C = new double[ldy*dimS];
    
    for (int i = 0; i < dimS; i+=3){
        C[i*ldy+ku] = C[(i+1)*ldy+ku] = C[(i+2)*ldy+ku] = S[i+1];
        if (i != 0) {C[(i-1)*ldy] = g;}
        if (i<dimS-1){C[(i+2)*ldy-1] = S[i+2];}
    }
    
    // Step 5: Evaluate value of function f(S)
    cblas_dgbmv(CblasColMajor, CblasNoTrans, dimS, dimS, kl, ku, -1.0, C, ldy, dSdy, 1, 1.0, k, 1);
    
    delete[] B;
    delete[] C;
    delete[] dSdx;
    delete[] dSdy;
}


void ShallowWater::EvaluateFuncBlas(const int& dimS, double* S, const double* dSdx, const double* dSdy, double* k){
    // k = F(S) = - B*d(S)/dx - C*d(S)/dy
    // Construct banded matrix B
    int kl = 2;
    int ku = 2;
    int ldy = 1 + kl + ku;
    double* B = new double[ldy*dimS];
    
    for (int i = 0; i < dimS; i+=3){ 
        if (i != 0) {B[(i-1)*ldy] = g;}
        B[i*ldy+ku] = B[(i+1)*ldy+ku] = B[(i+2)*ldy+ku] = S[i];
        if (i!=dimS-1){B[(i+1)*ldy-1] = S[i+2];}
    }
     B[(dimS-1)*ldy] = g;
     
    cblas_dgbmv(CblasColMajor, CblasNoTrans, dimS, dimS, kl, ku, -1.0, B, ldy, dSdx, 1, 0, k, 1);
          
    
//        PrintMatrix(10, B+(dimS-10)*ldy, ldy, 1);
//        std::cout << std::endl;
//        std::cout << std::endl;
//        std::cout << std::endl;
//        std::cout << std::endl;      
//          
          
    for (int i = 0; i< ldy*dimS; i++){
        B[i] = 0;
    }
  
    // Construct banded matrix C, overwrite B
    kl = 1;
    ku = 1;
    ldy = 1 + kl + ku;
    delete[] B;
    
    
    B = new double[ldy*dimS];
    for (int i = 0; i < dimS; i+=3){
        B[i*ldy+ku] = B[(i+1)*ldy+ku] = B[(i+2)*ldy+ku] = S[i+1];
        if (i != 0) {B[(i-1)*ldy] = g;}
        if (i<dimS-1){B[(i+2)*ldy-1] = S[i+2];}
        
    }
     
    cblas_dgbmv(CblasColMajor, CblasNoTrans, dimS, dimS, kl, ku, -1.0, B, ldy, dSdy, 1, 1.0, k, 1);
        
    for (int i = 0; i< ldy*dimS; i++){
        B[i] = 0;
    }
    delete[] B;
}


void ShallowWater::getDerivatives(const int& kl, const int& ku, const double* A, const int& lday, const double* S, const int& ldsy, double* dSdx, double* dSdy, const double* coeffs){
    
    // Iterate over rows
    for (int iy = 0; iy < ldsy; iy ++){
        cblas_dgbmv(CblasColMajor, CblasNoTrans, Ny, Nx, kl, ku, 1/dx, A, lday, S+iy, ldsy, 0.0, dSdx + iy, ldsy);
    }
    
    // Iterate over columns
    for (int ix = 0 ; ix < Nx; ix ++ ){
        cblas_dgbmv(CblasColMajor, CblasNoTrans, Ny, Nx, kl, ku, 1/dy, A, lday, S + ix*ldsy, 3, 0.0, dSdy+ix*ldsy, 3);
        cblas_dgbmv(CblasColMajor, CblasNoTrans, Ny, Nx, kl, ku, 1/dy, A, lday, S + 1 + ix*ldsy, 3, 0.0, dSdy + 1 + ix*ldsy, 3);
        cblas_dgbmv(CblasColMajor, CblasNoTrans, Ny, Nx, kl, ku, 1/dy, A, lday, S + 2 + ix*ldsy, 3, 0.0, dSdy + 2 + ix*ldsy, 3);
    }
    
    ShallowWater::ApplyPeriodicBC(Nx, S, ldsy, dSdx, dSdy, coeffs);    
}



void ShallowWater::ApplyPeriodicBC(const int& Nx, const double* S, const int& ldsy, double* dSdx, double* dSdy, const double* coeffs){    
    // Apply BC on x derivatives
    // Loop over rows
    for (int iy = 0; iy < ldsy; iy++){
        dSdx[iy] += coeffs[2]*S[iy+(Nx-1)*ldsy] + coeffs[1]*S[iy+(Nx-2)*ldsy] + coeffs[0]*S[iy+(Nx-3)*ldsy];
        dSdx[iy+ldsy] += coeffs[1]*S[iy+(Nx-1)*ldsy] + coeffs[0]*S[iy+(Nx-2)*ldsy];
        dSdx[iy+2*ldsy] += coeffs[0]*S[iy+(Nx-1)*ldsy];
        
        dSdx[iy + (Nx-1)*ldsy] += coeffs[3]*S[iy] + coeffs[4]*S[iy+ldsy] + coeffs[5]*S[iy + 2*ldsy];
        dSdx[iy + (Nx-2)*ldsy] += coeffs[4]*S[iy] + coeffs[5]*S[iy+ldsy];
        dSdx[iy + (Nx-3)*ldsy] += coeffs[5]*S[iy];
        
    }    
    
    // Apply BC on y derivatives
    // Loop over columns
    for (int ix = 0; ix< Nx; ix++){
        for (int j = 0; j<3; j++){
        dSdy[ix*ldsy+j] += coeffs[2]*S[(ix+1)*ldsy-3+j] + coeffs[1]*S[(ix+1)*ldsy-6+j] + coeffs[0]*S[(ix+1)*ldsy-9+j];
        dSdy[ix*ldsy+3+j] += coeffs[1]*S[(ix+1)*ldsy-3+j] + coeffs[0]*S[(ix+1)*ldsy-6+j];
        dSdy[ix*ldsy+6+j] += coeffs[0]*S[(ix+1)*ldsy-3+j];
        
        dSdy[(ix+1)*ldsy -3+j] += coeffs[3]*S[ix*ldsy+j] + coeffs[4]*S[ix*ldsy+3+j] + coeffs[5]*S[ix*ldsy +6+j];
        dSdy[(ix+1)*ldsy -6+j] += coeffs[4]*S[ix*ldsy+j] + coeffs[5]*S[ix*ldsy+3+j];
        dSdy[(ix+1)*ldsy -9+j] += coeffs[5]*S[ix*ldsy+j];
        }
    }
}


void ShallowWater::ConstructSVector(double* S){
    for (int i = 0; i<3*Nx*Ny; i+=3){
        S[i] =  u[i/3];
        S[i+1] = v[i/3];
        S[i+2] = h[i/3];
    }
    
}
void ShallowWater::PopulateA(const int& N, double* A, const int& lday, const double* coeffs){
    // Banded
    if (Nx == Ny){
        for (int i =0; i < N; i++){
            for (int j = 0; j <3; j++){
                A[i*lday + j] = coeffs[5-j];
                A[(i+1)*lday -1 -j] = coeffs[j];
            }
        }
    }
}

void ShallowWater::WriteFile(const double* S){
    int ldsy = 3*Ny;
    int incy = 3;
    
    std::ofstream myfile;
    myfile.open("Output.txt");
    for (int iy = 0; iy< Ny; iy++){
        for (int ix = 0; ix<Nx; ix++){
            myfile << ix*dx << "\t" << iy*dy << "\t" << S[iy*incy+ix*ldsy] << "\t" << S[iy*incy+ix*ldsy + 1] << "\t" << S[iy*incy+ix*ldsy + 2] << "\n"; 
        }
    }
    myfile << "Writing output to a file .\n";
    myfile.close();
    
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