#include "ShallowWater.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream> 
#include <cblas.h>
#include <cstdlib>

#include <omp.h>

#define g 9.81

// constructor definition
ShallowWater::ShallowWater(){
}   // Default Constructor

ShallowWater::ShallowWater(double dtt, double Tt, int Nxx, int Nyy, int icc, double dxx, double dyy, int typeAnalysis) : dt(dtt), T(Tt), Nx(Nxx), Ny(Nyy), ic(icc), dx(dxx), dy(dyy), analysis(typeAnalysis){
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
                    h[i*Nx + j] = (double) (10 + std::exp(-((i*dx-50)*(i*dx-50) + (j*dy-50)*(j*dy-50))/25.));
                }
            }
            break;
        case 4: 
            for (int i = 0; i<Nx; i++){
                for (int j = 0; j<Ny; j++){
                    h[i*Nx + j] = (double) (10 + std::exp(-((i*dx-25)*(i*dx-25) + (j*dy-25)*(j*dy-25))/25.) + std::exp(-((i*dx-75)*(i*dx-75) + (j*dy-75)*(j*dy-75))/25.));
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

void ShallowWater::TimeIntegrateForLoop(){
    std::cout << std::setprecision(16) << std::fixed;
    std::string str;
    
    int threadid;
    int NumThreads;
    
    int nodes = Nx*Ny;
    int local_n;
    int reminder =0;
    
    double* kutemp = new double[Nx*Ny];
    double* kvtemp = new double[Nx*Ny];
    double* khtemp = new double[Nx*Ny];    
    
    double* unew = new double[Nx*Ny];
    double* vnew = new double[Nx*Ny];
    double* hnew = new double[Nx*Ny];
    
    double* dhdx = new double[Nx*Ny];
    double* dudx = new double[Nx*Ny];
    double* dvdx = new double[Nx*Ny];

    double* dhdy = new double[Nx*Ny];
    double* dudy = new double[Nx*Ny];
    double* dvdy = new double[Nx*Ny];

    double* ku = new double[Nx*Ny];
    double* kv = new double[Nx*Ny];
    double* kh = new double[Nx*Ny];
    
    double coeffs[6] = {-0.016667, 0.15, -0.75, 0.75, -0.15, 0.016667};
    
    // Open branch of threads
    #pragma omp parallel default(shared) private(threadid) 
    {
        threadid = omp_get_thread_num();
        if (threadid == 0){
            NumThreads = omp_get_num_threads();
            local_n = nodes/NumThreads;
            reminder = nodes - local_n*NumThreads;
            std::cout << "Number of threads: " << NumThreads << "." << std::endl;
            std::cout << "Number of nodes: " << nodes << "." << std::endl;
            std::cout << "Number of nodes per thread: " << local_n << "." << std::endl;
            std::cout << "Number of reminding nodes: " << reminder << "." << std::endl;
            
        }
        #pragma omp barrier
        
        // Start integration loop 
        double t = dt;
        while (t < T + dt/2){
            
            // Calculate k1 and propagate Snew
            GetDerivativesForLoop(u, dudx, dudy, coeffs);
            GetDerivativesForLoop(v, dvdx, dvdy, coeffs);
            GetDerivativesForLoop(h, dhdx, dhdy, coeffs);
            
            for (int node = 0; node < Nx*Ny; node++){
                ku[node] = -u[node]*dudx[node] - v[node]*dudy[node] - g*dhdx[node];
                unew[node] = u[node] + dt/6 * ku[node];
                
                kv[node] =  -u[node]*dvdx[node] - v[node]*dvdy[node] - g*dhdy[node];
                vnew[node] = v[node] + dt/6 * kv[node];
                
                kh[node] = -h[node]*dudx[node] - u[node]*dhdx[node] - h[node]*dvdy[node] - v[node]*dhdy[node];
                hnew[node] = h[node] + dt/6 * kh[node];
                
                u[node] += dt/2*ku[node];
                v[node] += dt/2*kv[node];
                h[node] += dt/2*kh[node];
            }
            
            // Calculate k2 and propagate Snew
            GetDerivativesForLoop(u, dudx, dudy, coeffs);
            GetDerivativesForLoop(v, dvdx, dvdy, coeffs);
            GetDerivativesForLoop(h, dhdx, dhdy, coeffs);
            
            for (int node = 0; node < Nx*Ny; node++){
                kutemp[node] = ku[node];
                ku[node] = -u[node]*dudx[node] - v[node]*dudy[node] - g*dhdx[node];
                unew[node] += dt/3 * ku[node];
                
                kvtemp[node] = kv[node];
                kv[node] =  -u[node]*dvdx[node] - v[node]*dvdy[node] - g*dhdy[node];
                vnew[node] += dt/3 * kv[node];
                
                khtemp[node] = kh[node];
                kh[node] = -h[node]*dudx[node] - u[node]*dhdx[node] - h[node]*dvdy[node] - v[node]*dhdy[node];
                hnew[node] += dt/3 * kh[node];
                
                u[node] += dt/2*ku[node] - dt/2*kutemp[node];
                v[node] += dt/2*kv[node] - dt/2*kvtemp[node];
                h[node] += dt/2*kh[node] - dt/2*khtemp[node];
            }
            
            // Calculate k3 and propagate Snew
            GetDerivativesForLoop(u, dudx, dudy, coeffs);
            GetDerivativesForLoop(v, dvdx, dvdy, coeffs);
            GetDerivativesForLoop(h, dhdx, dhdy, coeffs);
            
            for (int node = 0; node < Nx*Ny; node++){
                kutemp[node] = ku[node];
                ku[node] = -u[node]*dudx[node] - v[node]*dudy[node] - g*dhdx[node];
                unew[node] += dt/3 * ku[node];
                
                kvtemp[node] = kv[node];
                kv[node] =  -u[node]*dvdx[node] - v[node]*dvdy[node] - g*dhdy[node];
                vnew[node] += dt/3 * kv[node];
                
                khtemp[node] = kh[node];
                kh[node] = -h[node]*dudx[node] - u[node]*dhdx[node] - h[node]*dvdy[node] - v[node]*dhdy[node];
                hnew[node] +=dt/3 * kh[node];
                
                u[node] += dt*ku[node] - dt/2*kutemp[node];
                v[node] += dt*kv[node] - dt/2*kvtemp[node];
                h[node] += dt*kh[node] - dt/2*khtemp[node];
            }
            
            // Calculate k4 and propagate Snew
            GetDerivativesForLoop(u, dudx, dudy, coeffs);
            GetDerivativesForLoop(v, dvdx, dvdy, coeffs);
            GetDerivativesForLoop(h, dhdx, dhdy, coeffs);
            
            for (int node = 0; node < Nx*Ny; node++){
                ku[node] = -u[node]*dudx[node] - v[node]*dudy[node] - g*dhdx[node];
                kv[node] =  -u[node]*dvdx[node] - v[node]*dvdy[node] - g*dhdy[node];
                kh[node] = -h[node]*dudx[node] - u[node]*dhdx[node] - h[node]*dvdy[node] - v[node]*dhdy[node];
                
                u[node] = unew[node] + dt/6 * ku[node];
                v[node] = vnew[node] + dt/6 * kv[node];
                h[node] = hnew[node] + dt/6 * kh[node];
            }
            std::cout << std::string(str.length(),'\b');
            str = "Time: " + std::to_string(t) + ". " + std::to_string((int) ((t)/dt)) + " time steps done out of " + std::to_string((int) (T/dt)) + ".";
            std::cout << str;
            
            t += dt; 
        }
    }
    
    delete[] kutemp;
    delete[] kvtemp;
    delete[] khtemp;
    
    delete[] unew;
    delete[] vnew;
    delete[] hnew;
    
    delete[] dhdx;
    delete[] dudx;
    delete[] dvdx;

    delete[] dhdy;
    delete[] dudy;
    delete[] dvdy;

    delete[] ku;
    delete[] kv;
    delete[] kh;
    
}

void ShallowWater::GetDerivativesForLoop(const double* var, double* dvardx, double* dvardy, const double* coeffs){
    // Calculate derivatives in direction x and y (ASSUME SQUARE)
    int ldy = Ny;    
    
    // X - DERIVATIVES
    for (int iy = 0; iy < Ny; iy++){
        // Boundary points for ix <3 and ix > Nx-3
        
        // Left most points
        dvardx[iy] = coeffs[0]*var[iy + (Nx-4)*ldy] + coeffs[1]*var[iy + (Nx-3)*ldy] + coeffs[2]*var[iy + (Nx-2)*ldy] + coeffs[3]*var[iy + (1)*ldy] + coeffs[4]*var[iy + (2)*ldy] + coeffs[5]*var[iy + (3)*ldy];
        dvardx[iy + 1*ldy] = coeffs[0]*var[iy + (Nx-3)*ldy] + coeffs[1]*var[iy + (Nx-2)*ldy] + coeffs[2]*var[iy + (0)*ldy] + coeffs[3]*var[iy + (2)*ldy] + coeffs[4]*var[iy + (3)*ldy] + coeffs[5]*var[iy + (4)*ldy];
        dvardx[iy + 2*ldy] = coeffs[0]*var[iy +  (Nx-2)*ldy] + coeffs[1]*var[iy + (0)*ldy] + coeffs[2]*var[iy + (1)*ldy] + coeffs[3]*var[iy + (3)*ldy] + coeffs[4]*var[iy + (4)*ldy] + coeffs[5]*var[iy + (5)*ldy];
        
        // Right most points
        dvardx[iy+(Nx-1)*ldy] = coeffs[0]*var[iy + (Nx-4)*ldy] + coeffs[1]*var[iy + (Nx-3)*ldy] + coeffs[2]*var[iy + (Nx-2)*ldy] + coeffs[3]*var[iy + (1)*ldy] + coeffs[4]*var[iy + (2)*ldy] + coeffs[5]*var[iy + (3)*ldy];
        dvardx[iy+(Nx-2)*ldy] = coeffs[0]*var[iy + (Nx-5)*ldy] + coeffs[1]*var[iy + (Nx-4)*ldy] + coeffs[2]*var[iy + (Nx-3)*ldy] + coeffs[3]*var[iy + (Nx-1)*ldy] + coeffs[4]*var[iy + (1)*ldy] + coeffs[5]*var[iy + (2)*ldy];
        dvardx[iy+(Nx-3)*ldy] = coeffs[0]*var[iy + (Nx-6)*ldy] + coeffs[1]*var[iy + (Nx-5)*ldy] + coeffs[2]*var[iy + (Nx-4)*ldy] + coeffs[3]*var[iy + (Nx-2)*ldy] + coeffs[4]*var[iy + (Nx-1)*ldy] + coeffs[5]*var[iy + (1)*ldy];
    
        // Inner points
        for (int ix = 3; ix<Nx-3; ix++){
            dvardx[iy+ldy*ix] = coeffs[0]*var[iy + (ix-3)*ldy] + coeffs[1]*var[iy + (ix-2)*ldy] + coeffs[2]*var[iy + (ix-1)*ldy] + coeffs[3]*var[iy + (ix+1)*ldy] + coeffs[4]*var[iy + (ix+2)*ldy] + coeffs[5]*var[iy + (ix+3)*ldy];
        }
    }
        
    // Y - DERIVATVES
    
    for (int ix = 0; ix < Nx; ix++){
        // Boundary points for iy <3 and iy > Nx-3
        
        // Top points
        dvardy[0+ix*ldy] = coeffs[0]*var[ix*ldy + Ny-4] + coeffs[1]*var[ix*ldy + Ny -3] + coeffs[2]*var[ix*ldy + Ny-2] + coeffs[3]*var[ix*ldy+1] + coeffs[4]*var[ix*ldy+2] + coeffs[5]*var[ix*ldy+3];
        dvardy[1+ix*ldy] = coeffs[0]*var[ix*ldy + Ny-3] + coeffs[1]*var[ix*ldy + Ny -2] + coeffs[2]*var[ix*ldy + 0] + coeffs[3]*var[ix*ldy+2] + coeffs[4]*var[ix*ldy+3] + coeffs[5]*var[ix*ldy+4];
        dvardy[2+ix*ldy] = coeffs[0]*var[ix*ldy + Ny-2] + coeffs[1]*var[ix*ldy + 0] + coeffs[2]*var[ix*ldy + 1] + coeffs[3]*var[ix*ldy+3] + coeffs[4]*var[ix*ldy+4] + coeffs[5]*var[ix*ldy+5];
        
        
        // Bottom points
        dvardy[Ny-1+ix*ldy] = coeffs[0]*var[ix*ldy + Ny-4] + coeffs[1]*var[ix*ldy + Ny -3] + coeffs[2]*var[ix*ldy + Ny-2] + coeffs[3]*var[ix*ldy+1] + coeffs[4]*var[ix*ldy +2] + coeffs[5]*var[ix*ldy +3];
        dvardy[Ny-2+ix*ldy] = coeffs[0]*var[ix*ldy + Ny-5] + coeffs[1]*var[ix*ldy + Ny -4] + coeffs[2]*var[ix*ldy + Ny-3] + coeffs[3]*var[ix*ldy+Ny-1] + coeffs[4]*var[ix*ldy +1] + coeffs[5]*var[ix*ldy +2];
        dvardy[Ny-3+ix*ldy] = coeffs[0]*var[ix*ldy + Ny-6] + coeffs[1]*var[ix*ldy + Ny -5] + coeffs[2]*var[ix*ldy + Ny-4] + coeffs[3]*var[ix*ldy+Ny-2] + coeffs[4]*var[ix*ldy + Ny-1] + coeffs[5]*var[ix*ldy + 1];
        
        // Inner points
        for (int iy = 3; iy<Ny-3; iy++){
            dvardy[iy+ldy*ix] = coeffs[0]*var[iy-3 + ix*ldy] + coeffs[1]*var[iy -2 + ix*ldy] + coeffs[2]*var[iy -1 + ix*ldy] + coeffs[3]*var[iy + 1 + ix*ldy] + coeffs[4]*var[iy + 2 + ix*ldy] + coeffs[5]*var[iy + 3 + ix*ldy];
        }
    }
}

void ShallowWater::TimeIntegrate(){ 
    
    std::string str;
    
    B = new double[5*3*Ny*Nx];
    C = new double[3*3*Nx*Ny];
    
    // Populate Differentiation matrix (Only Required by BLAS implementation)
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
    double t = dt;
    while (t < T + dt/2){
        
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
        
        std::cout << std::string(str.length(),'\b');
        str = "Time: " + std::to_string(t) + ". " + std::to_string((int) (t/dt)) + " time steps done out of " + std::to_string((int) (T/dt)) + ".\n";
        std::cout << str;
        t += dt;
    }  
    
    for (int i = 0; i<dimS; i++){
        u[i/3] = S[i];
        v[i/3] = S[i+1];
        h[i/3] = S[i+2];
    }
    
    delete[] B;
    delete[] C;
}

void ShallowWater::EvaluateFuncBlasV2(const int& kla, const int& kua, const double* A, const int& lday, double* S, const int& ldsy, const double* coeffs, double* k){
    // k = F(S) = - B*d(S)/dx - C*d(S)/dy
    int dimS = ldsy*Nx;
    
    double* dSdx = new double[dimS];
    double* dSdy = new double[dimS];
    
    // Step 1: Evaluate derivatives of State S
    GetDerivativesBlas(kla, kua, A, lday, S, ldsy, dSdx, dSdy, coeffs);
    
    // Step 2: Construct banded matrix B
    int kl = 2;
    int ku = 2;
    int ldy = 1 + kl + ku;
    
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
    
    for (int i = 0; i < dimS; i+=3){
        C[i*ldy+ku] = C[(i+1)*ldy+ku] = C[(i+2)*ldy+ku] = S[i+1];
        if (i != 0) {C[(i-1)*ldy] = g;}
        if (i<dimS-1){C[(i+2)*ldy-1] = S[i+2];}
    }
    
    // Step 5: Evaluate value of function f(S)
    cblas_dgbmv(CblasColMajor, CblasNoTrans, dimS, dimS, kl, ku, -1.0, C, ldy, dSdy, 1, 1.0, k, 1);
    
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

void ShallowWater::GetDerivativesBlas(const int& kl, const int& ku, const double* A, const int& lday, const double* S, const int& ldsy, double* dSdx, double* dSdy, const double* coeffs){
    
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
        dSdx[iy] += coeffs[2]*S[iy+(Nx-2)*ldsy] + coeffs[1]*S[iy+(Nx-3)*ldsy] + coeffs[0]*S[iy+(Nx-4)*ldsy];
        dSdx[iy+ldsy] += coeffs[1]*S[iy+(Nx-2)*ldsy] + coeffs[0]*S[iy+(Nx-3)*ldsy];
        dSdx[iy+2*ldsy] += coeffs[0]*S[iy+(Nx-2)*ldsy];
        
        dSdx[iy + (Nx-1)*ldsy] += coeffs[3]*S[iy+1*ldsy] + coeffs[4]*S[iy+2*ldsy] + coeffs[5]*S[iy + 3*ldsy];
        dSdx[iy + (Nx-2)*ldsy] += coeffs[4]*S[iy+1*ldsy] + coeffs[5]*S[iy+2*ldsy];
        dSdx[iy + (Nx-3)*ldsy] += coeffs[5]*S[iy+1*ldsy];
        
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

void ShallowWater::WriteFile(){
    std::ofstream myfile;
    myfile.open("Output.txt");
    for (int iy = 0; iy< Ny; iy++){
        for (int ix = 0; ix<Nx; ix++){
            myfile << ix*dx << "\t" << iy*dy << "\t" << u[iy+ix*Ny] << "\t" << v[iy+ix*Ny] << "\t" << h[iy+ix*Ny] << "\n"; 
        }
    }
    std::cout << "\nWriting output to file.\n" << std::endl;
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