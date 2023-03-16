#include <iostream>
#include <cmath>

#ifndef SHALLOWWATER_H
#define SHALLOWWATER_H

class ShallowWater
{
    // Default initialisation
    double dt = 0.1;
    double T = 25.1;
    int Nx = 11;
    int Ny = 11;
    int ic = 1;
    double dx = 1.;
    double dy = 1.; 
    int analysis = 1;
    \
    double* h = nullptr;
    double* u = nullptr;
    double* v = nullptr;
    
    void PopulateA(const int& N, double* A, const int& lda, const double* coeffs);
    void ConstructSVector(double* S);
    void GetDerivativesBlas(const int& kl, const int& ku, const double* A, const int& lda, const double* S, const int& ldx, double* dSdx, double* dXdy, const double* coeffs);
    void GetDerivativesForLoop(const double* var, double* dvardx, double* dvardy, const double* coeffs);
    void EvaluateFuncBlasV2(const int& kl, const int& ku, const double* A, const int& lday, double* S, const int& ldsy, const double* coeffs, double* k);
    void ApplyPeriodicBC(const int& Nx, const double* S, const int& ldx, double* dSdx, double* dSdy, const double* coeffs); 
    void EvaluateFuncBlas(const int& dimS, double* S, const double* dSdx, const double* dSdy, double* k);
    void WriteFile(const double* S);
public:
    // Constructors
    ShallowWater(); // Default Constructorn declaration
    ShallowWater(double dtt, double Tt, int Nxx, int Nyy, int icc, double dxx, double dyy, int typeAnalysis); // Constructor declaration
    
    // Methods
    void sayHello();
    void SetInitialCondition();
    void PrintMatrix(const int& N,  const double* A, const int& lda, const int& inc);
    void PrintVector(const int& N, const double* x);
    void TimeIntegrate();
    void TimeIntegrateForLoop();
    
    // 'Getter' functions
    double getTimeStep();
    double getIntegrationTime();
    int getNx();
    int getNy();
    int getIc();
    double getdx();
    double getdy();
    double* geth();
    double* getu();
    double* getv();
    
    ~ShallowWater(); // Destructor
    
};

struct Limits{
    double xl, xu, yl, yu;
public:
    Limits(double xa, double xb, double ya, double yb); // Constructor declaration
        
    ~Limits(); // Default destructor declaration
};


#endif
