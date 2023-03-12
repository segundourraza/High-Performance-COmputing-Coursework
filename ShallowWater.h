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
    \
    double* h = nullptr;
    double* u = nullptr;
    double* v = nullptr;
    double* f = nullptr;
    
public:
    // Constructors
    ShallowWater(); // Default Constructorn declaration
    ShallowWater(double dtt, double Tt, int Nxx, int Nyy, int icc, double dxx, double dyy); // Constructor declaration
    
    // Methods
    void sayHello();
    void SetInitialCondition();
    double GetDerivatives(double* vect, const double& step);
    void PrintMatrix(const int& N,  double* A, const int& lda);
    void EvaluateFuncBLAS(double* uu, double* vv, double* hh, double* f);
    void TimeIntegrate();
    
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

#endif