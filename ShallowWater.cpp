#include <ShallowWater.h>
#include <iostream>

// constructor definition
ShallowWater::ShallowWater(){
}   // Default Constructor

ShallowWater::ShallowWater(double dtt, double Tt, int Nxx, int Nyy, int icc, double dxx, double dyy) : dt(dtt), T(Tt), Nx(Nxx), Ny(Nyy), ic(icc), dx(dxx), dy(dyy){
}   // Constructor using initialization list to avoid calling default class constructor and then over writting


ShallowWater::~ShallowWater(){
std::cout << "Class destroyed" << std::endl;
}   // Default destructor

// Method definition
void ShallowWater::sayHello(){
std::cout << "Hello from class 'ShallowWater'!" << std::endl;
}


// 'Getter' function definition
double ShallowWater::getTimeStep(){ return dt;}
double ShallowWater::getIntegrationTime(){ return T;}
int ShallowWater::getNx(){return Nx;}
int ShallowWater::getNy(){return Ny;}
int ShallowWater::getIc(){return ic;}
double ShallowWater::getdx(){ return dx;}
double ShallowWater::getdy(){return dy;}



// Structure 'limits'

 Limits::Limits(double xa = 0, double xb = 0, double ya = 0, double yb = 0): xl{xa}, xu{xb}, yl{ya}, yu{yb} {
    if(xa>xb) std::swap(xl, xu);
    if(ya>yb) std::swap(yl, yu);
}   // Default Constructor of structure 'Limit'
    
Limits::~Limits(){
    std::cout << "Structure destroyed" << std::endl;
}   // Default destructor