struct Limits{
    double xl, xu, yl, yu;
public:
    Limits(double xa, double xb, double ya, double yb); // Constructor declaration
        
    ~Limits(); // Default destructor declaration
};

class ShallowWater
{
    double dt, T;
    int Nx, Ny, ic;
    
//    Generate
public:
    ShallowWater(); // Default Constructorn declaration
    ShallowWater(double dtt, double Tt, int Nxx, int Nyy, int icc); // Constructor declaration
    
//    SetInitialCondition(
    
    
    double getTimeStep();
    double getIntegrationTime();
    int getNx();
    int getNy();
    int getIc();
    
    void sayHello();
    
    ~ShallowWater(); // Destructor
    
};