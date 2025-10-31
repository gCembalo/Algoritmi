//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

void EulerStep (double, double *, void (*)(double, double *, double *), double, int);
void dYdt (double, double *,double *);

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata1, fdata2;
    fdata1.open("ode2Eu.dat");
    // 
    fdata2.open("ode2RK.dat");

    //set the significant digits and the scientific notation
    fdata1 << setprecision(7);
    fdata1 << setiosflags ( ios::scientific );
    fdata2 << setprecision(7);
    fdata2 << setiosflags ( ios::scientific );

    int neq = 2; // ordine della ode
    double Y[neq]; // definisco l'array delle soluzioni

    return 0;

}

void EulerStep (double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq){
    
    // Take one step dt using Euler method for the solution of dY/dt = rhs.
    // Here neq is the number of ODE (the dimensionality of Y[]) and *RHSFunc() is a
    // pointer to the actual function (in this case it points to dYdt()) that
    // calculates the right-hand side of the system of equations.
    
    int k; // variabile per visitare tutte le componenti di *Y
    double rhs[256]; // Make sure rhs[] is large enough (this implies neq < 256)
    // rhs[neq] is *NOT* standard C++ (forbidden as static declaration)

    RHSFunc (t, Y, rhs);

    for( k = 0 ; k < neq ; k++ ){

        Y[k] += dt*rhs[k]; // Update solution array

    }

}

void RK2Step(double t, double *Y, double (*RHS_Func)(double t, double *Y, double *R), double h, int neq){
    
    double Y1[neq], k1[neq], k2[neq];
    
    RHS_Func(t,Y,k1);

    for( int i = 0 ; i < neq-1 ; i++ ){
        
        Y1[i] = Y[i]+0.5*h*k1[i]

    }

    RHS_Func(t+0.5*h,Y1,k2);
    
    for (int i = 0 ; i < neq-1 ; i++){
        
        Y[i] += h*k2[i]

    }

}

void dYdt (double t, double *Y, double *R){

    // Compute the right-hand side of the ODE (2 equation)
    double y = Y[0];
    double x = 
    R[0] = y;
    R[1] = -x;

}