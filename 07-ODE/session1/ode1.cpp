// Consider the ODE
// dy/dt = -ty ; y(0) = 1
// This has analytical solution y = exp( -t*t/2 )
//
// Implement Euler’s method and try integrating from t = 0 up to t = 3 using different step sizes, e.g., h=0.5, 0.2, 0.1, 0.05, ..., 0.001.
//
// The following table shows the numerical solution and errors for h = 0.5.
//
// Write ASCII data file showing the absolute and relative errors for each step size as a function of position. Comment on this.
//
// Choose h = 1. What happens ?
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

void EulerStep (double, double *, void (*)(double, double *, double *), double, int);
void dYdt (double, double *,double *);
double ode1Sol(double);

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata1, fdata2;
    fdata1.open("ode1.dat");
    // 
    fdata2.open("ode1Ex.dat");

    //set the significant digits and the scientific notation
    fdata1 << setprecision(7);
    fdata1 << setiosflags ( ios::scientific );
    fdata2 << setprecision(7);
    fdata2 << setiosflags ( ios::scientific );

    // stampo i valori della soluzione analitica
    for( double ti = 0.0 ; ti <= 3 ; ti += 1.e-3 ){

        fdata2 << ti << "     " << ode1Sol(ti) << endl;

    }
    fdata2.close();

    // definisco la spaziatura
    double step[] = {0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001}; // poi metterò le componenti in un ciclo così da variarlo

    int neq = 1; // ordine della ode

    double t; // variabile tempo
    double tb; // l'estremo di integrazione
    double te; // l'estremo di integrazione
    double dt; // incremento
    double nstep; // numero di step

    // definisco l'array delle soluzioni
    double Y[neq];

    // definisco le variabili degli errori
    double abs_err;
    double rel_err;

    // definisco la variabile per la soluzione esatta
    double SolEx;

    for( int k = 0 ; k < 9 ; k ++ ){

        tb = 0.0;
        te = 3.0;

        dt = step[k];
        nstep = ceil( fabs(tb - te) / dt );

        // definisco la condizione iniziale
        t = tb;
        Y[0] = 1.0; // condizione iniziale del problema

        SolEx = ode1Sol(t); // soluzione esatta

        // calcolo gli errori
        abs_err = fabs( Y[0] - SolEx );
        rel_err = fabs( abs_err / SolEx );

        // fdata1 << t << "  " << Y[0] << "  " << fabs( Y[0] - exp( -t*t/2.0 ) ) << "   " << fabs( Y[0]/exp( -t*t/2.0 ) - 1.0) << endl;
        fdata1 << t << "  " << Y[0] << "  " << abs_err << "   " << rel_err << endl;

        for( int i = 0 ; i < nstep ; i++ ){

            EulerStep(t, Y, dYdt, dt, neq);
            //SolEx = ode1Sol(t);

            t += dt;

            SolEx = ode1Sol(t); // soluzione esatta

            // calcolo gli errori
            abs_err = fabs( Y[0] - SolEx );
            rel_err = fabs( abs_err/SolEx );

            // fdata1 << t << "  " << Y[0] << "  " << fabs( Y[0] - exp( -t*t/2.0 ) ) << "   " << fabs( Y[0]/exp( -t*t/2.0 ) - 1.0) << endl;
            fdata1 << t << "  " << Y[0] << "  " << abs_err << "   " << rel_err << endl;
        }

        /*
        // ciclo per h = 0.5
        for( t = 0.0 ; t <= te ; t += h[k] ){

            EulerStep(t, Y, dYdt, h[k], neq);
            SolEx = ode1Sol(t);

            // calcolo gli errori
            abs_err = fabs( Y[0] - SolEx );
            rel_err = abs_err/SolEx;

            // stampo a terminale
            cout << t << "  " << Y[0] << "  " << abs_err << "   " << rel_err << endl;

            // stampo nel file
            fdata1 << t << "     " << Y[0] << endl;

        }*/

        // separo i diversi incrementi
        fdata1 << endl << endl;

    }

    fdata1.close();

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

void dYdt (double t, double *Y, double *R){

    // Compute the right-hand side of the ODE dy/dt = -t*y
    double y = Y[0];
    R[0] = -t*y;

}

double ode1Sol(double t){

    return exp( -t*t/2 );

}