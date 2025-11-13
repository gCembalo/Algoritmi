// Consider the ODE
// dy/dt = -ty ; y(0) = 1
// This has analytical solution y = exp( -t*t/2 )
//
// Implement Euler’s method and try integrating from t = 0 up to t = 3 using
// different step sizes, e.g., h=0.5, 0.2, 0.1, 0.05, ..., 0.001.
//
// The following table shows the numerical solution and errors for h = 0.5.
//
// Write ASCII data file showing the absolute and relative errors for each step
// size as a function of position. Comment on this.
//
// Choose h = 1. What happens ?
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

void EulerStep(double, double *, void (*)(double, double *, double *), 
               double, int);
void RHSFuncOde1(double, double *,double *);
double ode1Sol(double);

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata1, fdata2;
    fdata1.open("ode1.dat"); // file per soluzione metodo di Eulero
    fdata2.open("ode1Ex.dat"); // file per soluzione esatta

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
    // poi metterò le componenti in un ciclo così da variarlo
    double step[] = {0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001}; 

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

    // ciclo che ogni volta varia gli step (vedi array step[] )
    for( int k = 0 ; k < 9 ; k ++ ){

        // tempo di integrazione
        tb = 0.0;
        te = 3.0;

        // definisco gli step degli intervalli
        dt = step[k];
        nstep = ceil( fabs(tb - te) / dt );

        // definisco la condizione iniziale
        t = tb;
        Y[0] = 1.0; // condizione iniziale del problema

        SolEx = ode1Sol(t); // soluzione esatta

        // calcolo gli errori
        abs_err = fabs( Y[0] - SolEx );
        rel_err = fabs( abs_err / SolEx );

        // stampo nel file (la condizione iniziale)
        fdata1 << t << "  " << Y[0] << "  " << abs_err << "   " 
               << rel_err << endl;

        // ciclo per determinare la soluzione con Eulero
        for( int i = 0 ; i < nstep ; i++ ){

            // richiamo il metodo di Eulero per trovare le soluzioni
            EulerStep(t, Y, RHSFuncOde1, dt, neq);

            t += dt; // implemento il tempo

            SolEx = ode1Sol(t); // soluzione esatta

            // calcolo gli errori
            abs_err = fabs( Y[0] - SolEx );
            rel_err = fabs( abs_err/SolEx );

            // stampo nel file
            fdata1 << t << "  " << Y[0] << "  " << abs_err << "   " 
                   << rel_err << endl;
        }

        // separo i diversi incrementi
        fdata1 << endl << endl;

    }

    fdata1.close();

    return 0;

}

// implemento il metodo di Eulero. Gli do in input la variabile di integrazione,
// un puntatore alle soluzioni, un puntatore al Right-Hand-Side-Function, lo
// step da utilizzare e la dimensionalità di Y, che non è altro che il numero
// di ODE di primo ordine che abbiamo.
void EulerStep(double t, double *Y, void (*RHSFunc)(double, double *, double *),
               double dt, int neq){
    
    // dt è lo step che utilizziamo per trovare la soluzione di dY/dt = rhs.
    // neq è il numero di ODE (dimensionalità di Y[])
    // *RHSFunc() punta al Right-Hand-Side-Function (in questo caso dYdt())
    
    int k; // variabile per visitare tutte le componenti di *Y
    double rhs[256]; // per assicurarsi che rhs[] sia grande
                     // abbastanza (neq < 256)

    // calcolo il lato destro dell'equazione
    RHSFunc (t, Y, rhs);

    for( k = 0 ; k < neq ; k++ ){

        Y[k] += dt*rhs[k]; // Update solution array

    }

}

// funzione problem dependent
void RHSFuncOde1(double t, double *Y, double *R){

    // Compute the right-hand side of the ODE dy/dt = -t*y
    double y = Y[0];
    R[0] = -t*y;

}

// soluzione esatta
double ode1Sol(double t){

    return exp( -t*t/2 );

}