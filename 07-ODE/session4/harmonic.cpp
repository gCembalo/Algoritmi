// solve the harmonic oscillator problem:
// \dv[2]{x}{t} = -\omega^2 x
// using position-Verlet (or velocity-Verlet)
// algorithm, choosing a constant time step h = 0.02T and initial condition
// x0 = 1, v0 = 0 (T=2π/ω). You can use ω =1.
//
// Evolve the system for 10 periods and plot the mechanical energy of the
// system. Compare it with RK2 midpoint.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#define OMEGA 1.0 // pulsazione del problema
#define NMAX_EQ 64 // numero massimo di eq (sicurezza)

using namespace std;

void RHSFuncOde4(double, double *,double *);
void EulerStep(double, double *, void (*)(double, double *, double *),
               double, int);
void RK2StepMid(double, double *, void (*)(double, double *, double *),
                double, int);
void RK4Step(double, double *, void (*)(double, double *, double *), double, int);
void acceleration(double *, double *, int);
void PositionVerletStep(double *, double *, int, double, void (*)(double *,
                        double *, int));
void VelocityVerletStep(double *, double *, int, double, void (*)(double *,
                        double *, int));

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata;
    fdata.open("harmonic.dat"); // file per le soluzioni

    //set the significant digits and the scientific notation
    fdata << setprecision(7);
    fdata << setiosflags ( ios::scientific );

    int neq = 2; // ordine della ode
    int np = 1;
    int nstep; // numero di step
    double T = 2.0 * M_PI / OMEGA; // periodo
    double dt = 0.02*T; // step fisso
    double t; // variabile tempo
    double t0 = 0.0; // tempo iniziale
    double tf = 10.0*T; // tempo finale
    double E; // variabile energia

    nstep = fabs(tf - t0) / dt;

    // setto le condizioni iniziali
    double x0 = 1.0;
    double v0 = 0.0;

    // definisco l'array delle soluzioni e imposto le condizioni iniziali
    double Y[neq];
    double x[neq] , v[neq]; // array posizione e velocità per Verlet
    Y[0] = x0;
    Y[1] = v0;

    // stampo nel file (la condizione iniziale)
    fdata << t << " " << Y[0] << " " << Y[1] << endl;

    // richiamo i diversi metodi per risolvere il problema

    // ciclo per determinare la soluzione con Eulero
    for( int i = 0 ; i < nstep ; i++ ){

        // richiamo il metodo di Eulero per trovare le soluzioni
        EulerStep(t, Y, RHSFuncOde4, dt, neq);

        t += dt; // implemento il tempo

        E = 0.5*Y[1]*Y[1] + 0.5*Y[0]*Y[0]; // calcolo l'energia

        // stampo nel file
        fdata << t << "  " << Y[0] << "  " << Y[1] << "  " << E << endl;

    }

    fdata << endl << endl; // separo i diversi metodi

    // devo reinizializzare i valori iniziali della ODE e 
    // l'estremo di integrazione
    t = t0;
    Y[0] = x0;
    Y[1] = v0;

    // ciclo per determinare la soluzione con RK2Mid
    for( int i = 0 ; i < nstep ; i++ ){

        // richiamo il metodo di Eulero per trovare le soluzioni
        RK2StepMid(t, Y, RHSFuncOde4, dt, neq);

        t += dt; // implemento il tempo

        E = 0.5*Y[1]*Y[1] + 0.5*Y[0]*Y[0]; // calcolo l'energia

        // stampo nel file
        fdata << t << "  " << Y[0] << "  " << Y[1] << "  " << E << endl;

    }

    fdata << endl << endl; // separo i diversi metodi

    // devo reinizializzare i valori iniziali della ODE e l'estremo
    // di integrazione
    t = t0;
    Y[0] = x0;
    Y[1] = v0;

    // ciclo per determinare la soluzione con RK4
    for( int i = 0 ; i < nstep ; i++ ){

        // richiamo il metodo di Eulero per trovare le soluzioni
        RK4Step(t, Y, RHSFuncOde4, dt, neq);

        t += dt; // implemento il tempo

        E = 0.5*Y[1]*Y[1] + 0.5*Y[0]*Y[0]; // calcolo l'energia

        // stampo nel file
        fdata << t << "  " << Y[0] << "  " << Y[1] << "  " << E << endl;

    }

    fdata << endl << endl; // separo i diversi metodi

    // devo reinizializzare i valori iniziali della ODE e l'estremo
    // di integrazione
    t = t0;
    x[0] = x0;
    v[0] = v0;

    // ciclo per determinare la soluzione con VerletPosition
    for( int i = 0 ; i < nstep ; i++ ){

        // richiamo il metodo di Eulero per trovare le soluzioni
        PositionVerletStep(x, v, np, dt, acceleration);

        t += dt; // implemento il tempo

        E = 0.5*x[0]*x[0] + 0.5*v[0]*v[0]; // calcolo l'energia

        // stampo nel file
        fdata << t << "  " << x[0] << "  " << v[0] << "  " << E << endl;

    }

    fdata << endl << endl; // separo i diversi metodi

    // devo reinizializzare i valori iniziali della ODE e l'estremo
    // di integrazione
    t = t0;
    x[0] = x0;
    v[0] = v0;

    // ciclo per determinare la soluzione con VerletVelocity
    for( int i = 0 ; i < nstep ; i++ ){

        // richiamo il metodo di Eulero per trovare le soluzioni
        VelocityVerletStep(x, v, np, dt, acceleration);

        t += dt; // implemento il tempo

        E = 0.5*x[0]*x[0] + 0.5*v[0]*v[0]; // calcolo l'energia

        // stampo nel file
        fdata << t << "  " << x[0] << "  " << v[0] << "  " << E << endl;

    }

    fdata.close();

    return 0;

}


// definisco il Right-Hand-Side-Function (è problem dependent). 
// Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore)
// mi faccio dare R
void RHSFuncOde4(double t, double *Y, double *R){

    // Compute the right-hand side of the ODE (2 equation)
    double x = Y[0];
    double y = Y[1];

    R[0] = y;
    R[1] = - OMEGA * OMEGA * x;

}

// implemento il metodo di Eulero. Gli do in input la variabile
// di integrazione, un puntatore alle soluzioni, un puntatore al
//Right-Hand-Side-Function, lo step da utilizzare e la dimensionalità
// di Y, che non è altro che il numero di ODE di primo ordine che abbiamo.
void EulerStep(double t, double *Y, void (*RHSFunc)(double, double *, double *),
               double dt, int neq){
    
    // dt è lo step che utilizziamo per trovare la soluzione di dY/dt = rhs.
    // neq è il numero di ODE (dimensionalità di Y[])
    // *RHSFunc() punta al Right-Hand-Side-Function (in questo caso dYdt())
    
    int k; // variabile per visitare tutte le componenti di *Y
    double rhs[NMAX_EQ]; // per assicurarsi che rhs[] sia grande
                     // abbastanza (neq < 256)

    // calcolo il lato destro dell'equazione
    RHSFunc (t, Y, rhs);

    for( k = 0 ; k < neq ; k++ ){

        Y[k] += dt*rhs[k]; // Update solution array

    }

}

// implemento il metodo Runge-Kutta del secondo ordine (midpoint).
// gli do in input la variabile di integrazione, il puntatore alle soluzioni,
// il puntatore alla funzione del Right-Hand-Side-Function, l'incremento e
// l'ordine della ODE.
void RK2StepMid(double t, double *Y, void (*RHSFunc)(double t, double *Y, 
                double *R), double h, int neq){
    
    // definisco i vettori per gli step intermedi
    double Y1[NMAX_EQ], k1[NMAX_EQ], k2[NMAX_EQ];
    
    RHSFunc(t,Y,k1); // calcolo k1 con il RSH con t_n e Y_n

    // scrivo il ciclo per determinare Y_n + k1*h/2
    for( int i = 0 ; i < neq ; i++ ){
        
        Y1[i] = Y[i] + 0.5*h*k1[i];

    }

    RHSFunc(t+0.5*h,Y1,k2); // calcolo k2 con il RSH con t_n+h/2 e Y_n+k1*h/2
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + k2*h
    for (int i = 0 ; i < neq ; i++){
        
        Y[i] += h*k2[i];

    }

}

// implemento il metodo Runge-Kutta del quarto ordine.
// gli do in input la variabile di integrazione, il puntatore alle soluzioni,
// il puntatore alla funzione del Right-Hand-Side-Function, l'incremento e
// l'ordine della ODE.
void RK4Step(double t, double *Y, void (*RHSFunc)(double t, double *Y, double *R),
             double h, int neq){
    
    // definisco i vettori per gli step intermedi
    double Y1[NMAX_EQ], k1[NMAX_EQ], k2[NMAX_EQ], k3[NMAX_EQ], k4[NMAX_EQ];
    
    RHSFunc(t,Y,k1); // calcolo k1 con il RSH con t_n e Y_n

    // scrivo il ciclo per determinare Y_n + k1*h/2
    for( int i = 0 ; i < neq ; i++ ){
        
        Y1[i] = Y[i] + 0.5*h*k1[i];

    }

    RHSFunc(t+0.5*h,Y1,k2); // calcolo k2 con il RSH con t_n+h/2 e Y_n+k1*h/2
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + k2*h/2
    for (int i = 0 ; i < neq ; i++){
        
        Y1[i] = Y[i] + h*k2[i]*0.5;

    }

    RHSFunc(t+0.5*h,Y1,k3); // calcolo k3 con il RSH con t_n+h/2 e Y_n+k2*h/2
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + k3*h
    for (int i = 0 ; i < neq ; i++){
        
        Y1[i] = Y[i] + h*k3[i];

    }

    RHSFunc(t+h,Y1,k4); // calcolo k4 con il RSH con t_n+h e Y_n+k3*h
    
    // scrivo il ciclo per calcolare 
    // Y_{n+1} = Y_n + h/6 * ( k1 + 2*k2 + 2*k3 + k4 )
    for (int i = 0 ; i < neq ; i++){
        
        Y[i] += h * ( k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i] ) / 6.0;

    }

}

// implemento la funzione accelerazione che mi serve per il metodo di Verlet
// gli do in input gli array di posizione e accelerazione, oltre ad il numero di
// punti
void acceleration(double *x, double *a, int n){

    for( int i = 0 ; i < n ; i++ ){

        a[i] = - OMEGA * OMEGA * x[i];
    }

}

// implemento il metodo di Verlet per la posizione
// gli do in input i puntatori ai vettori posizione e velocità, la dimensione
// degli array (ordine della EDO), l'incremento degli intervalli e la funzione 
// accelerazione
void PositionVerletStep(double *x, double *v, int neq, double dt, void (*acceleration)(double *, double *, int)){

    double a[NMAX_EQ]; // creo il vettore per l'accelerazione

    // calcolo x per mezzi step
    for( int i = 0 ; i < neq ; i++ ){

        x[i] += 0.5 * dt * v[i];

    }

    // calcolo l'accelerazione a t = i*t + t/2
    acceleration( x, a, neq );

    // calcolo la velocità per step interi
    for( int i = 0 ; i < neq ; i++ ){

        v[i] += dt * a[i];

    }

    // calcolo la posizione per step intero
    for( int i = 0 ; i < neq ; i++ ){

        x[i] += 0.5 * dt * v[i];

    }

}

// implemento il metodo di Verlet per la posizione
// gli do in input i puntatori ai vettori posizione e velocità, la dimensione
// degli array (ordine della EDO), l'incremento degli intervalli e la funzione 
// accelerazione
void VelocityVerletStep(double *x, double *v, int neq, double dt, void (*acceleration)(double *, double *, int)){

    double a[NMAX_EQ]; // creo il vettore per l'accelerazione

    // calcolo l'accelerazione a t = i*t
    acceleration( x, a, neq );

    // calcolo v per mezzi step
    for( int i = 0 ; i < neq ; i++ ){

        v[i] += 0.5 * dt * a[i];

    }

    // calcolo la posizione  per step interi
    for( int i = 0 ; i < neq ; i++ ){

        x[i] += dt * v[i];

    }

    // calcolo l'accelerazione per x^{n+1}
    acceleration( x, a, neq );

    // calcolo v per step intero
    for( int i = 0 ; i < neq ; i++ ){

        v[i] += 0.5 * dt * a[i];

    }
    
}