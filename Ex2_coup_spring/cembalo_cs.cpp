// Name: Gabriele, Cembalo
// Date: 27 Nov 2025
//
// Code output:
// *****************************************************
// Number of crossing[RK4] = 11
// Rmin[RK4]:                5.0210829e-01
// Number of crossing[Verlet] = 11
// Rmin[Verlet]:             5.0210977e-01
// *****************************************************

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

#define NMAX_EQ 64 // numero massimo di eq (sicurezza)

void RHSFuncOde(double, double *,double *);
void acceleration(double *, double *, int);
void RK4Step(double, double *, void (*)(double, double *, double *), double, int);
void PositionVerletStep(double *, double *, int, double, void (*)(double *,
                        double *, int));

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata;
    fdata.open("coupled_springs.dat"); // file per le soluzioni

    //set the significant digits and the scientific notation
    fdata << setprecision(7);
    fdata << setiosflags ( ios::scientific );

    double L0 = 1.0 , k = 6.0; // lunghezza riposo e costante elastica
    double kb = 1.0 , alphab = 0.1;
    int npoint = 800; // numero di punti
    double t , t0 = 0.0 , tf = 30.0;
    double dt = fabs( tf - t0 )/(double)npoint;
    double x0 = 1.0 , y0 = 0.5 , xb0 = 1.0; // condizioni iniziali
    double vx0 = 0.0 , vy0 = 0.0 , vxb0 = 0.0; // condizioni iniziali

    int neq = 6; // numero ODE
    double Y[neq];
    double x[3] , v[3];
    Y[0] = x0 , Y[1] = vx0;
    Y[2] = y0 , Y[3] = vy0;
    Y[4] = xb0 , Y[5] = vxb0;
    x[0] = x0 , x[1] = y0 , x[2] = xb0;
    v[0] = vx0 , v[1] = vy0 , v[2] = vxb0;

    double yold1 , yold2; // variabile per salvare la posizione vecchia
    int ycountRK = 0; // variabile per contare quante volte mA attraversa x
    int ycountV = 0;

    double Rrk , Rv;
    double RrkMin = 1000.0 , RvMin = 1000.0; // varibili che uso per salvare il minimo

    // stampo nel file (la condizione iniziale)
    fdata << t << " " << Y[0] << " " << Y[1] << " " 
          << Y[2] << " " << Y[3] 
          << Y[4] << " " << Y[5] << endl;

    // richiamo i diversi metodi per risolvere il problema
    for( int i = 0 ; i < npoint ; i++ ){

        // salvo le vecchie posizioni
        yold1 = Y[2];
        yold2 = x[1];

        RK4Step(t, Y, RHSFuncOde, dt, neq);
        PositionVerletStep(x, v, 3, dt, acceleration);

        t += dt; // implemento il tempo

        if( Y[2]*yold1 < 0 ) ycountRK ++;
        if( x[1]*yold2 < 0 ) ycountV ++;

        Rrk = sqrt( ( Y[0] - Y[4] )*( Y[0] - Y[4] ) + Y[2]*Y[2] );
        Rv = sqrt( ( x[0] - x[2] )*( x[0] - x[2] ) + x[1]*x[1] );

        // stampo nel file
        fdata << t << "  " << Y[0] << " " << Y[1] << " " << Y[2] << " "
                           << Y[3] << " " << Y[4] << " " << Y[5] << " "
                           //<< x[0] << " " << x[1] << " " << x[2] << " "
                           << Rrk << " " << Rv << endl;

        if( Rrk < RrkMin ) RrkMin = Rrk;
        if( Rv < RvMin ) RvMin = Rv;

    }

    cout << " *****************************************************" << endl;
    cout << "Number of crossing[RK4] = " << ycountRK << "\n"
         << "Rmin[RK4]:                " << RrkMin << "\n" 
         << "Number of crossing[Verlet] = " << ycountV << "\n"
         << "Rmin[Verlet]:             " << RvMin
         << endl;
    cout << " *****************************************************" << endl;

    return 0;

}


// ------------------------ Funzioni implementate ------------------------ //

// definisco il Right-Hand-Side-Function (è problem dependent). 
// Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore)
// mi faccio dare R
void RHSFuncOde(double t, double *Y, double *R){

    double L0 = 1.0 , k = 6.0; // lunghezza riposo e costante elastica
    double kb = 1.0 , alphab = 0.1;

    // x = Y[0]
    // vx = Y[1]
    // y = Y[2]
    // vy = Y[3]
    // xb = Y[4]
    // vxb = Y[5]

    // la definisco minuscola per non confonderla con *R
    double r = sqrt( ( Y[0] - Y[4] )*( Y[0] - Y[4] ) + Y[2]*Y[2] );

    R[0] = Y[1];
    R[1] = - k * ( r - L0 ) * ( Y[0] - Y[4] ) / r;
    R[2] = Y[3];
    R[3] = - k * (r - L0 ) * Y[2] / r;
    R[4] = Y[5];
    R[5] = - kb * Y[4] + alphab * Y[4] * Y[4];

}

// implemento la funzione accelerazione che mi serve per il metodo di Verlet
// gli do in input gli array di posizione e accelerazione, oltre ad il numero di
// punti
void acceleration(double *x, double *a, int n){

    double L0 = 1.0 , k = 6.0; // lunghezza riposo e costante elastica
    double kb = 1.0 , alphab = 0.1;

    // x = x[0]
    // y = x[1]
    // xb = x[3]

    // la definisco minuscola per non confonderla con *R
    double r = sqrt( ( x[0] - x[2] )*( x[0] - x[2] ) + x[1]*x[1] );

    a[0] = - k * ( r - L0 )*( x[0] - x[2] )/r;
    a[1] = - k * ( r - L0 )*( x[1] )/r;
    a[2] = - kb*x[2] + alphab*x[2]*x[2];

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

// implemento il metodo di Verlet per la posizione
// gli do in input i puntatori ai vettori posizione e velocità, la dimensione
// degli array (ordine della EDO), l'incremento degli intervalli e la funzione 
// accelerazione
void PositionVerletStep(double *x, double *v, int neq, double dt,
                        void (*acceleration)(double *, double *, int)){

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