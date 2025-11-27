// Name: Gabriele, Cembalo
// Date: 26 Nov 2025
//
// Code output:
// *****************************************************
// Loop break at nstep = ??; t = ??     # when you exit from the loop
//               eps1 = ??; eps2 = ??   # differences between the two solutions
//               ip1 = ??; ip2 = ??     # number of inversion points
// *****************************************************

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

#define NMAX_EQ 64 // numero massimo di eq (sicurezza)

// Definizione funzioni
void RHSFuncOde(double, double *,double *);
void acceleration(double *, double *, int);
void RK4Step(double, double *, void (*)(double, double *, double *), double, int);
void PositionVerletStep(double *, double *, int, double, void (*)(double *,
                        double *, int));
void VelocityVerletStep(double *, double *, int, double, void (*)(double,
                        double *, int));

int main(){

    //cout << setprecision(7);
    //cout << setiosflags ( ios::scientific );

    ofstream fdata;
    fdata.open("coupled_rotation.dat"); // file per le soluzioni

    // set the significant digits and the scientific notation
    fdata << setprecision(7);
    fdata << setiosflags ( ios::scientific );

    double m = 1.0; // massa
    double R1 = 1.0 , R2 = 2.0; // raggi circonferenze
    double k = 1.2; // costante elastica della molla
    double L0 = R2 - R1; // lunghezza a riposo

    int nloop = 0; // conto il numero di loop

    // Definisco i parametri della ODE
    int neq = 4; // equazioni
    double tol = 0.1; // tolleranza
    // Condizioni iniziali
    double theta10 = 0.0 , theta20 = M_PI/3.0;
    double omega10 = 0.0 , omega20 = 0.05;
    // Definisco l'array delle soluzioni e imposto le condizioni iniziali.
    double Y[neq];
    // array posizione e velocità per Verlet
    // metto entrambe le soluzioni nello stesso vettore
    double x[2] , v[2];
    Y[0] = theta10 , Y[1] = omega10 , Y[2] = theta20 , Y[3] = omega20;
    x[0] = theta10 , x[1] = theta20;
    v[0] = omega10 , v[1] = omega20;

    double Deltat = 0.2; // intervallo
    double t0 = 0.0;
    double tf = 50.0;
    double t = t0; // inizializzo il tempo
    int nstep = fabs( (tf - t0) )/Deltat;

    double vold1 , vold2; // variabile per salvare la vecchia velocità
    int count1 = 0 , count2 = 0; // conteggi per i punti di inversione

    // differenze tra le soluzioni (li inizializzo così da iniziare il ciclo)
    double epsilon1 = 1.0e-3 , epsilon2 = 1.0e-3;

    // stampo nel file (la condizione iniziale)
    fdata << t << " " << Y[0] << " " << Y[1] << " " 
          << Y[2] << " " << Y[3] << endl;

    // richiamo i diversi metodi per risolvere il problema e blocco il loop
    // quando supero la tolleranza
    for( int i = 0 ; i < nstep ; i++ ){

        while( max(epsilon1, epsilon2) < tol ){

            nloop ++; // aumento il conteggio

            // salvo le vecchie velocità
            vold1 = Y[1];
            vold2 = Y[3];

            // richiamo il metodo RK4 per trovare le soluzioni
            RK4Step(t, Y, RHSFuncOde, Deltat, neq);

            // richiamo il metodo di PositionVerlet per trovare le soluzioni
            PositionVerletStep(x, v, 2, Deltat, acceleration);

            // calcolo le differenze
            epsilon1 = abs( Y[0] - x[0] )/( M_PI );
            epsilon2 = abs( Y[2] - x[1] )/( M_PI );

            t += Deltat; // implemento il tempo

            // controllo e conteggio il turning point (se presente)
            if( Y[1]*vold1 < 0 ){

                count1 ++;

            }

            if( Y[3]*vold2 < 0 ){

                count2 ++;

            }

            // stampo nel file
            fdata << t << "  " << Y[0] << " " << Y[2] << " " 
                           << x[0] << " " << x[1] << endl;

        }

    }

    cout << " *****************************************************" << endl;
    cout << "Loop break at nstep = " << nloop << "; t = " << t << "\n"
         << "               eps1 = " << epsilon1 << "; eps2 = " << epsilon2
         << "\n" << "                ip1 = " << count1 << "; ip2 = " << count2 
         << endl;
    cout << " *****************************************************" << endl;


    // separo l'output così posso rifare andare il ciclo e graficare tra 0 e 50
    fdata << endl << endl;

    // risetto le condizioni iniziali
    Y[0] = theta10 , Y[1] = omega10 , Y[2] = theta20 , Y[3] = omega20;
    x[0] = theta10 , x[1] = theta20;
    v[0] = omega10 , v[1] = omega20;
    t = t0; // inizializzo il tempo

    // stampo nel file (la condizione iniziale)
    fdata << t << " " << Y[0] << " " << Y[1] << " " 
          << Y[2] << " " << Y[3] << endl;

    // richiamo i diversi metodi per risolvere il problema e blocco il loop
    // quando supero la tolleranza
    for( int i = 0 ; i < nstep ; i++ ){

        // richiamo il metodo RK4 per trovare le soluzioni
        RK4Step(t, Y, RHSFuncOde, Deltat, neq);

        // richiamo il metodo di PositionVerlet per trovare le soluzioni
        PositionVerletStep(x, v, 2, Deltat, acceleration);

        t += Deltat; // implemento il tempo

        // stampo nel file
        fdata << t << "  " << Y[0] << " " << Y[2] << " " 
                           << x[0] << " " << x[1] << endl;

    }

    return 0;

}


// ------------------------ Funzioni implementate ------------------------ //

// Definisco il Right-Hand-Side-Function (è problem dependent). 
// Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore)
// mi faccio dare R
void RHSFuncOde(double t, double *Y, double *R){

    // definisco i parametri delle circonferenze
    double R1 = 1.0 , R2 = 2.0;
    double k = 1.2; // costante elastica della molla
    double L0 = R2 - R1; // lunghezza a riposo

    // calcolo la lunghezza L
    double L = sqrt( R1*R1 + R2*R2 - 2.0*R1*R2*cos( Y[0] - Y[2] ) );

    R[0] = Y[1];
    R[1] = - k * ( 1 - L0/L ) * ( R2/R1 ) * ( sin(Y[0]) - sin(Y[2]) );
    R[2] = Y[3];
    R[3] = + k * ( 1 - L0/L ) * ( R1/R2 ) * ( sin(Y[0]) - sin(Y[2]) );

}

// Implemento la funzione accelerazione che mi serve per il metodo di Verlet
// gli do in input gli array di posizione e accelerazione, oltre ad il numero di
// punti
void acceleration(double *theta, double *a, int n){

    // il vettore con gli angoli
    double theta1 = theta[0];
    double theta2 = theta[1];

    // definisco i parametri delle circonferenze
    double R1 = 1.0 , R2 = 2.0;
    double k = 1.2; // costante elastica della molla
    double L0 = R2 - R1; // lunghezza a riposo

    // calcolo la lunghezza L
    double L = sqrt( R1*R1 + R2*R2 - 2.0*R1*R2*cos( theta1 - theta2 ) );

    // accelerazione per theta1
    a[0] = - k * ( 1 - L0/L ) * ( R2/R1 ) * ( sin(theta1) - sin(theta2) );
    // accelerazione per theta2
    a[1] = + k * ( 1 - L0/L ) * ( R1/R2 ) * ( sin(theta1) - sin(theta2) );

}

// Implemento il metodo Runge-Kutta del quarto ordine.
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

// Implemento il metodo di Verlet per la posizione
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

// Implemento il metodo di Verlet per la posizione
// gli do in input i puntatori ai vettori posizione e velocità, la dimensione
// degli array (ordine della EDO), l'incremento degli intervalli e la funzione 
// accelerazione
void VelocityVerletStep(double *x, double *v, int neq, double dt,
                        void (*acceleration)(double *, double *, int)){

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