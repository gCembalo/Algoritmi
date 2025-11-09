// using Euler, RK2 and RK4 methods, solve the system of ODEs
// 
// \dot{x} = \dv{x}{t} = y
// \dot{y} = \dv{y}{t} = -x
// 
// with initial condition x(0) = 1, y(0) = 0.
// 
// What is the analytical solution?
// 
// Try to integrate the system between 0 ≤ t ≤ 20π using 200 points and compare the solutions. Which method is the most accurate ?
//
// How can you choose the step size ?
//
// Is there any quantity that you expect to be conserved ? Are they ?
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double ode2Sol(double);
void EulerStep(double, double *, void (*)(double, double *, double *), double, int);
void RHSFuncOde2(double, double *,double *);
void RK2StepMid(double, double *, void (*)(double, double *, double *), double, int);
void RK2StepHeun(double, double *, void (*)(double, double *, double *), double, int);
void RK4Step(double, double *, void (*)(double, double *, double *), double, int);

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata1, fdata2;
    fdata1.open("ode2.dat"); // file per le soluzioni
    fdata2.open("ode2Conv.dat"); // file per la convergenza

    //set the significant digits and the scientific notation
    fdata1 << setprecision(7);
    fdata1 << setiosflags ( ios::scientific );
    fdata2 << setprecision(7);
    fdata2 << setiosflags ( ios::scientific );

    int neq = 2; // ordine della ode
    double Y[neq]; // definisco l'array delle soluzioni

    double t; // variabile tempo
    double tb = 0.0; // l'estremo di integrazione
    double te = 20.0*M_PI; // l'estremo di integrazione
    double dt; // incremento

    int npoint = 500; // numero di punti

    // definisco le condizioni iniziali del problema
    Y[0] = 1.0;
    Y[1] = 0.0;

    dt = fabs(tb - te) / (double)npoint; // calcolo gli incrementi
    t = tb; // inizializzo t al tempo iniziale

    // stampo nel file (la condizione iniziale)
    fdata1 << t << " " << Y[0] << " " << Y[1] << endl;

    // richiamo i diversi metodi per risolvere il problema

    // ciclo per determinare la soluzione con Eulero
    for( int i = 0 ; i < npoint ; i++ ){

        // richiamo il metodo di Eulero per trovare le soluzioni
        EulerStep(t, Y, RHSFuncOde2, dt, neq);

        t += dt; // implemento il tempo

        // stampo nel file
        fdata1 << t << "  " << Y[0] << "  " << Y[1] << endl;

    }

    fdata1 << endl << endl; // separo i diversi metodi

    // devo reinizializzare i valori iniziali della ODE e l'estremo di integrazione
    t = tb;
    Y[0] = 1.0;
    Y[1] = 0.0;

    // ciclo per determinare la soluzione con RK2Mid
    for( int i = 0 ; i < npoint ; i++ ){

        // richiamo il metodo di Eulero per trovare le soluzioni
        RK2StepMid(t, Y, RHSFuncOde2, dt, neq);

        t += dt; // implemento il tempo

        // stampo nel file
        fdata1 << t << "  " << Y[0] << "  " << Y[1] << endl;

    }

    fdata1 << endl << endl; // separo i diversi metodi

    // devo reinizializzare i valori iniziali della ODE e l'estremo di integrazione
    t = tb;
    Y[0] = 1.0;
    Y[1] = 0.0;

    // ciclo per determinare la soluzione con Heun
    for( int i = 0 ; i < npoint ; i++ ){

        // richiamo il metodo di Eulero per trovare le soluzioni
        RK2StepHeun(t, Y, RHSFuncOde2, dt, neq);

        t += dt; // implemento il tempo

        // stampo nel file
        fdata1 << t << "  " << Y[0] << "  " << Y[1] << endl;

    }

    fdata1 << endl << endl; // separo i diversi metodi

    // devo reinizializzare i valori iniziali della ODE e l'estremo di integrazione
    t = tb;
    Y[0] = 1.0;
    Y[1] = 0.0;

    // ciclo per determinare la soluzione con RK4
    for( int i = 0 ; i < npoint ; i++ ){

        // richiamo il metodo di Eulero per trovare le soluzioni
        RK4Step(t, Y, RHSFuncOde2, dt, neq);

        t += dt; // implemento il tempo

        // stampo nel file
        fdata1 << t << "  " << Y[0] << "  " << Y[1] << endl;

    }

    fdata1.close();


    // faccio lo studio della convergenza (rifaccio tutti i cicli visto che cambio dominio di integrazione e numero di punti)
    // uso la variabile tempo di prima
    // l'estremo di integrazione t_begin è lo stesso di sopra
    double te1 = 3; // l'estremo di integrazione
    double te2 = 2*M_PI; // l'estremo di integrazione
    double dtc; // incremento

    // definisco gli array delle soluzioni
    double Y1[neq], Y2[neq], Y3[neq], Y4[neq];

    double SolEx; // soluzione esatta
    double errEul, errRK2M, errRK2H, errRK4; // errore

    int npointC[10] = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048}; // numero di punti

    for( int k = 0 ; k < 10 ; k ++ ){

        // devo reinizializzare i valori iniziali della ODE e l'estremo di integrazione
        t = tb;
        Y1[0] = 1.0 , Y2[0] = 1.0 , Y3[0] = 1.0 , Y4[0] = 1.0;
        Y1[1] = 0.0 , Y2[1] = 0.0 , Y3[1] = 0.0 , Y4[1] = 0.0;

        // definisco gli step degli intervalli
        dtc = fabs(tb - te1) / ( (double)npointC[k] );

        //SolEx = ode2Sol(t); // soluzione esatta

        //err = fabs( Y[0] - SolEx ); // calcolo gli errori

        // stampo nel file (la condizione iniziale)
        //fdata2 << t << "  " << Y[0] << "  " << Y[1] << err << "   " << endl;

        for( int i = 0 ; i < npointC[k] ; i++ ){

            // richiamo i metodo per trovare le soluzioni
            EulerStep(t, Y1, RHSFuncOde2, dtc, neq);
            RK2StepMid(t, Y2, RHSFuncOde2, dtc, neq);
            RK2StepHeun(t, Y3, RHSFuncOde2, dtc, neq);
            RK4Step(t, Y4, RHSFuncOde2, dtc, neq);

            t += dtc; // implemento il tempo

        }

        // calcolo gli errori
        SolEx = ode2Sol(te1); // soluzione esatta

        errEul = fabs( Y1[0] - SolEx );
        errRK2M = fabs( Y2[0] - SolEx );
        errRK2H = fabs( Y3[0] - SolEx );
        errRK4 = fabs( Y4[0] - SolEx );

            // stampo nel file
            fdata2 << dtc << "  " << errEul << "  " << errRK2M << "  " << errRK2H << "   " << errRK4 << endl;

    }

    fdata2.close();


    return 0;

}

double ode2Sol(double t){

    return cos(t);

}

// implemento il metodo di Eulero. Gli do in input la variabile di integrazione, un puntatore alle soluzioni, un puntatore al Right-Hand-Side-Function, lo step da utilizzare e la dimensionalità di Y, che non è altro che il numero di ODE di primo ordine che abbiamo.
void EulerStep(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq){
    
    // dt è lo step che utilizziamo per trovare la soluzione di dY/dt = rhs.
    // neq è il numero di ODE (dimensionalità di Y[])
    // *RHSFunc() punta al Right-Hand-Side-Function (in questo caso dYdt())
    
    int k; // variabile per visitare tutte le componenti di *Y
    double rhs[256]; // per assicurarsi che rhs[] sia grande abbastanza (neq < 256)

    // calcolo il lato destro dell'equazione
    RHSFunc (t, Y, rhs);

    for( k = 0 ; k < neq ; k++ ){

        Y[k] += dt*rhs[k]; // Update solution array

    }

}

// definisco il Right-Hand-Side-Function (è problem dependent). Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore) mi faccio dare R
void RHSFuncOde2(double t, double *Y, double *R){

    // Compute the right-hand side of the ODE (2 equation)
    double x = Y[0];
    double y = Y[1];
    R[0] = y;
    R[1] = -x;

}

// implemento il metodo Runge-Kutta del secondo ordine (midpoint).
// gli do in input la variabile di integrazione, il puntatore alle soluzioni, il puntatore alla funzione del Right-Hand-Side-Function, l'incremento e l'ordine della ODE.
void RK2StepMid(double t, double *Y, void (*RHSFunc)(double t, double *Y, double *R), double h, int neq){
    
    // definisco i vettori per gli step intermedi
    double Y1[neq], k1[neq], k2[neq];
    
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

// implemento il metodo Runge-Kutta del secondo ordine (modified Eulero).
// gli do in input la variabile di integrazione, il puntatore alle soluzioni, il puntatore alla funzione del Right-Hand-Side-Function, l'incremento e l'ordine della ODE.
void RK2StepHeun(double t, double *Y, void (*RHSFunc)(double t, double *Y, double *R), double h, int neq){
    
    // definisco i vettori per gli step intermedi
    double Y1[neq], k1[neq], k2[neq];
    
    RHSFunc(t,Y,k1); // calcolo k1 con il RSH con t_n e Y_n

    // scrivo il ciclo per determinare Y_n + k1*h
    for( int i = 0 ; i < neq ; i++ ){
        
        Y1[i] = Y[i] + h*k1[i];

    }

    RHSFunc(t+h,Y1,k2); // calcolo k2 con il RSH con t_n+h e Y_n+k1*h
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + (k1+k2)*h/2
    for (int i = 0 ; i < neq ; i++){
        
        Y[i] += ( k1[i] + k2[i] )*h*0.5;

    }

}

// implemento il metodo Runge-Kutta del quarto ordine.
// gli do in input la variabile di integrazione, il puntatore alle soluzioni, il puntatore alla funzione del Right-Hand-Side-Function, l'incremento e l'ordine della ODE.
void RK4Step(double t, double *Y, void (*RHSFunc)(double t, double *Y, double *R), double h, int neq){
    
    // definisco i vettori per gli step intermedi
    double Y1[neq], k1[neq], k2[neq], k3[neq], k4[neq];
    
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
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + h/6 * ( k1 + 2*k2 + 2*k3 + k4 )
    for (int i = 0 ; i < neq ; i++){
        
        Y[i] += h * ( k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i] ) / 6.0;

    }

}