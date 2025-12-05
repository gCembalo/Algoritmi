// Compute the potential of an infinitely long charged cylinder by solving the Poisson equation
//
//          \nabla^2 \phi = -\rho
// (see slide)
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double SolCylinder(double, double, double);
void Boundary(double **, double *, double *, int, double);
int JacobiMethod(double **, double **, double *, double *, double, double, 
                    double, int, int);
int GaussSeidelMethod(double **, double **, double *, double *, double, double, 
                    double, int, int);
int SORMethod(double **, double **, double *, double *, double, double, 
                    double, int, int);

int main(){

    // Implementazione main()

    return 0;

}


// ------------------------ Funzioni implementate ------------------------ //

// Implemento la funzione soluzione
double SolCylinder(double x, double y, double a){

    double rho0 = 1.0;

    double r = sqrt( x*x + y*y ); // coordinata radiale

    if( r <= a ){
        return -rho0*r*r*0.25;
    }
    if( r > a){
        return -rho0*a*a*0.5*( log(r/a) + 0.5 );
    }

}

// Implemento una funzione che mi fissa le condizioni al contorno (alla griglia)
// Gli do in input il vettore soluzione, gli array x e y della griglia, il
// numero di punti ed il valore di S
void Boundary(double **phi, double *x, double *y, int nx, int ny, double a){

    // Costruisco il ciclo per fissare le condizioni (usando la soluzione)

    // bordo sinistro e destro
    for (int i = 0; i < nx; i++) {
        phi[i][0]     = SolCylinder(x[i], y[0], a);
        phi[i][ny-1]  = SolCylinder(x[i], y[ny-1], a);
    }

    // bordo inferiore e superiore
    for (int j = 0; j < ny; j++) {
        phi[0][j]     = SolCylinder(x[0], y[j], a);
        phi[nx-1][j]  = SolCylinder(x[nx-1], y[j], a);
    }

}

// Implemento il metodo di Jacobi
// gli do in input gli array bidimensionale delle soluzioni, la griglia,
// la soluzione, la spaziatura della griglia, la tolleranza e
// i punti della griglia.
// Restituisce il numero di iterazioni.
int JacobiMethod(double **phi0, double **phi1, double *x, double *y, 
                  double S, double h, double tol, int nx, int ny){

    double sum; // variabile somma in cui metto la soluzione ad ogni step

    double err = 1.e3; // variabile errore che utilizzo per bloccare il ciclo
    double deltax = 0. , deltay = 0. ; // variabili che mi servono per l'errore

    int iter = 0; // variabile per contare le iterazioni

    while( err > tol ){

        sum = 0.0; // azzero la variabile

        // setto le condizioni al bordo
        Boundary(phi0, x, y, nx, ny, S);

        // copio la vecchia soluzione phi0 nella nuova phi1
        for(int i = 0 ; i < nx ; i++ ){

            for( int j = 0 ; j < ny ; j++ ){

                phi0[i][j] = phi1[i][j];

            }

        }

        // implemento il metodo iterativo (ricordando che il primo e ultimo
        // elemento sono fissati dal bordo)
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                phi1[i][j] = 0.25*( phi0[i+1][j] + phi0[i-1][j] 
                             + phi0[i][j+1] + phi0[i][j-1] - h*h*S);

            }

        }

        // calcolo l'errore
        Boundary(phi1, x, y, nx, ny, S);
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                deltax = phi1[i+1][j] - 2.0*phi1[i][j] + phi1[i-1][j];
                deltay = phi1[i][j+1] - 2.0*phi1[i][j] + phi1[i][j-1];

                sum += fabs( deltax + deltay - h*h*S );

            }

        }

        err = sum;

        iter += 1;

    }

    return iter;
    
}

// Implemento il metodo di Gauss-Seidel.
// gli do in input gli array bidimensionale delle soluzioni, la griglia,
// la soluzione, la spaziatura della griglia, la tolleranza e
// i punti della griglia.
// Restituisce il numero di iterazioni.
int GaussSeidelMethod(double **phi0, double **phi1, double *x, double *y,
                      double S, double h, double tol, int nx, int ny){

    double sum; // variabile somma in cui metto la soluzione ad ogni step

    double err = 1.e3; // variabile errore che utilizzo per bloccare il ciclo
    double deltax = 0. , deltay = 0. ; // variabili che mi servono per l'errore

    int iter = 0; // variabile per contare le iterazioni

    while( err > tol ){

        sum = 0.0; // azzero la variabile

        // setto le condizioni al bordo
        Boundary(phi0, x, y, nx, ny, S);

        // copio la vecchia soluzione phi0 nella nuova phi1
        for(int i = 0 ; i < nx ; i++ ){

            for( int j = 0 ; j < ny ; j++ ){

                phi0[i][j] = phi1[i][j];

            }

        }

        // implemento il metodo iterativo (ricordando che il primo e ultimo
        // elemento sono fissati dal bordo)
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                phi1[i][j] = 0.25*( phi0[i+1][j] + phi1[i-1][j] 
                             + phi0[i][j+1] + phi1[i][j-1] - h*h*S);

            }

        }

        // calcolo l'errore
        Boundary(phi1, x, y, nx, ny, S);
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                deltax = phi1[i+1][j] - 2.0*phi1[i][j] + phi1[i-1][j];
                deltay = phi1[i][j+1] - 2.0*phi1[i][j] + phi1[i][j-1];

                sum += fabs( deltax + deltay - h*h*S );

            }

        }

        err = sum;

        iter += 1;

    }

    return iter;
    
}

// Implemento il metodo SOR.
// gli do in input gli array bidimensionale delle soluzioni, la griglia,
// la soluzione, la spaziatura della griglia, la tolleranza e
// i punti della griglia.
// Restituisce il numero di iterazioni.
int SORMethod(double **phi0, double **phi1, double *x, double *y,
              double S, double h, double tol, int nx, int ny){

    double sum; // variabile somma in cui metto la soluzione ad ogni step

    double err = 1.e3; // variabile errore che utilizzo per bloccare il ciclo
    double deltax = 0. , deltay = 0. ; // variabili che mi servono per l'errore
    double omega; // parametro di rilassamento

    omega = 2.0/( 1.0 + M_PI/(double)nx );

    int iter = 0; // variabile per contare le iterazioni

    while( err > tol ){

        sum = 0.0; // azzero la variabile

        // setto le condizioni al bordo
        Boundary(phi0, x, y, nx, ny, S);

        // copio la vecchia soluzione phi0 nella nuova phi1
        for(int i = 0 ; i < nx ; i++ ){

            for( int j = 0 ; j < ny ; j++ ){

                phi0[i][j] = phi1[i][j];

            }

        }

        // implemento il metodo iterativo (ricordando che il primo e ultimo
        // elemento sono fissati dal bordo)
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                phi1[i][j] = ( 1 - omega )*phi0[i][j] + 
                             omega*0.25*( phi1[i-1][j] + phi0[i+1][j] + 
                             phi0[i][j+1] + phi1[i][j-1] - h*h*S);

            }

        }

        // calcolo l'errore
        Boundary(phi1, x, y, nx, ny, S);
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                deltax = phi1[i+1][j] - 2.0*phi1[i][j] + phi1[i-1][j];
                deltay = phi1[i][j+1] - 2.0*phi1[i][j] + phi1[i][j-1];

                sum += fabs( deltax + deltay - h*h*S );

            }

        }

        err = sum;

        iter += 1;

    }

    return iter;
    
}