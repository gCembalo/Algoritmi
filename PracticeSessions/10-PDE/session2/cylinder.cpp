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

double Source(double, double);
double SolCylinder(double, double, double);
void Boundary(double **, double *, double *, int, int);
int JacobiMethod(double **, double **, double *, double *, double, double, 
                    double, int, int);
int GaussSeidelMethod(double **, double *, double *, double, double, 
                    double, int, int);
int SORMethod(double **, double *, double *, double, double, 
                    double, int, int);

int main(){

    ofstream fdata;
    fdata.open("cylinder.dat"); // file per le soluzioni

    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 2 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 10 );

    int nx , ny , n = 128; // numero di nodi della griglia ( nx = ny )
    nx = ny = n;

    double tol = 1.e-7; // tolleranza

    // definisco le variabili per le iterazioni
    int iterJ , iterGS , iterSOR;

    // Le condizioni sulla sorgente le metto con una funzione apposita
    double S; // Sorgente
    double r; // variabile raggio
    double a = 0.1;

    // definisco gli array 2d
    double **phi0 = new double *[nx];
    double **phi1 = new double *[nx];

    phi0[0] = new double [nx*ny];
    phi1[0] = new double [nx*ny];
    
    for( int i=1 ; i<nx ; i++ ){

        phi0[i] = phi0[i-1] + ny;
        phi1[i] = phi1[i-1] + ny;

    }

    // definisco il dominio di integrazione
    double xi = -1.0 , yi = -1.0;
    double xf = 1.0 , yf = 1.0;
    double x[nx], y[ny];

    double h = fabs(xf-xi)/(double)(nx-1); // intervallo

    // inizializzo la griglia
    for( int i = 0 ; i < nx ; i++ ){
        x[i] = xi + i*h;
        y[i] = yi + i*h;
    }

    // inizializzo le soluzioni
    for( int i = 1 ; i < nx-1 ; i++ ){
            
        // inizializzo i vettori
        for( int j = 1 ; j < ny-1 ; j++ ){

            phi0[i][j] = phi1[i][j] = 0.;

        }

    }

    // calcolo le iterazioni invocando il metodo
    iterJ = 0;
    iterJ = JacobiMethod(phi0, phi1, x, y, S, h, tol, nx, ny);
        
    for( int i = 0 ; i < n ; i++ ){

        for( int j = 0 ; j < n ; j++ ){

            fdata << x[i] << " " << y[j] << " " << phi1[i][j] << endl;

        }
            
        fdata << endl;

    }

    cout << "Jacobi: \t k = " << iterJ << endl;

    // costruisco un ciclo in cui uso il metodo di GS
    for( int i = 1 ; i < nx-1 ; i++ ){
            
        // inizializzo i vettori
        for( int j = 1 ; j < ny-1 ; j++ ){

            phi0[i][j] = phi1[i][j] = 0.;

        }

    }

    // calcolo le iterazioni invocando il metodo
    iterGS = 0;
    iterGS = GaussSeidelMethod(phi1, x, y, S, h, tol, nx, ny);
        
    for( int i = 0 ; i < n ; i++ ){

        for( int j = 0 ; j < n ; j++ ){

            fdata << x[i] << " " << y[j] << " " << phi1[i][j] << endl;

        }
            
        fdata << endl;

    }

    cout << "Gauss Seidel: \t k = " << iterGS << endl;

    // costruisco un ciclo in cui uso il metodo di SOR
    for( int i = 1 ; i < nx-1 ; i++ ){
            
        // inizializzo i vettori
        for( int j = 1 ; j < ny-1 ; j++ ){

            phi0[i][j] = phi1[i][j] = 0.;

        }

    }

    // calcolo le iterazioni invocando il metodo
    iterSOR = 0;
    iterSOR = SORMethod(phi1, x, y, S, h, tol, nx, ny);
        
    for( int i = 0 ; i < n ; i++ ){

        for( int j = 0 ; j < n ; j++ ){

            fdata << x[i] << " " << y[j] << " " << phi1[i][j] << endl;

        }
            
        fdata << endl;

    }

    cout << "SOR: \t k = " << iterSOR << endl;

    fdata.close();

    delete[] phi0[0];
    delete[] phi1[0];
    delete[] phi0;
    delete[] phi1;

    return 0;

}


// ------------------------ Funzioni implementate ------------------------ //

// Implemento la funzione sorgente
double Source(double x, double y){

    double a = 0.1;
    double rho0 = 1.0;

    double r = sqrt( x*x + y*y ); // coordinata radiale

    if( r >= 0 && r <= a ){
        return rho0;
    }
    if( r > a){
        return 0;
    }

    return 0.0;

}


// Implemento la funzione soluzione
double SolCylinder(double x, double y){

    double rho0 = 1.0;
    double a = 0.1;

    double r = sqrt( x*x + y*y ); // coordinata radiale

    if( r >= 0 && r <= a ){
        return -rho0*r*r*0.25;
    }
    if( r > a){
        return -rho0*a*a*0.5*( log(r/a) + 0.5 );
    }

    return 0.0;

}

// Implemento una funzione che mi fissa le condizioni al contorno (alla griglia)
// Gli do in input il vettore soluzione, gli array x e y della griglia, il
// numero di punti
void Boundary(double **phi, double *x, double *y, int nx, int ny){

    // Costruisco il ciclo per fissare le condizioni (usando la soluzione)

    // bordo sinistro e destro
    for (int i = 0; i < nx; i++) {
        phi[i][0]     = SolCylinder(x[i], y[0]);
        phi[i][ny-1]  = SolCylinder(x[i], y[ny-1]);
    }

    // bordo inferiore e superiore
    for (int j = 0; j < ny; j++) {
        phi[0][j]     = SolCylinder(x[0], y[j]);
        phi[nx-1][j]  = SolCylinder(x[nx-1], y[j]);
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
        Boundary(phi0, x, y, nx, ny);

        // implemento il metodo iterativo (ricordando che il primo e ultimo
        // elemento sono fissati dal bordo)
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                S = Source(x[i],y[j]);

                phi1[i][j] = 0.25*( phi0[i+1][j] + phi0[i-1][j] 
                             + phi0[i][j+1] + phi0[i][j-1] - h*h*S);

            }

        }

        // calcolo l'errore
        Boundary(phi1, x, y, nx, ny);
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                S = Source(x[i], y[j]);

                deltax = phi1[i+1][j] - 2.0*phi1[i][j] + phi1[i-1][j];
                deltay = phi1[i][j+1] - 2.0*phi1[i][j] + phi1[i][j-1];

                sum += fabs( deltax + deltay - h*h*S );

            }

        }

        err = sum;

        iter += 1;

        // copio la soluzione
        for(int i=0;i<nx;i++){

            for(int j=0;j<ny;j++){

                phi0[i][j] = phi1[i][j];

            }

        }

    }

    return iter;
    
}

// Implemento il metodo di Gauss-Seidel.
// gli do in input gli array bidimensionale delle soluzioni, la griglia,
// la soluzione, la spaziatura della griglia, la tolleranza e
// i punti della griglia.
// Restituisce il numero di iterazioni.
int GaussSeidelMethod(double **phi, double *x, double *y,
                      double S, double h, double tol, int nx, int ny){

    double sum; // variabile somma in cui metto la soluzione ad ogni step

    double err = 1.e3; // variabile errore che utilizzo per bloccare il ciclo
    double deltax = 0. , deltay = 0. ; // variabili che mi servono per l'errore

    int iter = 0; // variabile per contare le iterazioni

    while( err > tol ){

        sum = 0.0; // azzero la variabile

        // setto le condizioni al bordo
        Boundary(phi, x, y, nx, ny);

        // implemento il metodo iterativo (ricordando che il primo e ultimo
        // elemento sono fissati dal bordo)
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                S = Source(x[i],y[j]);

                phi[i][j] = 0.25*( phi[i+1][j] + phi[i-1][j] 
                             + phi[i][j+1] + phi[i][j-1] - h*h*S);

            }

        }

        // calcolo l'errore
        Boundary(phi, x, y, nx, ny);
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                S = Source(x[i], y[j]);

                deltax = phi[i+1][j] - 2.0*phi[i][j] + phi[i-1][j];
                deltay = phi[i][j+1] - 2.0*phi[i][j] + phi[i][j-1];

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
int SORMethod(double **phi, double *x, double *y,
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
        Boundary(phi, x, y, nx, ny);

        // implemento il metodo iterativo (ricordando che il primo e ultimo
        // elemento sono fissati dal bordo)
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                S = Source(x[i], y[j]);
                
                phi[i][j] = ( 1 - omega )*phi[i][j] + 
                             omega*0.25*( phi[i-1][j] + phi[i+1][j] + 
                             phi[i][j+1] + phi[i][j-1] - h*h*S);

            }

        }

        // calcolo l'errore
        Boundary(phi, x, y, nx, ny);
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                S = Source(x[i], y[j]);

                deltax = phi[i+1][j] - 2.0*phi[i][j] + phi[i-1][j];
                deltay = phi[i][j+1] - 2.0*phi[i][j] + phi[i][j-1];

                sum += fabs( deltax + deltay - h*h*S );

            }

        }

        err = sum;

        iter += 1;

    }

    return iter;
    
}