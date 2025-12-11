// Find the steady-state temperature distribution of a rectangular plate
// 0 ≤x ≤ 2, 0 ≤ y ≤1, insulated at x=0, with temperature fixed to 0 and
// to 2-x on the lower and upper sides and constant heat flux = 3 at
// the right side. For constant thermal conductivity, this entails the
// solution of the Laplace equation:
//          \nabla^2 \phi = 0
// with Neumann boundary conditions.
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
    fdata.open("conductivity.dat"); // file per le soluzioni

    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 2 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 10 );

    int nx , ny; // numero di nodi della griglia
    nx = 129;
    ny = 65;

    double tol = 1.e-7; // tolleranza

    // definisco le variabili per le iterazioni
    int iterJ , iterGS , iterSOR;

    // Le condizioni sulla sorgente le metto con una funzione apposita
    double S = 0.0; // Sorgente (inizializzata)
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
    double xi = 0.0 , yi = 0.0;
    double xf = 2.0 , yf = 1.0;
    double x[nx], y[ny];

    double h = fabs(xf-xi)/(double)(nx-1); // intervallo

    // inizializzo la griglia
    for( int i = 0 ; i < nx ; i++ ){
        x[i] = xi + i*h;
    }
    for( int j = 0 ; j < ny ; j++ ){
        y[j] = yi + j*h;
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
        
    for( int i = 0 ; i < nx ; i++ ){

        for( int j = 0 ; j < ny ; j++ ){

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
        
    for( int i = 0 ; i < nx ; i++ ){

        for( int j = 0 ; j < ny ; j++ ){

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
        
    for( int i = 0 ; i < nx ; i++ ){

        for( int j = 0 ; j < ny ; j++ ){

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

    // Costruisco il ciclo per fissare le condizioni
    double sigma;
    int i,j;

    double dx = x[1] - x[0];

    // condizioni bottom: Dirichlet b.c.
    j = 0;
    for( i = 0 ; i < nx ; i++ ){
        phi[i][j] = 0.0;
    }

    // condizioni top: Dirichlet b.c.
    j = ny - 1;
    for( i = 0 ; i < nx ; i++ ){
        phi[i][j] = 2.0 - x[i];
    }

    // condizioni left: Neumann b.c.
    sigma = 0.0;
    i = 0;
    for( j = 0 ; j < ny ; j++ ){

        //phi[i][j] = phi[i+1][j] - dx*sigma; // primo ordine

        // secondo ordine
        phi[i][j] = ( -phi[i+2][j] + 4*phi[i+1][j] - 2.0*dx*sigma )/3.0;

    }

    // condizioni right: Neumann b.c.
    sigma = 3.0;
    i = nx - 1;
    for( j = 0 ; j < ny ; j++ ){

        //phi[i][j] = phi[i-1][j] - dx*sigma; // primo ordine

        // secondo ordine
        phi[i][j] = ( -phi[i-2][j] + 4*phi[i-1][j] - 2.0*dx*sigma )/3.0;

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
    double diff; // variabile differenza in cui metto la soluzione ad ogni step

    double err = 1.e3; // variabile errore che utilizzo per bloccare il ciclo
    double deltax = 0. , deltay = 0. ; // variabili che mi servono per l'errore

    int iter = 0; // variabile per contare le iterazioni

    while( err > tol ){

        //sum = 0.0; // azzero la variabile
        diff = 0.0;

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
        err = 0.0;

        Boundary(phi1, x, y, nx, ny);
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                S = Source(x[i], y[j]);

                diff = fabs( ( phi1[i][j] - phi0[i][j] ) * h * h );

                //deltax = phi1[i+1][j] - 2.0*phi1[i][j] + phi1[i-1][j];
                //deltay = phi1[i][j+1] - 2.0*phi1[i][j] + phi1[i][j-1];
                //sum += fabs( deltax + deltay - h*h*S );

                err += diff;

            }

        }

        //err = sum;

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

    double phi0[nx][ny]; // matrice soluzioni step precedente (serve per l'errore)
    // inizializzo le soluzioni
    for( int i = 1 ; i < nx-1 ; i++ ){
        // inizializzo i vettori
        for( int j = 1 ; j < ny-1 ; j++ ){
            phi0[i][j] = 0.0;
        }
    }

    double sum; // variabile somma in cui metto la soluzione ad ogni step
    double diff; // variabile differenza in cui metto la soluzione ad ogni step

    double err = 1.e3; // variabile errore che utilizzo per bloccare il ciclo
    double deltax = 0. , deltay = 0. ; // variabili che mi servono per l'errore

    int iter = 0; // variabile per contare le iterazioni

    while( err > tol ){

        //sum = 0.0; // azzero la variabile
        diff = 0.0;

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
        err = 0.0;

        Boundary(phi, x, y, nx, ny);
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                S = Source(x[i], y[j]);

                diff = fabs( ( phi[i][j] - phi0[i][j] ) * h * h );

                //deltax = phi1[i+1][j] - 2.0*phi1[i][j] + phi1[i-1][j];
                //deltay = phi1[i][j+1] - 2.0*phi1[i][j] + phi1[i][j-1];
                //sum += fabs( deltax + deltay - h*h*S );

                err += diff;

            }

        }

        //err = sum;

        iter += 1;

        // copio la soluzione
        for(int i=0;i<nx;i++){

            for(int j=0;j<ny;j++){

                phi0[i][j] = phi[i][j];

            }

        }

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

    double phi0[nx][ny]; // matrice soluzioni step precedente (serve per l'errore)
    // inizializzo le soluzioni
    for( int i = 1 ; i < nx-1 ; i++ ){
        // inizializzo i vettori
        for( int j = 1 ; j < ny-1 ; j++ ){
            phi0[i][j] = 0.;
        }
    }

    double sum; // variabile somma in cui metto la soluzione ad ogni step
    double diff; // variabile differenza in cui metto la soluzione ad ogni step

    double err = 1.e3; // variabile errore che utilizzo per bloccare il ciclo
    double deltax = 0. , deltay = 0. ; // variabili che mi servono per l'errore
    double omega; // parametro di rilassamento

    omega = 2.0/( 1.0 + M_PI/(double)nx );

    int iter = 0; // variabile per contare le iterazioni

    while( err > tol ){

        //sum = 0.0; // azzero la variabile
        diff = 0.0;

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
        err = 0.0;

        Boundary(phi, x, y, nx, ny);
        for( int i = 1 ; i < nx-1 ; i++ ){

            for( int j = 1 ; j < ny-1 ; j++ ){

                S = Source(x[i], y[j]);

                diff = fabs( ( phi[i][j] - phi0[i][j] ) * h * h );

                //deltax = phi1[i+1][j] - 2.0*phi1[i][j] + phi1[i-1][j];
                //deltay = phi1[i][j+1] - 2.0*phi1[i][j] + phi1[i][j-1];
                //sum += fabs( deltax + deltay - h*h*S );

                err += diff;

            }

        }

        //err = sum;

        iter += 1;

        // copio la soluzione
        for(int i=0;i<nx;i++){

            for(int j=0;j<ny;j++){

                phi0[i][j] = phi[i][j];

            }

        }

    }

    return iter;
    
}