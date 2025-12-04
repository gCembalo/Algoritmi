// solve the Poisson equation
//      \pdv[2]{\phi}{x} + \pdv{\phi}{y} - S(x,y) = 0
// on the unit square 0 ≤ x,y ≤ 1with S = const. and b.c. given by 
// the exact solution
// \phi(x,y) = e^{-\pi x}\sin{(-\pi y)} + S/4 (x^2 + y^2)
//
// Set NX = NY = 32and try S = 0and then S = 2using Jacobi, Gauss-Seidel and
// SOR. Compare the number of iterations necessary to achieve convergence,
//  using the residual and a tolerance of 10-7. Results are given by the
// following table.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double SolElliptic(double, double, double);
void Boundary(double **, double *, double *, int, double);

int main(){

    ofstream fdata;
    fdata.open("elliptic.dat"); // file per le soluzioni

    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 2 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 10 );

    int n = 32; // numero di nodi della griglia ( nx = ny )

    // definisco il dominio di integrazione
    double xi = 0.0 , yi = 0.0;
    double xf = 1.0 , yf = 1.0;

    double tol = 1.e-7; // tolleranza

    double S; // RHSfunction

    double h = fabs(xf-xi)/(double)(n-1); // intervallo

    // definisco gli array 2d
    double **phi0 , **phi1;
    
    for ( int i=1; i<n; i++ )
    {
        phi0[i] = phi0[i-1]+n;
    }

    phi1 = new double*[n];
    phi1[0] = new double[n*n];

    for ( int i=1; i<n; i++ )
    {
        phi1[i] = phi1[i-1]+n;
    }


    return 0;

}


// ------------------------ Funzioni implementate ------------------------ //

// Implemento la funzione soluzione
double SolElliptic(double x, double y, double S){

    return exp( -M_PI*x )*sin( -M_PI*y ) + S*( x*x + y*y )*0.25;

}

// Implemento una funzione che mi fissa le condizioni al contorno (alla griglia)
// Gli do in input il vettore soluzione, gli array x e y della griglia, il
// numero di punti ed il valore di S
void Boundary(double **phi, double *x, double *y, int n, double S){

    // Costruisco il ciclo per fissare le condizioni (usando la soluzione)
    for( int i = 0 ; i < n ; i++ ){

        phi[0][i] = SolElliptic(x[0], y[i], S);
        phi[i][0] = SolElliptic(x[i], y[0], S);
        phi[n-1][i] = SolElliptic(x[n-1], y[i], S);
        phi[i][n-1] = SolElliptic(x[i], y[n-1], S);
        
    }

}