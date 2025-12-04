// implement the tridiagonal matrix solver and test it on the following system
// (nelle slide)
//
// Solution: {2,-3,4,-2,1}.
//

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void TridiagSolver(double *, double *, double *, double *, double *, int);
void PrintVector(double *, int);

int main(){

    // non definisco la matrice ma direttamente i vettori (diagonali) e
    // x e b
    int n = 5; // dimensione matrici e vettori
    double *a, *b, *c, *r, *x;

    a = new double [n];
    b = new double [n];
    c = new double [n];
    r = new double [n];
    x = new double [n];

    // metto gli elementi sotto, sopra e dei termini noti
    // i termini a[0] e c[4] non verranno utilizzati (vedi eq. slide)
    a[0] = 0.;   a[1] = 1.;   a[2] = 1.;   a[3] = 1.;   a[4] = 1.;
    b[0] = 2.;   b[1] = 2.;   b[2] = 2.;   b[3] = 2.;   b[4] = 2.;
    c[0] = 1.;   c[1] = 1.;   c[2] = 1.;   c[3] = 1.;   c[4] = 0.;
    r[0] = 1.;   r[1] = 0.;   r[2] = 3.;   r[3] = 1.;   r[4] = 0.;
    // commentato in TridiagSolver c'è un ciclo per una matrice generica
    // per il calcolo di r

    // richiamo il metodo e stampo la soluzione
    TridiagSolver(a, b, x, r, c, n);
    PrintVector(x, n);

    // pulisco
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] r;
    delete[] x;

    return 0;

}


// ------------------------ Funzioni implementate ------------------------ //

// Implemento la funzione di TridiagSolver
// gli do in input i vettori contenenti rispettivamente: gli elementi sotto
// la diagonale, i termini noti, le soluzioni, r, gli elementi sopra la diagonale
// e la dimensione
void TridiagSolver(double *a, double *b, double *x, double *r, double *c, int n){

    // definisco i vettori h e p
    double *h = new double [n];
    double *p = new double [n];

    // calcolo gli elementi h[n] e p[n] per il metodo Tridiag
    // separo i termini patologici
    h[0] = c[0]/b[0];
    p[0] = r[0]/b[0];

    for( int i = 1 ; i < n ; i++ ){

        h[i] = c[i] / ( b[i] - a[i]*h[i-1] );

        p[i] = ( r[i] - a[i]*p[i-1] ) / ( b[i] - a[i]*h[i-1] );

    }

    // applico il metodo di risoluzione
    x[n-1] = p[n-1]; // termine patologico non avendo definito x[n+1]
    
    // applichiamo back-substitution
    for( int i = n-1 ; i >= 0 ; i-- ){

        x[i] = p[i] - h[i]*x[i+1];

    }

    // pulisco
    delete[] h;
    delete[] p;

}

// implemento la funzione per stampare un vettore (dinamico) di n dimensioni
void PrintVector(double *v, int n){

    cout << fixed << setprecision(4);

    for(int j = 0 ; j < n ; j++ ){

        cout << setw(10) << right << v[j] << endl;
        
    }

}