// solve the Boundary Value Problem using finite diﬀerence scheme:
//      \dv[2]{y}{x} = 1
//      y(0) = 0       ,        y(1) = 0
// Use 32 points (31 intervals) in the range 0 ≤ x ≤ 1.
//
// Now change your b.c. and solve
//      \dv[2]{y}{x} = 1
//      y(0) = 1       ,        y(1) = 0.9
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double RHSFuncBVP(double);
void TridiagSolver(double *, double *, double *, double *, double *, int);
void PrintVector(double *, int);

int main(){

    ofstream fdata;
    fdata.open("bvp.dat"); // file per le soluzioni
    
    int n = 32; // dimensione matrice = numero di punti

    // definisco i vettori degli elementi sopra e sotto la diagonale
    double *a, *b, *c, *r, *y;

    a = new double [n];
    b = new double [n];
    c = new double [n];
    r = new double [n];
    y = new double [n];

    // setto le condizioni iniziali
    double x, h; // variabile e incremento derivata
    double xi = 0., xf = 1.;
    double alpha, beta;

    int bvp; // variabile per scegliere qualc condizione al bordo mettere
    cout << "Quale condizioni al contorno vuoi?\n"
         << "n = 1 :    y(0) = 0    ,   y(1) = 0\n"
         << "n = 2 :    y(0) = 1    ,   y(1) = 0.9\n"
         << "n = ";
    cin >> bvp;

    if( bvp == 1 ){
        alpha = 0.;
        beta = 0.;
    }

    if( bvp == 2 ){
        alpha = 1.;
        beta = 0.9;
    }

    h = fabs( xf - xi )/(double)( n - 1 );

    y[0] = alpha;
    y[n-1] = beta;

    // riempio i vettori a, b, c ed r
    for( int i = 0 ; i < n ; i++ ){

        x = xi + i*h;

        a[i] = 1.; // sotto la diagonale (y_{i-1})
        b[i] = -2.; // diagonale (y_{i})
        c[i] = 1.; // sopra la diagonale (y_{i+1})

        r[i] = h*h*RHSFuncBVP(x);

    }

    // metto le condizioni iniziali nel primo e ultimo elemento di r
    r[1] -= y[0];
    r[n-2] -= y[n-1];

    // richiamo il TridiagSolver per trovare il vettore delle soluzioni
    // schiftiamo di 1 l'argomento, poiché l'indice 0 lo fissiamo con la
    // condizione al bordo e quindi il primo elemento sarà y_2
    TridiagSolver(a+1, b+1, y+1, r+1, c+1, n-2);
    // PrintVector(y, n);

    // stampo il vettore nel file
    for ( int i = 0 ; i < n ; i++ ){

        x = xi+i*h;

        fdata << x << " " << y[i] << endl;

    }

    fdata.close ();

    // pulisco
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] r;
    delete[] y;

    return 0;

}


// ------------------------ Funzioni implementate ------------------------ //

// definisco la RHS function della ODE
double RHSFuncBVP(double x){

    return 1.;

}

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