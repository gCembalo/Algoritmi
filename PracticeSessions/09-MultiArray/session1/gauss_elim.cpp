// Implement the Gaussian elimination algorithm for a general linear system
// of n x nequations.
//
// Test the algorithm on 3x3 and 4x4 system
//

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void GaussElimination(double **, double *, double *, int);
void PrintMatrix(double **, int);
void PrintVector(double *, int);

int main(){

    double **A, *b, *x; // definisco matrici e vettori che servono

    // Vediamo il caso di 3x3
    int n;

    cout << "+--------------------------------------+\n" << endl;
    cout << "Vuoi il sistema per le matrici 3x3 o 4x4?" << endl;
    cout << "n = ";
    cin >> n;
    cout << "\n" << endl;

    A = new double *[n]; // array delle righe
    A[0] = new double [n*n];
    b = new double [n];
    x = new double [n];

    // copy static matrix As to dynamic Ad (ciclo nelle slide)
    // stessa cosa per bs e bd
    for( int i = 1 ; i < n ; i++ ){

        A[i] = A[i-1] + n;

    }

    // metto gli elementi di matrice e vettore (a seconda di n)
    if( n == 3 ){

        A[0][0] = 2.;   A[0][1] = -1.;  A[0][2] = 0.;
        A[1][0] = -1.;  A[1][1] = 2.;   A[1][2] = -1.;
        A[2][0] = 0.;   A[2][1] = -1.;  A[2][2] = 2.;

        b[0] = 1.;   
        b[1] = 2.;   
        b[2] = 1.;
    
    }

    if( n == 4 ){
        A[0][0] = 1.;   A[0][1] = 2.;   A[0][2] = 1.;   A[0][3] = -1.;
        A[1][0] = 3.;   A[1][1] = 2.;   A[1][2] = 4.;   A[1][3] = 4.;
        A[2][0] = 4.;   A[2][1] = 4.;   A[2][2] = 3.;   A[2][3] = 4.;
        A[3][0] = 2.;   A[3][1] = 0.;   A[3][2] = 1.;   A[3][3] = 5.; 

        b[0] = 5.;   
        b[1] = 16.;   
        b[2] = 22.;   
        b[3] = 15.;
    }

    // richiamo la funzione di eliminazione gaussiana
    GaussElimination(A, b, x, n);

    cout << endl;
    cout << "Sistema n = " << n << endl;

    PrintVector(x, n); // soluzione del sistema

    delete[] A[0]; // pulisco la memoria
    delete[] A;
    delete[] b;
    delete[] x;

    return 0;

}


// ------------------------- Funzioni implementate ------------------------ //

// implemento la funzione di eliminazione di Gauss.
// gli do in input la matrice da rendere diagonale superiore, il vettore di
// termini noti, il vettore delle incognite e la dimensione
void GaussElimination(double **A, double *b, double *x, int n){

    double tmp; // variabile che mi servirà per la risoluzione
    double g; // varaibile per invertire

    for( int k = 0 ; k < n-1 ; k++ ){   // loop over the Gk's

        for( int i = k+1 ; i < n ; i++ ){   // loop sulle righe

            g = A[i][k] / A[k][k];

            for( int j = k+1 ; j < n ; j++ ){

                A[i][j] -= g*A[k][j];

            }

            A[i][k] = 0.0;
            b[i] -= g*b[k];

        }

    }

    // risoluzione sistema (back-substitution)
    // vedi la formula nelle slide
    for( int i = n-1 ; i >= 0 ; i-- ){

        tmp = b[i];

        for( int j = n-1 ; j > i ; j-- ){

            tmp -= x[j] * A[i][j];

        }

        x[i] = tmp/A[i][i];

    }
    
}

// implemento la funzione per stampare una matrice (dinamica) di n dimensioni
void PrintMatrix(double **A, int n){

    cout << fixed << setprecision(4);

    for( int i = 0 ; i < n ; i++ ){

        for(int j = 0 ; j < n ; j++ ){

            cout << setw(10) << right << A[i][j] << " ";

        }

        cout << endl;

    }

}

// implemento la funzione per stampare un vettore (dinamico) di n dimensioni
void PrintVector(double *v, int n){

    cout << fixed << setprecision(4);

    for(int j = 0 ; j < n ; j++ ){

        cout << setw(10) << right << v[j] << endl;
        
    }

}