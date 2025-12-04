// Modify the previous program by introducing partial pivoting. 
// The modification should search for the row jmax such that
//               |A[jmax][k]| >= |A[j][k]|for j > k.
//

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void Pivot(double **, double *, double *, int);
void SwapRows(double **, double *, int, int, int);
void PrintMatrix(double **, int);
void PrintVector(double *, int);
void MVmult_dinamic(double **, double *, double *, int, int);

int main(){

    double **A, *b, *x; // definisco matrici e vettori che servono
    int n = 4; // dimensione matrici e vettori

    A = new double *[n]; // array delle righe
    A[0] = new double [n*n];
    b = new double [n];
    x = new double [n];

    // copy static matrix As to dynamic Ad (ciclo nelle slide)
    for( int i = 1 ; i < n ; i++ ){

        A[i] = A[i-1] + n;

    }

    // stessa matrice di "gauss_elim.cpp" ma con l'elemento [1][1] cambiato
    A[0][0] = 1.;   A[0][1] = 2.;   A[0][2] = 1.;   A[0][3] = -1.;
    A[1][0] = 3.;   A[1][1] = 6.;   A[1][2] = 4.;   A[1][3] = 4.;
    A[2][0] = 4.;   A[2][1] = 4.;   A[2][2] = 3.;   A[2][3] = 4.;
    A[3][0] = 2.;   A[3][1] = 0.;   A[3][2] = 1.;   A[3][3] = 5.; 

    b[0] = 5.;   
    b[1] = 16.;   
    b[2] = 22.;   
    b[3] = 15.;

    // salvo la matrice originale per poi controllare la soluzione
    double **A0;
    A0 = new double *[n]; // array delle righe
    A0[0] = new double [n*n];

    for( int i = 1 ; i < n ; i++ ){

        A0[i] = A0[i-1] + n;

    }

    for( int i = 0 ; i < n ; i++ ){

        for( int j = 0 ; j < n ; j++){

            A0[i][j] = A[i][j];

        }

    }

    // richiamo la funzione di Pivot
    Pivot(A, b, x, n);

    cout << endl;
    cout << "x = " << endl;
    PrintVector(x, n); // soluzione del sistema
    cout << endl;

    cout << "Controllo Ax = b :" << endl;
    MVmult_dinamic(A0, x, b, n, n);
    PrintVector(b, n);

    delete[] A[0]; // pulisco la memoria
    delete[] A;
    delete[] A0[0]; // pulisco la memoria
    delete[] A0;
    delete[] b;
    delete[] x;

    return 0;

}


// ------------------------- Funzioni implementate ------------------------ //

// implemento la funzione pivot (evoluzione della GaussElimination)
// gli do in input la matrice, il vettore dei termini noti, il vettore
// soluzione e la dimensione
void Pivot(double **A, double *b, double *x, int n){

    double tmp , g; // variabili temporanea e per eseguire l'eliminazione
    // variabili per salvare elementi di matrice (diagonale e non)
    double Akk , Aik;
    int m;

    for( int k = 0 ; k < n-1 ; k++ ){   // loop over the Gk's

        Akk = fabs( A[k][k] ); // salvo l'elemento diagonale
        m = k; // salvo l'indice

        // ricerco la riga con |A_{ik}| massimo (e maggiore) di |A_{kk}|
        for( int i = k+1 ; i < n ; i++ ){ // loop sulle righe

            Aik = fabs( A[i][k] ); // salvo l'elemento

            if( Aik > Akk ){

                Akk = Aik;
                m = i;

            }

        }

        // scambio le righe solo se i due elementi non coincidono
        if( k != m ){

            SwapRows(A, b, m, k, n);
        
        }

        // implemento l'eliminazione gaussiana (uguale all'altra funzione)
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

// implemento la funzione per scambiare righe
// gli do in input la matrice, il vettore, le due righe da scambiare e la
// dimensione
void SwapRows(double **A, double *b, int i1, int i2, int n){

    double tmp; // variabile temporanea per memorizzare gli elementi

    // scambio elementi della matrice
    for( int j = 0 ; j < n ; j++ ){

        tmp = A[i1][j];

        A[i1][j] = A[i2][j];

        A[i2][j] = tmp;

    }

    // scambio elementi del vettore
    tmp = b[i1];

    b[i1] = b[i2];

    b[i2] = tmp;

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

// implemento la funzione per la moltiplicazione matrice vettore in modo
// DINAMICO.
// do in input la matrice, il vettore, il vettore prodotto, e il numero di 
// righe e colonne
void MVmult_dinamic(double **Ad, double *bd, double *Abd, int nrow, int ncol){

    for( int i = 0 ; i < nrow ; i++ ){

        Abd[i] = 0.0; // inizializzo il valore di Ab

        for( int j = 0 ; j < ncol ; j++ ){

            Abd[i] += Ad[i][j] * bd[j];

        }

    }
    
}