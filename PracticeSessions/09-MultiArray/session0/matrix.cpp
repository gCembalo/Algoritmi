// write a simple program that creates a NxN matrix A using either standard
// native or dynamical method (use, e.g. N = 4).
// Print the matrix on screen using nice-formatted output. 
//

#include <iostream>
#include <iomanip>
#include <cmath>

// definizione di variabili globali statiche
#define Nrow 4
#define Ncol 4

using namespace std;

void MVmult_static(double [Nrow][Ncol], double *, double *, int, int);
void MVmult_dinamic(double **, double *, double *, int, int);
void PrintMatrix(double **, int);
void PrintVector(double *, int);

int main(){

    cout << fixed << setprecision(4);

    // definisco le variabili per righe e colonne
    int nrow = Nrow , ncol = Ncol;

    // definisco matrici e vettori
    double As[Nrow][Ncol]; // definizione matrice 4x4 per caso STATICO
    double bs[Nrow]; // vettore caso STATICO
    double Abs[Nrow]; // vettore moltiplicazione Ab caso STATICO
    double **Ad , *bd , *Abd; // definizione puntatore per caso DINAMICO

    cout << "+--------------------------------------+\n" << endl;
    cout << "Allocazione statica:\n" << endl;

    // assegno i valori alla matrice
    As[0][0] = 1.0;   As[0][1] = 3.0;   As[0][2] = 2.0;   As[0][3] = -4.0;
    As[1][0] = 7.0;   As[1][1] = 2.0;   As[1][2] = 4.0;   As[1][3] = 1.0;
    As[2][0] = 0.0;   As[2][1] = -1.0;  As[2][2] = 2.0;   As[2][3] = 2.0;
    As[3][0] = 6.0;   As[3][1] = 3.0;   As[3][2] = 0.0;   As[3][3] = 1.0; 

    // assegno i valori al vettore
    bs[0] = 1.0;
    bs[1] = 0.0;
    bs[2] = 3.0;
    bs[3] = 2.0;

    // stampo a video la matrice (codice nelle slide)
    cout << "matrice A" << endl;
    for( int i = 0 ; i < nrow ; i++){

        for( int j = 0 ; j < ncol ; j++){

            cout << setw(10) << right << As[i][j] << " ";

        }

        cout << endl;

    }

    // stampo a video il vettore
    cout << "\nvettore b" << endl;
    for( int j = 0 ; j < nrow ; j++ ){

        cout << setw(10) << right << bs[j] << endl;

    }

    // richiamo la funzione di moltiplicazione STATICA e stampo Ab
    MVmult_static(As, bs, Abs, nrow, ncol);

    cout << "\nvettore Ab" << endl;
    for( int j = 0 ; j < nrow ; j++ ){

        cout << setw(10) << right << Abs[j] << endl;

    }


    cout << "\n+--------------------------------------+\n" << endl;
    cout << "Allocazione dinamica:\n" << endl;

    // creo il multi array A[nrow][nrow]
    Ad = new double *[nrow]; // array delle righe
    Ad[0] = new double [ncol*nrow];
    bd = new double [nrow];
    Abd = new double [nrow];

    // copy static matrix As to dynamic Ad (ciclo nelle slide)
    // stessa cosa per bs e bd
    for( int i = 1 ; i < nrow ; i++ ){

        Ad[i] = Ad[i-1] + ncol;

    }

    // inserisco i vari elementi nella matrice e nel vettore
    for( int i = 0 ; i < nrow ; i++ ){

        for( int j = 0 ; j < ncol ; j++ ){

            Ad[i][j] = As[i][j];

        }

    }

    for( int j = 0 ; j < ncol ; j++ ){

        bd[j] = bs[j];

    }

    // stampo a video
    cout << "matrice A" << endl;
    PrintMatrix(Ad, Nrow);

    // stampo a video il vettore
    cout << "\nvettore b" << endl;
    PrintVector(bd, Nrow);

    // richiamo la funzione di moltiplicazione Dinamica e stampo Ab
    MVmult_dinamic(Ad, bd, Abd, nrow, ncol);

    cout << "\nvettore Ab" << endl;
    PrintVector(Abd, Nrow);

    // libero la memoria
    delete[] Ad[0];
    delete[] Ad;
    delete[] bd;
    delete[] Abd;

    return 0;

}


// ------------------------ Funzioni implementate ------------------------ //

// implemento la funzione per la moltiplicazione matrice vettore in modo
// STATICO.
// do in input la matrice, il vettore, il vettore prodotto, e il numero di 
// righe e colonne
void MVmult_static(double As[Nrow][Ncol], double* bs, double* Abs, 
                   int nrow, int ncol){
    
    for( int i = 0 ; i < nrow ; i++ ){

        Abs[i] = 0.0; // inizializzo il valore di Ab

        for( int j = 0 ; j < ncol ; j++ ){

            Abs[i] += As[i][j] * bs[j];

        }

    }

}


void MVmult_dinamic(double **Ad, double *bd, double *Abd, int nrow, int ncol){

    for( int i = 0 ; i < nrow ; i++ ){

        Abd[i] = 0.0; // inizializzo il valore di Ab

        for( int j = 0 ; j < ncol ; j++ ){

            Abd[i] += Ad[i][j] * bd[j];

        }

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