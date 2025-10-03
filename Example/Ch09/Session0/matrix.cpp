#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

#define Nrow 4
#define Ncol 4

using namespace std;

void MVmult_static ( double[Nrow][Ncol], double*, double*, int, int );
void MVmult_dinamic ( double**, double*, double*, int, int );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << fixed << setprecision ( 4 );

    int i, j;
    int nrow = Nrow, ncol = Ncol;

    double As[Nrow][Ncol], bs[ncol], Abs[ncol];
    double **Ad, *bd, *Abd;

    cout << endl;
    cout << "ALLOCAZIONE STATICA" << endl;
    cout << "----------------------------------------" << endl;
    cout << "matrice A" << endl;

    As[0][0] = 1.;   As[0][1] = 3.;   As[0][2] = 2.;   As[0][3] = -4.;
    As[1][0] = 7.;   As[1][1] = 2.;   As[1][2] = 4.;   As[1][3] = 1.;
    As[2][0] = 0.;   As[2][1] = -1.;  As[2][2] = 2.;   As[2][3] = 2.;
    As[3][0] = 6.;   As[3][1] = 3.;   As[3][2] = 0.;   As[3][3] = 1.; 

    for ( i=0; i<nrow; i++ )
    {
        for ( j=0; j<ncol; j++ )
        {
            cout << setw(10) << right << As[i][j];
        }

        cout << endl;
    }

    cout << "----------------------------------------" << endl;
    cout << "vettore b" << endl;

    bs[0] = 1.;   
    bs[1] = 0.;   
    bs[2] = 3.;   
    bs[3] = 2.;  

    for ( j=0; j<ncol; j++ )
    {
        cout << setw(10) << right << bs[j] << endl;
    }
    
    MVmult_static ( As, bs, Abs, nrow, ncol );

    cout << "----------------------------------------" << endl;
    cout << "vettore Ab" << endl;

    for ( i=0; i<nrow; i++ )
    {
        cout << setw(10) << right << Abs[i] << endl;
    }
    
    cout << endl;
    cout << "ALLOCAZIONE DINAMICA" << endl;

    Ad = new double *[nrow];
    Ad[0] = new double [ncol*nrow];
    bd = new double [ncol];
    Abd = new double [nrow];

    for ( i=1; i<nrow; i++ )
    {
        Ad[i] = Ad[i-1]+ncol;
    }

    for ( i=0; i<nrow; i++ )
    {
        for ( j=0; j<ncol; j++ )
        {
            Ad[i][j] = As[i][j];
        }
    }

    for ( j=0; j<ncol; j++ )
    {
        bd[j] = bs[j];
    }

    cout << "----------------------------------------" << endl;
    cout << "matrice A" << endl;

    for ( i=0; i<nrow; i++ )
    {
        for ( j=0; j<ncol; j++ )
        {
            cout << setw(10) << right << Ad[i][j];
        }

        cout << endl;
    }

    cout << "----------------------------------------" << endl;
    cout << "vettore b" << endl;

    for ( j=0; j<ncol; j++ )
    {
        cout << setw(10) << right << bd[j] << endl;
    }
    
    MVmult_dinamic ( Ad, bd, Abd, nrow, ncol );

    cout << "----------------------------------------" << endl;
    cout << "vettore Ab" << endl;

    for ( i=0; i<nrow; i++ )
    {
        cout << setw(10) << right << Abd[i] << endl;
    }
    
    cout << endl;

    delete[] Ad[0];
    delete[] Ad;
    delete[] bd;
    delete[] Abd;

    return 0;
}

void MVmult_static ( double As[Nrow][Ncol], double* bs, double* Abs, int nrow, int ncol )
{
    int i, j;

    for ( i=0; i<nrow; i++ )
    {
        Abs[i] = 0.;

        for ( j=0; j<ncol; j++ )
        {
            Abs[i] += As[i][j]*bs[j];
        }
    }
}

void MVmult_dinamic ( double** Ad, double* bd, double* Abd, int nrow, int ncol )
{
    int i, j;

    for ( i=0; i<nrow; i++ )
    {
        Abd[i] = 0.;

        for ( j=0; j<ncol; j++ )
        {
            Abd[i] += Ad[i][j]*bd[j];
        }
    }
}