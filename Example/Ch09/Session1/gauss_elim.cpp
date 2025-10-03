#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

void Elimination ( double**, double*, double*, int );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << fixed << setprecision ( 4 );

    int i, n;

    double **A, *b, *x;

    n = 4;

    A = new double *[n];
    A[0] = new double [n*n];
    b = new double [n];
    x = new double [n];

    for ( i=1; i<n; i++ )
    {
        A[i] = A[i-1]+n;
    }

    A[0][0] = 1.;   A[0][1] = 2.;   A[0][2] = 1.;   A[0][3] = -1.;
    A[1][0] = 3.;   A[1][1] = 2.;   A[1][2] = 4.;   A[1][3] = 4.;
    A[2][0] = 4.;   A[2][1] = 4.;   A[2][2] = 3.;   A[2][3] = 4.;
    A[3][0] = 2.;   A[3][1] = 0.;   A[3][2] = 1.;   A[3][3] = 5.; 

    b[0] = 5.;   
    b[1] = 16.;   
    b[2] = 22.;   
    b[3] = 15.;

    Elimination ( A, b, x, n );

    cout << endl;
    cout << "Sistema n = " << n << endl;
    cout << "-------------" << endl;

    for ( i=0; i<n; i++ )
    {
        cout << setw(10) << right << x[i] << endl;
    }

    delete[] A[0];
    delete[] A;
    delete[] b;
    delete[] x;

    n = 3;

    A = new double *[n];
    A[0] = new double [n*n];
    b = new double [n];
    x = new double [n];

    for ( i=1; i<n; i++ )
    {
        A[i] = A[i-1]+n;
    }

    A[0][0] = 2.;   A[0][1] = -1.;  A[0][2] = 0.;
    A[1][0] = -1.;  A[1][1] = 2.;   A[1][2] = -1.;
    A[2][0] = 0.;   A[2][1] = -1.;  A[2][2] = 2.;

    b[0] = 1.;   
    b[1] = 2.;   
    b[2] = 1.;

    Elimination ( A, b, x, n );

    cout << endl;
    cout << "Sistema n = " << n << endl;
    cout << "-------------" << endl;

    for ( i=0; i<n; i++ )
    {
        cout << setw(10) << right << x[i] << endl;
    }

    cout << endl;

    delete[] A[0];
    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}

void Elimination ( double** A, double* b, double* x, int n )
{
    int i, j, k;
    
    double g, tmp;

    for ( k=0; k<n-1; k++ )                                                          // Loop over the Gk's
    {
        for ( i=k+1; i<n; i++ )                                                      // Loop over rows
        { 
            g = A[i][k]/A[k][k];
            
            for ( j=k+1; j<n; j++ ) 
            {
                A[i][j] -= g*A[k][j];
            }
            
            A[i][k] = 0.;
            b[i] -= g*b[k];
        }
    }

    for ( i=n-1; i>=0; i-- )                                                         // system solution
    {
        tmp = b[i];
        
        for ( j=n-1; j>i; j-- ) 
        {
            tmp -= x[j]*A[i][j];
        }
        
        x[i] = tmp/A[i][i];
    }
}