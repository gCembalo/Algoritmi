#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

void Tridiagonal ( double*, double*, double*, double*, double*, int );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << fixed << setprecision ( 4 );

    int i, n = 5;

    double *a, *b, *c, *r, *x;

    a = new double [n];
    b = new double [n];
    c = new double [n];
    r = new double [n];
    x = new double [n];

    a[0] = 0.;   a[1] = 1.;   a[2] = 1.;   a[3] = 1.;   a[4] = 1.;                                   // a[0] qualsiasi (non verrà utilizzato)
    b[0] = 2.;   b[1] = 2.;   b[2] = 2.;   b[3] = 2.;   b[4] = 2.;
    c[0] = 1.;   c[1] = 1.;   c[2] = 1.;   c[3] = 1.;   c[4] = 0.;                                   // c[4] qualsiasi (non verrà utilizzato)
    r[0] = 1.;   r[1] = 0.;   r[2] = 3.;   r[3] = 1.;   r[4] = 0.;                               

    Tridiagonal ( a, b, c, r, x, n );

    cout << endl;
    cout << "Sistema n = " << n << endl;
    cout << "-------------" << endl;

    for ( i=0; i<n; i++ )
    {
        cout << setw(10) << right << x[i] << endl;
    }

    cout << endl;

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] r;
    delete[] x;

    return 0;
}

void Tridiagonal ( double* a, double* b, double* c, double* r, double* x, int n )
{
    int i;

    double *h, *p;
    
    h = new double [n];
    p = new double [n];

    h[0] = c[0]/b[0];
    p[0] = r[0]/b[0]; 

    for ( i=1; i<n; i++ )
    {
        h[i] = c[i]/(b[i]-a[i]*h[i-1]);
        p[i] = (r[i]-a[i]*p[i-1])/(b[i]-a[i]*h[i-1]);
    }

    x[n-1] = p[n-1];

    for ( i=n-2; i>=0; i-- )
    {
        x[i] = p[i]-h[i]*x[i+1];
    }

    delete[] h;
    delete[] p;
}