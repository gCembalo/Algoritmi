#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

ofstream fdata;

void Tridiagonal ( double*, double*, double*, double*, double*, int );

double Func ( double );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << fixed << setprecision ( 4 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 4 );
    
    static int n = 32;

    int imax, i, j, k;

    double *a, *b, *c, *r, *y;
    double x, h;
    double xi = 0., xf = 1.;
    double alpha, beta;

    a = new double [n];
    b = new double [n];
    c = new double [n];
    r = new double [n];
    y = new double [n];

    h = fabs(xf-xi)/(double)(n-1);

    y[0] = alpha;
    y[n-1] = beta;

    for ( i=0; i<n; i++ )
    {
        x = xi+i*h;
        a[i] = 1;
        b[i] = -2;
        c[i] = 1;
        r[i] = h*h*Func(x);
    }

    r[1] -= y[0];
    r[n-2] -= y[n-1];
    
    Tridiagonal ( a+1, b+1, c+1, r+1, y+1, n-2 );

    fdata.open ( "bvp.dat" );

    for ( i=0; i<n; i++ )
    {
        x = xi+i*h;
        fdata << x << " " << y[i] << endl;
    }

    fdata.close ();

    cout << endl;
    cout << " (vedi file associato) " << endl;
    cout << endl;

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] r;
    delete[] y;

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

double Func ( double x )
{
    return 1.;
}