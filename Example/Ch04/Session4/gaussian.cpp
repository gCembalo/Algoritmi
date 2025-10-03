#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

double FunzGauss ( double );

int main ()
{
    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 6 );
    
    int k, n = 1e5;
    
    double  x, y, c = 1.;

    srand48 ( time ( NULL ) );                                              

    fdata.open ( "gaussian.dat" );

    for ( k=0; k<n; k++ )
    {
        x = 10*drand48()-5;
        y = c*drand48();

        if ( y<=FunzGauss(x) )
        {
            fdata << x << " " << y << endl;
        }
    }

    fdata.close ();

    cout << endl;
    cout << " (vedi file associato) " << endl;
    cout << endl;

    return 0;
}

double FunzGauss ( double x )
{
    return 2/sqrt(2*M_PI)*exp(-0.5*(2*x)*(2*x));
}