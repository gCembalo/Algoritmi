#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

int main ()
{
    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 6 );
    
    int i, n_in, n = 4;
    
    double x, y, tol = 1e-4;
    double I, err, s;
  
    srand48 ( time ( NULL ) );

    fdata.open ( "pi.dat" );
    
    do
    {
        n_in = 0;
        I = 0.;
        err = 0.;
        s = 0.;
        
        for ( i=1; i<=n; i++ )
        {
            x = (drand48()-0.5)*2;
            y = (drand48()-0.5)*2;
            
            if ( x*x+y*y<=1. )
            {
                n_in += 1;
            }
        }

        I = 4*((double)n_in/(double)n);
        err = fabs((I/M_PI)-1.0);
        s = sqrt(1/(double)(n-1)*((double)n_in/(double)n-((double)n_in/(double)n)*((double)n_in/(double)n)));       
    
        fdata << n << " " << s << endl;

        n = 2*n;
    }
    while ( err>tol );
    
    fdata.close ();

    cout << endl;
    cout << " (vedi file associato) " << endl;
    cout << endl;

    return 0;
}