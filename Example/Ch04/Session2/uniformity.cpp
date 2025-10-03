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
    
    int i;
  
    srand48 ( time ( NULL ) ); 

    fdata.open ( "uniformity.dat" );

    for ( i=0; i<1000; i++ )  
    {
        fdata << i << " " << drand48() << endl;
    }
    
    fdata << endl << endl;

    int n;
    
    double w, m1, m2, err1, err2;
    
    for ( n=4; n<1e6; n=2*n )
    {
        m1 = 0.;
        m2 = 0.;
        w = 0.;
        
        for ( i=0; i<n; i++ )
        {
            w = drand48();
            m1 += w;
            m2 += w*w;
        }
    
        m1 = m1/(double)n;
        m2 = m2/(double)n;

        err1 = fabs(m1-0.5);                                                             // differenza tra sommatoria e funzione approssimata con k=1
        err2 = fabs(m2-1./3.);                                                           // differenza tra sommatoria e funzione approssimata con k=2

        fdata << n << " " << err1 << " " << err2 << endl;
    }

    fdata.close ();

    cout << endl;
    cout << " (vedi file associato) " << endl;
    cout << endl;
  
    return 0;
}