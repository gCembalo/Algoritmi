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
    fdata << setprecision ( 5 );

    int t, i;
    int ni = 1000, nr = 1000, nd;                                             // rispettivamente atomi iniziali, rimasti, decaduti
    
    double x, lambda = 0.01;
  
    srand48 ( time ( NULL ) );                                                // meglio lanciare una volta sola fuori 

    fdata.open ( "decay.dat" );

    for ( t=0; t<500; t++ )                                 
    {
        fdata << t << " " << ni << endl;

        nd = 0;

        for ( i=0; i<ni; i++ )
        {
            x = drand48();
            
            if ( x<=lambda ) 
            {
                nd += 1;
            }
        } 

        nr = ni-nd;
        ni = nr; 

        if ( nr==0 ) break;
    }

    fdata.close ();

    cout << endl;
    cout << " (vedi file associato) " << endl;
    cout << endl;

    return 0;
}