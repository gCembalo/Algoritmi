#include <iostream>
#include <cmath>
#include <cstdlib>

#define NSIZE 1000                                                   // è una macro cioè una sostituzione che vale in tutto il codice

using namespace std;

void Array ( double*, int, double&, double&, double& );

int main ()
{
    int i, n = NSIZE;

    double a[n], m, v, d;                                            // sempre definire bene la lunghezza massima dell'array con un valore dato
  
    srand48 ( time ( NULL ) );                                       // chiama un seme diverso ogni volta
  
    for ( i=0; i<n; i++ )
    {
        a[i] = drand48();
    }

    Array ( a, n, m, v, d );

    cout << endl;
    cout << "media: \t\t\t" << m << endl;
    cout << "varianza: \t\t" << v << endl;
    cout << "deviazione standard: \t" << d << endl;
    cout << endl;

    return 0;  
}

void Array ( double* a, int n, double& m, double& v, double& d )
{
    int i;
    
    double s = 0., e = 0.;                                           //     s = somma delle n componenti       e = somma degli errori delle n componenti

    for( i=0; i<n; i++ )
    {
        s = s+a[i];
    }

    m = s/n;

    for( i=0; i<n; i++ )
    {
        e = e+(a[i]-m)*(a[i]-m);
    }

    v = e/n;
    d = sqrt(v);
}