#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

int main ()
{
    double x1, x2, x3, x4, x5, x6;
    double a, b, c, delta;
    
    cout << endl;                                                           // inserisco le radici che voglio trovare
    cout << "Inserisci x1 e x2" << endl;

    cin >> x1 >> x2;

    a = 1.;                                                                 // meglio specificare il decimale
    b = -(x1+x2);
    c = x1*x2;
    delta = sqrt(b*b-4*a*c);                                                // meglio calcolare il delta solo una volta

    cout << endl;                                                           // test della formula standard
    cout << "Risultati con solo la formula standard" << endl;       
    
    x3 = 0.5*(-b+delta)/a;                                                  // meglio moltiplicare per 0.5 che dividere per 2
    x4 = 0.5*(-b-delta)/a;

    cout << "x+ = " << x3 << endl;
    cout << "x- = " << x4 << endl;

    cout << endl;                                                           // test della formula implementata
    cout << "Risultati con anche la formula razionalizzata" << endl;

    if ( b<0 )                                                              // con b=0 è totalmente indifferente
    {
        x5 = 2*c/(-b+delta);
        x6 = 0.5*(-b+delta)/a;

        cout << "x+ = " << x6 << endl;
        cout << "x- = " << x5 << endl;
    }
    else if ( b>0 )
    {
        x5 = 0.5*(-b-delta)/a;
        x6 = 2*c/(-b-delta);

        cout << "x+ = " << x6 << endl;
        cout << "x- = " << x5 << endl; 
    }
    else                                                                    // solo per avere tutte le combinazioni possibili, teoricamente non necessario
    {
        x5 = 2*c/delta;
        x6 = -2*c/delta;

        cout << "x+ = " << x6 << endl;
        cout << "x- = " << x5 << endl;
    }

    cout << endl;

    return 0;
}