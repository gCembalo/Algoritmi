#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 15 );

    int i = 1;
    
    double S, x0;
    double x, y = 0., err = 1.;

    cout << endl;
    cout << "Inserisci un numero reale positivo" << endl;

    cin >> S;

    cout << "prova a indovinarne la radice quadrata" << endl;

    cin >> x0;
    
    cout << "------------------------------------------------------------------------------" << endl;

    x = x0;

    while ( err!=0 )
    {
        y = 0.5*(x+S/x);
        err = abs(y-x);

        cout << "iterazione # " << i << "\t x = " << y << "\t err = " << err << endl;

        x = y;
        i++;
    }

    cout << "------------------------------------------------------------------------------" << endl;
    cout << "radice quadrata calcolata di " << S << " : \t " << y << endl;
    cout << "vera radice quadrata di " << S << " : \t " << sqrt(S) << endl;
    cout << endl;

    return 0;
}