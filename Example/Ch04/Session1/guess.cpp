#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

int main ()
{
    int g, n, i = 0;
    int inf = 1, sup = 100;

    srand ( time ( NULL ) );                                                          // inizializza una sequenza diversa a ogni avvio

    n = rand()%100+1;                                                                 // genera un numero casuale in [1,100]
    
    cout << endl;

    while ( g!=n )
    {
        cout << "Indovina n intero in [" << inf << "," << sup << "]" << endl;
        
        cin >> g;
        
        if ( g<n && g>=inf )
        { 
            inf = g+1;
        }
        else if ( g>n && g<=sup )
        {
            sup = g-1;
        }
        else if ( g>sup || g<inf )
        {
            cout << "Ritenta attentamente!" << endl;
        }
     
        i++; 
    }

    cout << "Hai indovinato in " << i << " tentativi!" << endl;
    cout << endl;

    return 0;
}