#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int main(){

    cout << setprecision(16);
    cout << setiosflags ( ios::scientific );

    double S, x, err, y; // y variabile aggiuntiva per il ciclo
    cout << "Dammi un numero reale: ";
    cin >> S;
    cout << endl << "Dammi la tua guess: ";
    cin >> x;

    //Implemento la funzione
    err = 10; // inizializzo l'errore per poter mettere un controllo sul ciclo

    for(int i=0; i<100; i++){
        if(err > pow(10,-16)){
            y = 0.5*( x + S/x );
            err = fabs(x-y);
            x=y;
        } else{
            break;
        }
        cout << "Iterazione " << i << ": " << y << "; err = " << err << endl;
    }

    //Stampo a video
    cout << endl;
    cout << "La radice quadrata di " << S << " è " << x << endl;
    cout << "con un errore di " << err << endl;
    cout << "La soluzione esatta è: " << sqrt(S) << endl;

    return 0;
}