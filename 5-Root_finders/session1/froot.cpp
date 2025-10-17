#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double func(double &);
void Bisection(double (*F)(double &), double, double, double, double&);
void false_position(double (*F)(double &), double, double, double, double&);

int main(){

    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    double tol = 1e-7; // tolleranza
    double a = -1. , b = 1.; // gli estremi dell'intervallo in cui sappiamo esserci lo zero

    double x0; // il nostro zero della funzione
    cout << "\n+-----------------------------------------------------------------------------+\n              Metodo di Bisezione\n" << endl;

    Bisection(func, a, b, tol, x0);

    cout << "\nLo zero è: " << x0 << " con una tolleranza: " << tol << endl;

    cout << "\n+-----------------------------------------------------------------------------+\n" << endl;

    double x1; // il nostro zero della funzione
    cout << "\n+-----------------------------------------------------------------------------+\n              Metodo di False position\n" << endl;

    false_position(func, a, b, tol, x1);

    cout << "\nLo zero è: " << x1 << " con una tolleranza: " << tol << endl;

    cout << "\n+-----------------------------------------------------------------------------+\n" << endl;

    return 0;
}




double func(double &x){
    return exp(-x) - x;
}

void Bisection(double (*F)(double &), double a, double b, double tol, double &zero){
    double x; // la mia guess dello zero che aggiorno ad ogni iterazione
    int n=0; // la variabile che mi permette di contare le iterazioni

    while( fabs(a-b) > tol ){
        
        n++;
        x = ( a+b ) / 2;

        // controllo dove si trova x rispetto gli estremi a e b
        if( F(a)*F(x) < 0 ){
            b = x;
        }
        else if ( F(a)*F(x) > 0 ){
            a = x;
        }

        // creo l'output voluto
        cout << "n = " << n << ";   [a,b] = [" << a << ", " << b << "];    xm = " << x << ";   Deltax = " << fabs(a-b) << ";   f(xm) = " << F(x) << endl;
    }

    zero = x;

}

void false_position(double (*F)(double &), double a, double b, double tol, double &zero){
    double x = 3; // la guess di zero della funzione
    double xk = 0; // variabile che mi serve per la tolleranza
    double m, q; // parametri della retta
    int n = 0; // variabile per contare

    while( fabs( x - xk ) > tol ){
        
        n++;
        xk = x;
        
        // trovo la retta (con conti sul quaderno)
        m = ( F(a) - F(b) ) / (a - b);
        q = ( F(b)*a - F(a)*b ) / (a - b);
        // trovo lo zero della retta e lo chiamo x
        x = - q / m;

        // controllo dove si trova x rispetto gli estremi a e b
        if( F(a)*F(x) < 0 ){
            b = x;
        }
        else if ( F(a)*F(x) > 0 ){
            a = x;
        }

        // creo l'output voluto
        cout << "n = " << n << ";   [a,b] = [" << a << ", " << b << "];    xm = " << x << ";   Deltax = " << fabs(a-b) << ";   f(xm) = " << F(x) << endl;
    }

    zero = x;
}