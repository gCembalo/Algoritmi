#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double func(double &);
void bisection(double (*F)(double &), double, double, double, double&);
void false_position(double (*F)(double &), double, double, double, double&);
void secant_method(double (*F)(double &), double, double, double, double &);
double derfunc(double &);
void newton_method(double (*F)(double &), double (*derF)(double &), double, double, double &);

int main(){

    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    double tol = 1e-7; // tolleranza
    double a = -1. , b = 1.; // gli estremi dell'intervallo in cui sappiamo esserci lo zero

    cout << "\n+-----------------------------------------------------------------------------+\n              e^(-x) - x\n\nLo zero è a x = " << 0.56714329 << endl;

    double x0; // il nostro zero della funzione
    cout << "\n+-----------------------------------------------------------------------------+\n              Metodo di Bisezione\n" << endl;

    bisection(func, a, b, tol, x0);

    cout << "\nLo zero è: " << x0 << " con una tolleranza: " << setprecision(0) << tol << endl;

    //cout << "\n+-----------------------------------------------------------------------------+\n" << endl;

    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    double x1; // il nostro zero della funzione
    cout << "\n+-----------------------------------------------------------------------------+\n              Metodo di False position\n" << endl;

    false_position(func, a, b, tol, x1);

    cout << "\nLo zero è: " << x1 << " con una tolleranza: " << setprecision(0) << tol << endl;

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    double x2; // il nostro zero della funzione
    cout << "\n+-----------------------------------------------------------------------------+\n              Metodo secant\n" << endl;

    secant_method(func, a, b, tol, x2);

    cout << "\nLo zero è: " << x2 << " con una tolleranza: " << setprecision(0) << tol << endl;

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    double x3; // il nostro zero della funzione
    cout << "\n+-----------------------------------------------------------------------------+\n              Metodo Newton\n" << endl;

    newton_method(func, derfunc, a, tol, x3);

    cout << "\nLo zero è: " << x3 << " con una tolleranza: " << setprecision(0) << tol << endl;

    cout << "\n+-----------------------------------------------------------------------------+\n" << endl;

    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    cout << "\nVediamo che il risultato del metodo di Bisezione si discosta di: " << fabs(x0 - 0.56714329) << "\n\nPer il metodo di false position: " << fabs(x1 - 0.56714329) << "\n\nPer il metodo secant: " << fabs(x2 - 0.56714329) << "\n\nPer il metodo Newton: " << fabs(x3 - 0.56714329) << endl;
    cout << "\n+-----------------------------------------------------------------------------+\n" << endl;

    return 0;
}




double func(double &x){
    return exp(-x) - x;
}

double derfunc(double &x){
    return -exp(-x) - 1;
}

void bisection(double (*F)(double &), double a, double b, double tol, double &zero){
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

    // metto nel ciclo la condizione sia sulla tolleranza che sul numero di cicli
    while( fabs( x - xk ) > tol and n < 100 ){
        
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

void secant_method(double (*F)(double &), double a, double b, double tol, double &zero){

    // definisco le variabili che mi servono per tenere traccia delle varie iterazioni di x
    double xk1 = a, xk = b, xk2 = xk + 1; // dove uso xk come x_k, xk1 come x_{k-1} e xk2 come x_{k+1} ; inizializzo gli zeri sugli estremi dell'intervallo in cui ricaviamo la retta
    int n = 0; // variabile per contare
    double xp = 0; // una variabile di controllo per vedere di quanto miglioriamo la guess

    // metto nel ciclo la condizione sia sulla tolleranza che sul numero di cicli
    while( fabs( xk2 - xp ) > tol and n < 100 ){

        n++;
        xp = xk2;

        // calcolo lo zero x_{k+1}
        xk2 = xk - F(xk)*( xk - xk1 )/( F(xk) - F(xk1) );

        xk1 = xk;
        xk = xk2;

        // creo l'output voluto
        cout << "n = " << n << ";   [a,b] = [" << xk1 << ", " << xk << "];    x0 = " << xk2 << ";   Deltax = " << fabs(xk2-xk1) << ";   f(x0) = " << F(xk2) << endl;

    }

    zero = xk2;
}

void newton_method(double (*F)(double &), double (*derF)(double &), double a, double tol, double &zero){

    // definisco le variabili che mi servono per tenere traccia delle varie iterazioni di x
    double xk = a, xk1 = xk + 1; // dove uso xk come x_k, xk1 come x_{k+1}
    int n = 0; // variabile per contare

    double xp = 0; // una variabile di controllo per vedere di quanto miglioriamo la guess
    
    // metto nel ciclo la condizione sia sulla tolleranza che sul numero di cicli
    while( fabs( xp - xk ) > tol and n < 100 ){

        n++;
        xp = xk;

        // calcolo lo zero x_{k+1}
        xk1 = xk - F(xk)/derF(xk);

        // creo l'output voluto
        cout << "n = " << n << ";    xc = " << xk1 << ";   Deltax = " << fabs(xk1-xk) << ";   f(x0) = " << F(xk1) << endl;

        xk = xk1;

    }

    zero = xk;
}