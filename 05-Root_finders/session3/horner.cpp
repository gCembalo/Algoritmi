// Consider the function f(x) = x^3 - 3x^2 + x + 5
// This function has a root in x = -1. Repeat the search over the interval [-5,0]. The following table gives the number of iterations obtained with the different methods using a tolerance 10-8:
// Bisection: 30
// False position: 80
// Secant: 12
// Newton: 6
//
// Can you explain why False Position performs so badly ? What happens when the initial interval is reduced to [-2,0] ?
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double pol(double);
double derpol(double);
int bisection(double (*F)(double), double, double, double, double&);
int false_position(double (*F)(double), double, double, double, double&);
int secant_method(double (*F)(double), double, double, double, double &);
int newton_method(double (*F)(double), double (*derF)(double), double, double, double &);

int main(){

    // definisco gli estremi dell'intervallo
    double a = -5, b = 0;
    // defiisco una tolleranza
    double tol = 1.e-8;
    //double facc = 1e-6;
    // definisco gli zeri
    double x0, x1, x2, x3;

    cout << "\n+-----------------------------------------------------------------------------+\nLo zero è a x = " << -1 << endl;

    cout << "prendendo l'intervallo [ " << a << " , " << b << " ]\n" << endl;

    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    cout << "\nCon # iterazioni:\n" << endl;
     // stampo gli algoritmi
    bisection(pol, a, b, tol, x0);
    false_position(pol, a, b, tol, x1);
    secant_method(pol, a, b, tol, x2);
    newton_method(pol, derpol, a, tol, x3);

    cout << "\nI diversi metodi restituiscono:\n" << endl;

    cout << "(Bisezione): " << x0 << "\n(False position): " << x1 << "\n(Secant): " << x2 << "\n(Newton): " << x3 << endl;

    cout << "\n+-----------------------------------------------------------------------------+\n" << endl;

    return 0;
}




double pol(double x){
    
    // definisco l'array e il polinomio
    double a[4] = {5, 1, -3, 1};
    // inizializzo p come a_n
    double p = a[3];

    // valuto il polinomio con il metodo di horner
    for( int j = 3-1 ; j >= 0 ; j-- ){

        // moltiplico per x il termine in p
        p = a[j] + p*x;
    }

    return p;
}

double derpol(double x){

    // definisco l'array e il polinomio
    double a[3] = {1, -6, 3};
    //double a[4] = {5, 1, -3, 1};
    // inizializzo p come a_n
    double p = a[2];
    double dp = 0;

    // valuto il polinomio con il metodo di horner
    for( int j = 2-1 ; j >= 0 ; j-- ){

        // moltiplico per x il termine in p
        //dp = dp*x + p; // devi mettere p[2] e int j = 2-1 e usare a[4]
        p = a[j] + p*x;
    }

    return p;
}

int bisection(double (*F)(double), double a, double b, double tol, double &zero){

    double x; // la mia guess dello zero che aggiorno ad ogni iterazione
    int n = 0; // la variabile che mi permette di contare le iterazioni

    // metto i controlli di non avere già uno zero
    if( F(a) == 0.0 ){
        zero = a;
        // creo l'output delle iterazioni
        cout << "(Bisection) # = " << n << "    (l'estremo " << a << " è già lo zero)"<< endl;
        return 0;
    }
    else if( F(b) == 0.0 ){
        zero = b;
        // creo l'output delle iterazioni
        cout << "(Bisection) # = " << n << "(l'estremo " << b << " è già lo zero)"<< endl;
        return 0;
    }
    else{
        while( fabs(a-b) > tol and n<100 ){
        
            n++;
            x = ( a+b ) / 2;

            if( F(x) == 0 ){
                zero = x;
                // creo l'output delle iterazioni
                cout << "(Bisection) # = " << n << endl;
                return 0;
            }
            // controllo dove si trova x rispetto gli estremi a e b
            else if( F(a)*F(x) < 0 ){
                b = x;
            }
            else if ( F(a)*F(x) > 0 ){
                a = x;
            }

            // creo l'output voluto
            //cout << "n = " << n << ";   [a,b] = [" << a << ", " << b << "];    xm = " << x << ";   Deltax = " << fabs(a-b) << ";   f(xm) = " << F(x) << endl;

        }

        // creo l'output delle iterazioni
        cout << "(Bisection) # = " << n << endl;

        zero = x;
        return 0;

    }

}

int false_position(double (*F)(double), double a, double b, double tol, double &zero){
    double x = 3; // la guess di zero della funzione
    double xk = 0; // variabile che mi serve per la tolleranza
    double m, q; // parametri della retta
    int n = 0; // variabile per contare

    // metto i controlli di non avere già uno zero
    if( F(a) == 0.0 ){

        zero = a;
        // creo l'output delle iterazioni
        cout << "(False position) # = " << n << "    (l'estremo " << a << " è già lo zero)"<< endl;
        return 0;

    }
    else if( F(b) == 0.0 ){

        zero = b;
        // creo l'output delle iterazioni
        cout << "(False position) # = " << n << "(l'estremo " << b << " è già lo zero)"<< endl;
        return 0;

    }
    else{
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
            //cout << "n = " << n << ";   [a,b] = [" << a << ", " << b << "];    xm = " << x << ";   Deltax = " << fabs(a-b) << ";   f(xm) = " << F(x) << endl;

        }

        // creo l'output delle iterazioni
        cout << "(False position) # = " << n << endl;

        zero = x;
        return 0;

    }
}

int secant_method(double (*F)(double), double a, double b, double tol, double &zero){

    // definisco le variabili che mi servono per tenere traccia delle varie iterazioni di x
    double xk1 = a, xk = b, xk2 = xk + 1; // dove uso xk come x_k, xk1 come x_{k-1} e xk2 come x_{k+1} ; inizializzo gli zeri sugli estremi dell'intervallo in cui ricaviamo la retta
    int n = 0; // variabile per contare
    double xp = 0; // una variabile di controllo per vedere di quanto miglioriamo la guess

    // metto i controlli di non avere già uno zero
    if( F(a) == 0.0 ){

        zero = a;
        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n << "    (l'estremo " << a << " è già lo zero)" << endl;
        return 0;

    }
    else if( F(b) == 0.0 ){

        zero = b;
        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n << "(l'estremo " << b << " è già lo zero)"<< endl;
        return 0;

    }
    else{
        
        // metto nel ciclo la condizione sia sulla tolleranza che sul numero di cicli
        while( fabs( xk2 - xp ) > tol and n < 100 ){

            n++;
            xp = xk2;

            // calcolo lo zero x_{k+1}
            xk2 = xk - F(xk)*( xk - xk1 )/( F(xk) - F(xk1) );

            xk1 = xk;
            xk = xk2;

            // creo l'output voluto
            //cout << "n = " << n << ";   [a,b] = [" << xk1 << ", " << xk << "];    x0 = " << xk2 << ";   Deltax = " << fabs(xk2-xk1) << ";   f(x0) = " << F(xk2) << endl;

        }

        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n << endl;

        zero = xk2;
        return 0;

    }

}

int newton_method(double (*F)(double), double (*derF)(double), double a, double tol, double &zero){

    // definisco le variabili che mi servono per tenere traccia delle varie iterazioni di x
    double xk = a, xk1 = xk + 1; // dove uso xk come x_k, xk1 come x_{k+1}
    int n = 0; // variabile per contare

    double xp = 0; // una variabile di controllo per vedere di quanto miglioriamo la guess

    // metto i controlli di non avere già uno zero
    if( F(a) == 0.0 ){

        zero = a;
        // creo l'output delle iterazioni
        cout << "(Newton) # = " << n << "    (l'estremo " << a << " è già lo zero)"<< endl;
        return 0;

    }
    else{
        
        // metto nel ciclo la condizione sia sulla tolleranza che sul numero di cicli
        while( fabs(F(xk)) > tol and n < 100 ){

            // controllo di non avere una derivata nulla
            if (fabs(derF(xk)) < 1e-12) {

            cout << "Errore: derivata troppo piccola.\n" << endl;
            return 0;
            }

            n++;
            xp = xk;

            // calcolo lo zero x_{k+1}
            xk1 = xk - F(xk)/derF(xk);

            // creo l'output voluto
            //cout << "n = " << n << ";    xc = " << xk1 << ";   Deltax = " << fabs(xk1-xk) << ";   f(x0) = " << F(xk1) << endl;

            xk = xk1;

        }

        // creo l'output delle iterazioni
        cout << "(Newton) # = " << n << endl;

        zero = xk;
        return 0;
    }
    
}