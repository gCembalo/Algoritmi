// Name: Gabriele Cembalo
// Date: 28/10/2025
// Code output:
// ************************************************************
// Bisection, results:
// s = 1; Root: 7.123375e+00; ntry = 32
// s = 2; Root: 2.815275e+01; ntry = 32
// s = 3; Root: 6.143373e+01; ntry = 32
// Secant, results:
// s = 1; Root: 7.123375e+00; ntry = 8
// s = 2; Root: 2.815275e+01; ntry = 6
// s = 3; Root: 6.143373e+01; ntry = 5
// ************************************************************

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// definisco come variabili globali i potenziali
double V1 = 250.0;
double V2 = 80.0;
// definisco come variabile globale il quanto di energia
int s;

double Energy(double);
int bisection(double (*F)(double), double, double, double, double&, int &);
int secant_method(double (*F)(double), double, double, double, double &, int &);

int main(){

    cout << setprecision(6);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    // definisco la tolleranza
    double tol = 1.e-8;

    // definisco lo zero e le itezazioni dei metodi
    double x1, x2;
    int l1, l2;

    double a,b;

    // stampo a video gli zeri trovati con i due metodi
    cout << "\n+-----------------------------------------+" << endl;
    cout << "Bisection, results:" << endl;
    for( s = 1 ; s <= 3 ; s++ ){

        // intervalli che definisco per ridurre leggermente 0<E<80; da una prima compilazione vedo che il secondo 0 è a 28, per cui il primo sarà sicuramente minore.
        double a = ( (s - 1) * 20.0 + 1.0 );
        double b = s * 25.0;


        // richiamo la funzione
        bisection(Energy, a, b, tol, x1, l1);

        // stampo l'output
        cout << "s = " << s << "; Root: " << x1 << "; ntry = " << l1 << endl;

    }

    cout << "Secant, results:" << endl;
    for( s = 1 ; s <= 3 ; s++ ){

        // intervalli che definisco per ridurre leggermente 0<E<80; da una prima compilazione vedo che il secondo 0 è a 28, per cui il primo sarà sicuramente minore.
        double a = ( (s - 1) * 20.0 + 1.0 );
        double b = s * 25.0;


        // richiamo la funzione
        secant_method(Energy, a, b, tol, x2, l2);
        
        // stampo l'output
        cout << "s = " << s << "; Root: " << x2 << "; ntry = " << l2 << endl;

    }
    cout << "+-----------------------------------------+" << endl;

    return 0;
}

// definisco la funzione dei livelli energetici
double Energy(double E){

    return sqrt(E) - (double)s*M_PI + asin( sqrt( E/V1 ) ) + asin( sqrt( E/V2 ) );

}

// metodo della bisezione
// gli do in input la funzione, gli estremi a e b, la tolleranza su x, uno zero per riferimento e il numero di iterazioni
int bisection(double (*F)(double), double a, double b, double tol, double &zero, int &l){

    double x; // la mia guess dello zero che aggiorno ad ogni iterazione
    int n = 0; // la variabile che mi permette di contare le iterazioni

    // definisco le variabili della funzione valutata
    double fa = F(a);
    double fb = F(b);
    double fx;

    // metto i controlli di non avere già uno zero
    if( fa == 0.0 ){

        zero = a;

        // creo l'output delle iterazioni
        cout << "(Bisection) # = " << n << "    (l'estremo " << a << " è già lo zero)"<< endl;

        l = n;
        return 0;

    }
    else if( fb == 0.0 ){

        zero = b;

        // creo l'output delle iterazioni
        cout << "(Bisection) # = " << n << "(l'estremo " << b << " è già lo zero)"<< endl;

        l = n;
        return 0;

    }
    else{

        while( fabs(a-b) > tol ){
        
            n++;

            // metto il controllo sul numero di iterazioni
            if( n == 100 ){

                cout << "(Bisection) Troppe iterazioni." << endl;

                l = n;
                return 0;

            }

            // calcolo la prima stima dello zero
            x = ( a+b ) / 2;

            // definisco le variabili delle funzioni valutate
            fa = F(a);
            fb = F(b);
            fx = F(x);

            // controllo se è uno zero
            if( fx == 0 ){

                zero = x;

                // creo l'output delle iterazioni
                //cout << "(Bisection) # = " << n << endl;

                l = n;
                return 0;

            }
            // controllo dove si trova x rispetto gli estremi a e b
            else if( fa*fx < 0 ){

                b = x;

            }
            else if ( fa*fx > 0 ){

                a = x;

            }

            // creo l'output voluto (esercizio froot.cpp)
            //cout << "n = " << n << ";   [a,b] = [" << a << ", " << b << "];    xm = " << x << ";   Deltax = " << fabs(a-b) << ";   f(xm) = " << F(x) << endl;

        }

        // creo l'output delle iterazioni
        //cout << "(Bisection) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = x;
        return 0;

    }

}

// metodo della secante
// gli do in input la funzione, gli estremi a e b, la tolleranza su x, uno zero per riferimento e il numero di iterazioni
int secant_method(double (*F)(double), double a, double b, double tol, double &zero, int &l){

    // definisco le variabili che mi servono per tenere traccia delle varie iterazioni di x
    double xk1 = a, xk = b, xk2 = xk + 1; // dove uso xk come x_k, xk1 come x_{k-1} e xk2 come x_{k+1} ; inizializzo gli zeri sugli estremi dell'intervallo in cui ricaviamo la retta
    int n = 0; // variabile per contare
    double xp = 0; // una variabile di controllo per vedere di quanto miglioriamo la guess

    // definisco le variabili della funzione valutata
    double fa = F(a);
    double fb = F(b);
    double fxk;
    double fxk1;

    // metto i controlli di non avere già uno zero
    if( fa == 0.0 ){

        zero = a;

        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n << "    (l'estremo " << a << " è già lo zero)" << endl;

        l = n; // sono le iterazioni
        return 0;

    }
    else if( fb == 0.0 ){

        zero = b;

        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n << "(l'estremo " << b << " è già lo zero)"<< endl;

        l = n; // sono le iterazioni
        return 0;

    }
    else{
        
        // metto nel ciclo la condizione sia sulla tolleranza che sul numero di cicli (messo dopo)
        while( fabs( xk2 - xp ) > tol ){

            n++;

            if( n == 100 ){

                cout << "(Secant) Troppe iterazioni." << endl;

                l = n; // sono le iterazioni
                return 0;

            }

            xp = xk2;
            fxk = F(xk);
            fxk1 = F(xk1);

            // calcolo lo zero x_{k+1}
            xk2 = xk - fxk*( xk - xk1 )/( fxk - fxk1 );

            xk1 = xk;
            xk = xk2;

            // creo l'output voluto (esercizio froot.cpp)
            //cout << "n = " << n << ";   [a,b] = [" << xk1 << ", " << xk << "];    x0 = " << xk2 << ";   Deltax = " << fabs(xk2-xk1) << ";   f(x0) = " << F(xk2) << endl;

        }

        // creo l'output delle iterazioni
        //cout << "(Secant) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = xk2;
        return 0;

    }

}