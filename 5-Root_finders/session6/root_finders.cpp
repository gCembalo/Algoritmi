#include "root_finders.h"

using namespace std;

int bisection(double (*F)(double), double a, double b, double tol, double &zero, int &l){

    double x; // la mia guess dello zero che aggiorno ad ogni iterazione
    int n = 0; // la variabile che mi permette di contare le iterazioni

    // metto i controlli di non avere già uno zero
    if( F(a) == 0.0 ){
        zero = a;
        // creo l'output delle iterazioni
        cout << "(Bisection) # = " << n << "    (l'estremo " << a << " è già lo zero)"<< endl;
        l = n;
        return 0;
    }
    else if( F(b) == 0.0 ){
        zero = b;
        // creo l'output delle iterazioni
        cout << "(Bisection) # = " << n << "(l'estremo " << b << " è già lo zero)"<< endl;
        l = n;
        return 0;
    }
    else{
        while( fabs(a-b) > tol ){
        
            n++;

            if( n == 100 ){
                cout << "(Bisection) Troppe iterazioni." << endl;
                l = n;
                return 0;
            }

            x = ( a+b ) / 2;

            if( F(x) == 0 ){
                zero = x;
                // creo l'output delle iterazioni
                //cout << "(Bisection) # = " << n << endl;
                l = n;
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
        //cout << "(Bisection) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = x;
        return 0;

    }

}

int false_position(double (*F)(double), double a, double b, double tol, double &zero, int &l){
    double x = 3; // la guess di zero della funzione
    double xk = 0; // variabile che mi serve per la tolleranza
    double m, q; // parametri della retta
    int n = 0; // variabile per contare

    // metto i controlli di non avere già uno zero
    if( F(a) == 0.0 ){

        zero = a;
        // creo l'output delle iterazioni
        cout << "(False position) # = " << n << "    (l'estremo " << a << " è già lo zero)"<< endl;
        l = n;
        return 0;

    }
    else if( F(b) == 0.0 ){

        zero = b;
        // creo l'output delle iterazioni
        cout << "(False position) # = " << n << "(l'estremo " << b << " è già lo zero)"<< endl;
        l = n;
        return 0;

    }
    else{
        // metto nel ciclo la condizione sia sulla tolleranza che sul numero di cicli
        while( fabs( x - xk ) > tol ){
        
            n++;

            if( n == 100 ){
                cout << "(False position) Troppe iterazioni." << endl;
                l = n;
                return 0;
            }

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
        //cout << "(False position) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = x;
        return 0;

    }
}

int secant_method(double (*F)(double), double a, double b, double tol, double &zero, int &l){

    // definisco le variabili che mi servono per tenere traccia delle varie iterazioni di x
    double xk1 = a, xk = b, xk2 = xk + 1; // dove uso xk come x_k, xk1 come x_{k-1} e xk2 come x_{k+1} ; inizializzo gli zeri sugli estremi dell'intervallo in cui ricaviamo la retta
    int n = 0; // variabile per contare
    double xp = 0; // una variabile di controllo per vedere di quanto miglioriamo la guess

    // metto i controlli di non avere già uno zero
    if( F(a) == 0.0 ){

        zero = a;
        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n << "    (l'estremo " << a << " è già lo zero)" << endl;
        l = n; // sono le iterazioni
        return 0;

    }
    else if( F(b) == 0.0 ){

        zero = b;
        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n << "(l'estremo " << b << " è già lo zero)"<< endl;
        l = n; // sono le iterazioni
        return 0;

    }
    else{
        
        // metto nel ciclo la condizione sia sulla tolleranza che sul numero di cicli
        while( fabs( xk2 - xp ) > tol ){

            n++;

            if( n == 100 ){
                cout << "(Secant) Troppe iterazioni." << endl;
                l = n; // sono le iterazioni
                return 0;
            }

            xp = xk2;

            // calcolo lo zero x_{k+1}
            xk2 = xk - F(xk)*( xk - xk1 )/( F(xk) - F(xk1) );

            xk1 = xk;
            xk = xk2;

            // creo l'output voluto
            //cout << "n = " << n << ";   [a,b] = [" << xk1 << ", " << xk << "];    x0 = " << xk2 << ";   Deltax = " << fabs(xk2-xk1) << ";   f(x0) = " << F(xk2) << endl;

        }

        // creo l'output delle iterazioni
        //cout << "(Secant) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = xk2;
        return 0;

    }

}

int newton_method(double (*F)(double), double (*derF)(double), double a, double tol, double &zero, int &l){

    // definisco le variabili che mi servono per tenere traccia delle varie iterazioni di x
    double xk = a, xk1 = xk + 1; // dove uso xk come x_k, xk1 come x_{k+1}
    int n = 0; // variabile per contare

    double xp = 0; // una variabile di controllo per vedere di quanto miglioriamo la guess

    // metto i controlli di non avere già uno zero
    if( F(a) == 0.0 ){

        zero = a;
        // creo l'output delle iterazioni
        cout << "(Newton) # = " << n << "    (l'estremo " << a << " è già lo zero)"<< endl;
        l = n; // sono le iterazioni
        return 0;

    }
    else{
        
        // metto nel ciclo la condizione sia sulla tolleranza che sul numero di cicli
        while( fabs(F(xk)) > tol ){

            // controllo di non avere una derivata nulla
            if (fabs(derF(xk)) < 1.e-12) {

            cout << "(Newton) Errore: derivata troppo piccola.\n" << endl;
            l = n; // sono le iterazioni
            return 0;
            
            }

            n++;

            if( n == 100 ){
                cout << "(Newton) Troppe iterazioni." << endl;
                l = n; // sono le iterazioni
                return 0;
            }

            xp = xk;

            // calcolo lo zero x_{k+1}
            xk1 = xk - F(xk)/derF(xk);

            // creo l'output voluto
            //cout << "n = " << n << ";    xc = " << xk1 << ";   Deltax = " << fabs(xk1-xk) << ";   f(x0) = " << F(xk1) << endl;

            xk = xk1;

        }

        // creo l'output delle iterazioni
        //cout << "(Newton) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = xk;
        return 0;
    }

}


int Bracket(double (*F)(double), double a, double b, double n, double *xL, double *xR, int &nroots){

    // definisco le variabili che uso per contare
    int count = 0, i; // count mi dice quanti zeri ho
    double dx = ( b - a )/(double)n; // spacing
    double aL, aR; // estremi di ogni sotto intervallo

    double fL, fR; // valori delle funzioni valutate agli estremi
    fL = F(a);

    // faccio il loop su tutti gli intervalli
    for( i = 0 ; i < n ; i++ ){

        aL = a + i*dx;
        aR = a + (i+1)*dx;

        // metto la condizione in cui abbiamo un cambiamento di segno e quindi con la quale ci ricordiamo il valore
        fR = F(aR);
        if( fL*fR <= 0.0 ){

            // metto gli estremi nell'array
            xL[count] = aL;
            xR[count] = aR;
            count++; // aggiorno il numero di roots

        }
        
        fL = fR;

    }

    nroots = count;
    return 0;

}