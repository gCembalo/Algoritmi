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

int newton_method(double (*F)(double), double (*derF)(double), double a, double b, double xtol, double ytol, double &zero, int &l){

    int n = 0; // variabile per contare

    // controllo che gli estremi non siano già degli zeri
    double fa = F(a), fb = F(b);
    if( fabs(fa) < ytol ){

        zero = a;
        // creo l'output delle iterazioni
        cout << "(Newton) # = " << n << "    (l'estremo " << a << " è già lo zero)"<< endl;

        return 0;

    }

    if( fabs(fb) < ytol ){
        
        zero = b;
        // creo l'output delle iterazioni
        cout << "(Newton) # = " << n << "    (l'estremo " << b << " è già lo zero)"<< endl;

        return 0;

    }

    // definisco delle variabili utili all'algoritmo
    double Deltax = fabs( b-a ); // l'ampiezza dell'intervallo inziale
    double deltax = Deltax * 0.5; // la nuova ampiezza (inizializzata)
    double x = ( a + b )*0.5; // la mia prima guess, che prendo a metà dell'intervallo
    double fx = F(x); // funzione valutata nella guess

    // faccio il ciclo per il metodo di Newton, con i controlli su x, su y e sul numero di iterazioni (nel mezzo del while)
    while( fabs(deltax) > xtol && fabs(fx) > ytol ){

        n++;

        // blocco il ciclo se raggiungo troppe iterazioni
        if( n == 100 ){

            cout << "(Newton) Troppe iterazioni. (" << l << ")" << endl;
            l = n; // sono le iterazioni
            return 0;

        }

        // metto un controllo sugli intervalli, che mi blocca il ciclo se la nuova ampiezza deltax è maggiore di quella precedente, quindi mi blocca se cominciamo a non convergere
        if( fabs(deltax) > fabs(Deltax) ){

            zero = x;
            l = n;
            return 0;

        }

        // controllo anche di non uscire dal mio intervallo di estremi a e b
        if( x < a || x > b ){

            cout << "(Newton) il metodo non converge." << endl;
            return 0;

        }

        // calcolo la derivata e controllo che non sia nulla (troppo piccola)
        double df = derF(x);
        if( fabs( df ) < 1.e-15 ){

            cout << "(Newton) Errore: derivata troppo piccola.\n" << endl;
            l = n; // sono le iterazioni
            return 0;
            
        }

        // superati tutti i controlli implemento il metodo
        Deltax = deltax;
        deltax = fx / df;
        // calcolo la nuova guess
        x = x - deltax;
        // ricalcolo la funzione nella guess
        fx = F(x);

    }

    // creo l'output delle iterazioni
    //cout << "(Newton) # = " << n << endl;

    l = n; // sono le iterazioni
    zero = x;
    return 0;

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