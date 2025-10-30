// Name: Gabriele, Cembalo
// Date: 30 Oct 2025
//
// Code output:
// *****************************************************
// xlo = 0;       xhi = 0;              nint = 0;       nfv = 0
// xlo = 0;       xhi = 6;              nint = 12;      nfv = 36
// xlo = 0;       xhi = 3;              nint = 6;       nfv = 54
// xlo = 0;       xhi = 3.26345;        nint = 7;       nfv = 75
// xlo = 0;       xhi = 3.26347;        nint = 7;       nfv = 96
// Particle position (t=5) = 3.26347
// *****************************************************

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// definisco variabili globali
int ng = 3; // punti gaussiani
double tstar = 5.0; // tempo t^{*}
double xlo = 0.0;

double Time(double);
double TimeZero(double);
double Inv_Velocity(double);
int newton_method(double (*F)(double), double (*derF)(double), double, double, double, double, double &, int &);
double Gauss(double (*)(double), double, double, int, int);

int main(){

    // definisco la tolleranza
    double tol = 1.e-6;
    // definisco le iterazioni (che sono un parametro che inserisco all'interno del metodo di Newton)
    int l = 0;

    // definisco il punto che mi interessa trovare
    double xstar;

    // richiamo Newton per trovare lo zero della funzione t(x^*) - t^*
    double xhi = 6.0; // definisco e inizializzo l'estremo superiore

    // In Newton devo dare: Funzione, Derivata, estremi dell'intervallo, tolleranza su x e y, zero della funzione che mi interessa e le iterazioni
    newton_method(TimeZero, Inv_Velocity, xlo, xhi, tol, tol, xstar, l);

    // stampo a terminale il risultato di x^*
    cout << "Particle position (t=5) = " << xstar << endl;

    return 0;

}


double Time(double x){

    static int nfv = 0; // Cumulative number of function evaluations

    // definisco le costanti che mi servono
    double xhi = x; // gli estremi di integrazione (xlo è una variabile globale)

    // questo è il numero di intervalli
    int nint = ceil( 2.0*fabs( xlo - xhi ) );

    // faccio l'integrale \int_{xlo}^{xhi} dt/v(x) , richiamando la funzione inversa della velocità
    double integral = Gauss(Inv_Velocity, xlo, xhi, nint, ng);

    // modifico il numero di cumulativo
    nfv += ng*nint;

    // stampo a terminale
    cout << "xlo = " << xlo << ";       xhi = " << xhi << ";        nint = " << nint << ";       nfv = " << nfv << endl;

    // ritorno il valore di t(x) ( calcolo la funzione t(x^*)-t^* nella funzione TimeZero )
    return integral;

}

double TimeZero(double x){

    // ritorno il valore di t(x) - t^* che sarà la funzione di cui ricerco lo zero
    return Time(x) - tstar;

}

double Inv_Velocity(double x){

    return sqrt( 2 + sin(x)*sin(x)*sin(x)*sin(x) );

}

// metodo di Newton
// gli do in input la funzione, la derivata, gli estremi a e b, la tolleranza su x, la tolleranza su y, uno zero per riferimento e il numero di iterazioni
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

        // metto un controllo sugli intervalli, che mi blocca il ciclo se la nuova ampiezza deltax è maggiore di quella precedente. Mi blocca se cominciamo a non convergere
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

    // creo l'output delle iterazioni (esercizio froot.cpp)
    //cout << "(Newton) # = " << n << endl;

    l = n; // sono le iterazioni
    zero = x;
    return 0;

}

// implemento il metodo di Gauss in cui inserisco degli if per selezionare il grado del polinomio di Legendre voluto. Puoi vedere il capitolo 5 in cui calcoliamo zeri e pesi dei polinomi di Legendre di qualsiasi ordine.
// gli do in input la funzione, gli estremi di integrazione, il numero di intervalli e il numero di punti
double Gauss(double (*F)(double), double a, double b, int N, int ng){

    // Implemento la funzione per funzionare fino a 5 punti

    // definisco gli array di pesi e zeri del polinomio di Legendre
    double w[ng], x[ng];
    // definisco la variabile che mi somma l'integrale
    double sum = 0.0;

    // Riempio gli array in base al numero di punti che abbiamo (per questo esercizio è superfluo controllare che valore si ha di ng, ma ho preferito lasciarlo)
    if(ng == 1){
        x[0] = 0.;
        w[0] = 2.;
    }
    else if(ng == 2){
        x[1] = sqrt(1./3.);
        x[0] = -sqrt(1./3.);
        w[0] = 1.;
        w[1] = 1.;
    }
    else if(ng == 3){
        x[0] = -sqrt(3./5.);
        x[1] = 0.;
        x[2] = sqrt(3./5.);
        w[0] = 5./9.;
        w[1] = 8./9.;
        w[2] = 5./9.;
    }
    else if(ng == 4){
        x[0] = -sqrt( 3./7. - 2./7.*sqrt(6./5.) );
        x[1] = -sqrt( 3./7. + 2./7.*sqrt(6./5.) );
        x[2] = sqrt( 3./7. - 2./7.*sqrt(6./5.) );
        x[3] = sqrt( 3./7. + 2./7.*sqrt(6./5.) );
        w[0] = ( 18. + sqrt(30.) )/( 36. );
        w[1] = ( 18. - sqrt(30.) )/( 36. );
        w[2] = ( 18. + sqrt(30.) )/( 36. );
        w[3] = ( 18. - sqrt(30.) )/( 36. );
    }
    else if(ng == 5){
        x[0] = -1./3.*sqrt( 5. - 2.*sqrt(10./7.) );
        x[1] = -1./3.*sqrt( 5. + 2.*sqrt(10./7.) );
        x[2] = 0.;
        x[3] = 1./3.*sqrt( 5. - 2.*sqrt(10./7.) );
        x[4] = 1./3.*sqrt( 5. + 2.*sqrt(10./7.) );

        w[0] = ( 322. + 13.*sqrt(70.) )/( 900. );
        w[1] = ( 322. - 13.*sqrt(70.) )/( 900. ); 
        w[2] = 128./225.;
        w[3] = ( 322. + 13.*sqrt(70.) )/( 900. );
        w[4] = ( 322. - 13.*sqrt(70.) )/( 900. );
    }
    else{
        cout << "codice non implementato per Ng>5." << endl;
        return 0.;
    }

    // Calcolo l'integrale. Definisco l'ampiezza dell'intervallo e una somma ausiliaria per trovare l'integrale di Gauss
    double h = fabs(b-a)/(double)N;
    double sumj = 0.0;

    for( int i=0 ; i<N ; i++ ){

            // Devo fare il cambio di variabili che poi viene iterato sui vari intervalli successivi
            double x0 = a + i*h;
            double x1 = x0 + h;
            double s1 = (x1-x0)/2;
            double s2 = (x1+x0)/2;

            // azzero il conteggio di sumj
            sumj = 0.0;

            // calcolo l'integrale su ng punti di Gauss
            for( int j=0 ; j<ng ; j++ ){

                sumj += s1 * w[j] * F( s1*x[j] + s2 );

            }

            sum += sumj;

        }
    
    // Restituisco il risultato
    return sum;

}