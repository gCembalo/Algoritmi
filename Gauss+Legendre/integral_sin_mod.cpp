// write a program to compute
// int_0^x sinx/x dx
// at x = 0.8 using intervals h = 0.8, 0.4, 0.2, 0.1. The correct value, to ten decimals, is Si(0.8) = 0.77209 57855.
//
// Next, using gnuplot, produce a plot like the one in the figure for 0 < x < 25. Try to not sacrifice accuracy as x increases !
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;
ofstream fdata;

// Funzione
double func(double);
// Diversi metodi
double Trapezoidal(double (*)(double),double,double,int);
double Simpson(double (*)(double),double,double,int);
double Gauss(double (*)(double),double,double,int,int);

int newton_method(double (*F)(double), double (*derF)(double), double, double, double, double, double &, int &);
int Bracket(double (*F)(double), double a, double b, double n, double *xL, double *xR, int &nroots);
double polLegendre(double);
double derpolLegendre(double);
double wi(double);

int ng;

int main(){
    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    fdata << setprecision(7);  //set the significant digits
    fdata << setiosflags ( ios::scientific );  //set the scientific notation

    double x0 = 0.; // L'estremo di integrazione
    double x = 0.8; // Il nostro punto
    // I diversi intervalli li seleziono dando diversi valori di N, così che si divida (x-x0) in N parti e si ottengano effettivamente 0.8, 0.4, 0.2, 0.1

    cout << "\n+-----------------------------------+" << endl;
    cout << "Che grado del polinomio di Legendre vuoi?" << endl;
    cin >> ng;

    cout << "\n+-----------------------------------+" << endl;
    cout << "\nValore esatto Si(0.8) = 0.77209 57855\n" << endl;
    cout << "Valutazione con i diversi metodi:\n" << endl;
    cout << "Gauss (n = 1, ng = " << ng << ") = " << Gauss(func,x0,x,1,ng) << endl;
    cout << "Trapezoidal (n = 1)   = " << Trapezoidal(func,x0,x,1) << endl;
    cout << "Trapezoidal (n = 2)   = " << Trapezoidal(func,x0,x,2) << endl;
    cout << "Trapezoidal (n = 4)   = " << Trapezoidal(func,x0,x,4) << endl;
    cout << "Trapezoidal (n = 8)   = " << Trapezoidal(func,x0,x,8) << endl;
    cout << "Simpson     (n = 2)   = " << Simpson(func,x0,x,2) << endl;
    cout << "Simpson     (n = 4)   = " << Simpson(func,x0,x,4) << endl;
    cout << "Simpson     (n = 8)   = " << Simpson(func,x0,x,8) << endl;


    // Implemento per stampare la funzione 0<x<25
    ofstream fdata; // declare Output stream class to operate on files
    fdata.open("integral_sin.dat"); // open output file

    double gauss = 0.0;
    for( double x=0.0; x<=25.; x += 0.1 ){
        // Scrivo il file
        gauss += Gauss(func,x0,x,1,ng);
        fdata << x << "    " << gauss << endl;
        x0 = x;
    }
    fdata.close();  // Chiudo il file

    cout << "\n+-----------------------------------+\n" << endl;
    double asint = 0.0;
    double x1 = 0.0;
    for( double i = 0.0 ; i<=1000 ; i += 0.1 ){
        // Scrivo il file
        asint += Gauss(func,x1,i,1,ng);
        x1 = i;
    }
    cout << "Puoi vedere il plot di Si(x) su gnuplot.\nL'asintoto è = " << asint << endl;
    cout << "\n+-----------------------------------+\n" << endl;

    return 0;
}



// Funzioni
double func(double x){
    if( x == 0 ){ // Le singolarità vanno rimosse manualmente
        return 1;
    }
    return sin(x)/x;
}

double Trapezoidal(double (*F)(double), double a, double b, int N){
    // Definisco l'ampiezza di ogni intervallo
    double h = fabs(a-b)/(double)N;
    // Definisco la variabile che somma i rettangoli
    double sum = 0.0;
    // Calcolo dell'integrale
    for(int i=0; i<=N-1; i++){
        sum += 0.5 * F(a + i*h) * h + 0.5 * F(a + (i+1)*h) * h;
    }
    return sum;
}

double Simpson(double (*F)(double), double a, double b, int N){
    int w = 4; // Definizione del peso
    double h, sum;
    // Definisco l'ampiezza di ogni intervallo
    h = fabs(a-b)/(double)N;
    sum = (F(a)+F(b))*h/3.;
    for (int i=1; i<N; i++ )
    {
        sum += w*F(a + i*h)*h/3.;
        w = 6.-w; // Calcolo il peso ogni volta
    }
    return sum;
}

/*
double Gauss(double (*F)(double), double a, double b, int N, int ng){
    // Implemento la funzione per funzionare fino a 5 punti
    double w[ng], x[ng]; // array di pesi e zeri del polinomio di Legendre
    double sum = 0.0;

    // Riempio gli array in base al numero di punti che abbiamo
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

    // Calcolo l'integrale
    double h = fabs(b-a)/(double)N;
    double sumj = 0.0;

    for(int i=0; i<N; i++){
            // Devo fare il cambio di variabili che poi viene iterato sui vari intervalli successivi
            double x0 = a + i*h;
            double x1 = x0 + h;
            double s1 = (x1-x0)/2;
            double s2 = (x1+x0)/2;
            // azzero il conteggio di sumj
            sumj = 0.0;
            for(int j=0; j<ng; j++){
                sumj += s1 * w[j] * F( s1*x[j] + s2 );
            }
            sum += sumj;
        } 
    
    // Restituisco il risultato
    return sum;
}
    */




// metodo di Gauss in cui calcolo i plinomi di Legendre di grado qualsiasi
double Gauss(double (*F)(double), double a, double b, int N, int ng){

    // definisco gli arrai degli estremi
    double xL[ng], xR[ng];

    // ng (il grado del polinomio) sarà il numero di root che avremo nel polinomio e lo possiamo trovare con il braket
    Bracket(polLegendre, a, b, 2*ng+1, xL, xR, ng);

    // utilizzo il metodo di Newton per trovare gli zeri
    double w[ng], x[ng]; // array di pesi e zeri del polinomio di Legendre
    int l = 0; // iterazioni
    double tol = 1.e-8;

    // gli do in pasto polinomi e sue derivate, gli estremi degli intervalli, una tolleranza di convergenza (sia x che y), un vettore in cui mette gli zeri, e un numero di iterazioni
    for( int k = 0 ; k < ng ; k++ ){

        // riempio il vettore di zeri
        newton_method(polLegendre, derpolLegendre, xL[k], xR[k], tol, tol, x[k], l);

        // riempio il vettore dei pesi
        w[k] = wi(x[k]);

    }



    // ora possiamo calcolare l'integrale con il metodo di Gauss
    // definisco la variabile che mi restituisce l'integrale
    double sum = 0.0;

    // Calcolo l'integrale
    double h = fabs(b-a)/(double)N;
    double sumj = 0.0; // somma ausiliaria che uso dentro il ciclo

    for( int i=0 ; i<N ; i++ ){

            // Devo fare il cambio di variabili che poi viene iterato sui vari intervalli successivi
            double x0 = a + i*h;
            double x1 = x0 + h;
            double s1 = (x1-x0)/2;
            double s2 = (x1+x0)/2;

            // azzero il conteggio di sumj
            sumj = 0.0;

            for(int j=0; j<ng; j++){

                sumj += s1 * w[j] * F( s1*x[j] + s2 );

            }

            sum += sumj;

        } 
    
    // Restituisco il risultato
    return sum;
}








//                  Polinomi di Legendre

double polLegendre(double x){

    if (ng == 0){
        return 1.;
    }

    if (ng == 1){
        return x;
    }

    // definisco i polinomi 0 e 1
    double P0 = 1., Pi = x, Pi1;  // P0 è P_0 ; Pi è P_n ; Pi1 è P_{n+1}

    for( int i = 1 ; i < ng ; i++ ){

        Pi1 = ( (2.0 * i + 1.0) * x * Pi - i * P0 ) / ( i+1.0 );

        // Aggiorno i valori
        P0 = Pi;
        Pi = Pi1;
       
    }

    return Pi1;

}


double derpolLegendre(double x){

    if (ng == 0.0){
        return 0.0;
    }

    if (ng == 1.0){
        return 1.0;
    }

    double P0 = 1.0; // P0(x)
    double P1 = x; // P1(x)
    double Pi1; // P_{n+1}

    // stesso ciclo della funzione del polinomio
    for( int i = 1 ; i < ng ; i++ ){

        Pi1 = ( (2.0 * i + 1.0) * x * P1 - i * P0 ) / ( i + 1.0 );

        // aggiorno i valori
        P0 = P1;
        P1 = Pi1;

    }

    // uso la formula per la derivata
    double dPi = ( ng * (x * P1 - P0) ) / ( x * x - 1.0 );

    return dPi;
    
}

double wi(double x){

    return 2.0/( ( 1.0-x*x )*( derpolLegendre(x)*derpolLegendre(x) ) );

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