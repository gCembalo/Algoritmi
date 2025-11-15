//
// Vedi le slide ch05 per vedere come organizzare il tuo codice.
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// includo l'header
#include "function_alg.h"

// definisco le macro
#define NMAX_EQ 64 // numero massimo di eq (lo inseriamo per sicurezza e per non
// lasciare array con dimensione variabile)
#define OMEGA 1.0 // pulsazione del problema (vedi harmonic.cpp Ch07)

// definizione di variabili globali. UTILIZZARE SOLO SE NECESSARIE
int g_LegendreN; // variabile che serve nel capitolo 5 per i polinomi di Legendre
double alpha = 10.0; // definita da "trajectory.cpp" (cap. 06)

// Dichiarazione delle funzioni in base alla lezione. Potrebbero esserci 
// funzioni definite più volte, quindi è sempre meglio controllare, 
// non solo il funzionamento, ma anche la sintassi con cui le si chiama


//------------------------------Per Gnuplot---------------------------------//

// se vuoi stampare in un file .dat dei fati devi includere la libreria:
//          #include <fstream>
// poi devi definire e aprire il file:
//          ofstream fdata; // declare Output stream class to operate on files
//          fdata.open("integral_sin.dat"); // open output file
// quando vuoi scriverci dentro devi fare:
//          fdata << x << "    " << gauss << endl;
// e una volta terminato devi chiudere il file:
//          fdata.close();
//
//
// conviene anche scrivere tutte le caratteristiche in un file .gp
// c'è un esempio.




//------------------------- 1-Introduzione --------------------------------//

// ------------------ practice2.cpp ------------------- //

// somma due numeri
double sum(double x, double y){
  return x+y;
}

// aggiunge 1 ad un numero dato
double addone(double x){
  return x + 1;
}

// calcola il quoziente e il resto di due numeri
int Quotient(int a, int b, int& q, int& r){
    if(b==0) return 1; //means failure
    q = a/b;
    r = a%b;
    return 0; //means success
}



//-------------------------- 2-Precision ---------------------------------//


// ------------------ quadratic.cpp ------------------- //

// funzione con cui prendo i coefficienti di una eq. di secondo grado
void datiQuad(double& x, double& y, double& z){

    // Prendo i coefficienti
    cout << "Abbiamo l'equazione:" << endl;
    cout << "    ax^2 + bx + c = 0      " << endl;
    cout << "Dammi a: ";
    cin >> x;
    cout << "Dammi b: ";
    cin >> y;
    cout << "Dammi c: ";
    cin >> z;

}

// funzione con formula standard per eq. secondo grado
void normQuad(const double& x, const double& y, const double& z, double& s1,
              double& s2){

    s1 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);
    s2 = ( -y - sqrt(y*y - 4*x*z) )/(2*x);

    cout << "La soluzione con (+) con il metodo standard: " << s1 << endl;
    cout << "La soluzione con (-) con il metodo standard: " << s2 << endl;

}

// funzione con formula alternativa per eq. secondo grado
void altQuad(const double& x, const double& y, const double& z, double& s3,
             double& s4){

    s3 = -(2*z)/( y + sqrt(y*y - 4*x*z) );
    s4 = -(2*z)/( y - sqrt(y*y - 4*x*z) );

    cout << "La soluzione con (+) con il metodo alternativa: " << s3 << endl;
    cout << "La soluzione con (-) con il metodo alternativa: " << s4 << endl;

}

// funzione che in base al segno del coefficiente b seleziona un metodo
// risolutivo piuttosto che l'altro
void signbQuad(const double& x, const double& y, const double& z, double& s1,
               double& s2){

    if(y>=0){

        s1  = ( -y - sqrt(y*y - 4*x*z) )/(2*x);
        s2 = -(2*z)/( y + sqrt(y*y - 4*x*z) );

    } else if (y<0){

        s1 = -(2*z)/( y - sqrt(y*y - 4*x*z) );
        s2 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);

    }

    cout << "+--------------------------------------------------------+" << endl;
    cout << "Si vede che a seconda del segno di b bisogna utilizzare una forma "
         << "o l'altra. Le soluzioni sono:" << endl;
    cout << s1 << endl;
    cout << s2 << endl;
    cout << "+--------------------------------------------------------+" << endl;

}



// ------------------ quadratic2.cpp ------------------- //

// richiedo le soluzioni di una certa eq di secondo grado
void solQuad(double& x1, double& x2){

    // Prendo i coefficienti
    cout << "Abbiamo l'equazione:" << endl;
    cout << "    ax^2 + bx + c = 0      " << endl;
    cout << "Dammi la soluzione 1: ";
    cin >> x1;
    cout << "Dammi la soluzione 2: ";
    cin >> x2;

}

// funzione che ordina le due input in modo crescente
void ordsolQuad(double& x1, double& x2){

    // Ordino le soluzioni
    double xtemp;
    if(x2 >= x1){
        xtemp = x1;
        x1 = x2;
        x2 = xtemp;
    }

}

// funzioni "norm" , "alt" e "singb" erano identiche a "quadratic.cpp" 
// ma senza l'output



// ------------------ roundoff.cpp ------------------- //


// Definizione funzione sqrt
double fx1r(const double& x){

    return sqrt(x*x + 1.) - x;

}

// Definizione funzione cos
double fx1c(const double& x){

    return 1. - cos(x);

}

// Definizione funzione sqrt razionalizzata
// moltiplico sopra e sotto per la funzione cambiata di segno
double fx2r(const double& x){

    return 1./(sqrt(x*x + 1.) + x);

}

// Definizione funzione cos razionalizzata
// moltiplico sopra e sotto per la funzione cambiata di segno
double fx2c(const double& x){

    return sin(x)/(1.+cos(x));

}

// Definizione sviluppo di Taylor della sqrt
double frTaylor(const double& x){

    return x + 1./(2.*x) - 1./(8.*x*x*x) - x;

}

// Definizione sviluppo di Taylor del cos
double fcTaylor(const double& x){

    return 1. - (1. - (x*x)/2. + (x*x*x*x)/(24.));

}



// ------------------ heron.cpp ------------------- //


// sviluppo il metodo di Herone per calcolare le radici quadrate
// gli do in input il numero di cui fare la radice, la guess, l'errore per 
// il controllo e il riferimento al risultato
double heron(const double &S, double &x, double &err, double &y){

    // sviluppo la funzione mettendo come controllo l'errore (il cui limite è
    // dato in input) sulla stima del metodo
    for( int i=0 ; i<100 ; i++ ){

        if(err > pow(10,-16)){

            y = 0.5*( x + S/x );
            err = fabs(x-y);
            x=y;

        }
        else{

            break;

        }

        cout << "Iterazione " << i << ": " << y << "; err = " << err << endl;

    }

    return x;

}



//-------------------------- 3-Quadrature ---------------------------------//


// ------------------ quadrature1.cpp ------------------- //

// Funzione che vogliamo integrare
double fExp(double x){

    return exp(-x);

}

// Regola rettangolo
// gli do in input la funzione, gli estremi di integrazione e il 
// numero di intervalli
double RectangularRule(double (*F)(double), double a, double b, int N){

    // Definisco l'ampiezza di ogni intervallo
    double h = fabs(a-b)/(double)N;

    // Definisco la variabile che somma i rettangoli
    double sum = 0.0;

    // Calcolo dell'integrale
    for(int i=0; i<=N-1; i++){

        sum += F(a + i*h) * h; // Integrale

    }

    return sum;

}

// Regola rettangolo ma con punto medio
// gli do in input la funzione, gli estremi di integrazione e il 
// numero di intervalli
double MidPointRule(double (*F)(double), double a, double b, int N){

    // Definisco l'ampiezza di ogni intervallo
    double h = fabs(a-b)/(double)N;

    // Definisco la variabile che somma i rettangoli
    double sum = 0.0;

    // Calcolo dell'integrale
    for(int i=0; i<=N-1; i++){

        sum += F(a + i*h + h/2) * h;

    }

    return sum;
}

// Regola trapezio
// gli do in input la funzione, gli estremi di integrazione e il 
// numero di intervalli
double TrapezoidalRule(double (*F)(double), double a, double b, int N){

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

// Regola Simpson estesa
// gli do in input la funzione, gli estremi di integrazione e il
// numero di intervalli
double ExtSimpsonRule(double (*F)(double), double a, double b, int N){

    int w = 4; // Definizione del peso
    // definisco le variabili di ampiezza degli intervalli e di somma
    double h , sum;

    h = fabs(a-b)/(double)N;
    sum = (F(a)+F(b))*h/3;

    // calcolo l'integrale
    for( int i=1 ; i<N ; i++ ){

        // vedi la formula sulle slide
        sum += w*F(a + i*h)*h/3;
        w = 6-w; // Calcolo il peso ogni volta

    }

    return sum;
}

// Funzione di convergenza per il metodo rettangolo
// gli do in input la funzione, gli estremi di integrazione e la tolleranza
double ConvergenceRectangular(double (*F)(double), double a, double b, double tol){

    // Controllo quanti intervalli servono
    int i = 2;

    // faccio il controllo (dimezzando l'intervallo) finché il risultato 
    // ottenuto non migliora meno della tolleranza
    while(fabs( RectangularRule(F,a,b,i) - RectangularRule(F,a,b,i*0.5) ) > tol){

        i = 2*i;

    }

    // stampo a terminale il risultato della regola insieme al numero 
    // delle iterazioni
    cout << "Rectangular: " << RectangularRule(F,a,b,i) << " iter: " << i << endl;

    return 0.;

}

// Funzione di convergenza per il metodo del punto medio
// gli do in input la funzione, gli estremi di integrazione e la tolleranza
double ConvergenceMidPoint(double (*F)(double), double a, double b, double tol){

    // Controllo quanti intervalli servono
    int i = 2;

    // faccio il controllo (dimezzando l'intervallo) finché il risultato
    // ottenuto non migliora meno della tolleranza
    while ( fabs( MidPointRule(F,a,b,i) - MidPointRule(F,a,b,i*0.5) ) > tol ){

        i = 2*i;

    }

    // stampo a terminale il risultato della regola insieme al numero 
    // delle iterazioni
    cout << "Mid Point: " << MidPointRule(F, a, b, i) << " iter: " << i << endl;

    return 0.;

}

// Funzione di convergenza per il metodo del trapezio
// gli do in input la funzione, gli estremi di integrazione e la tolleranza
double ConvergenceTrapezoidal(double (*F)(double), double a, double b, double tol){

    // Controllo quanti intervalli servono
    int i = 2;

    // faccio il controllo (dimezzando l'intervallo) finché il risultato ottenuto
    // non migliora meno della tolleranza
    while(fabs( TrapezoidalRule(F,a,b,i) - TrapezoidalRule(F,a,b,i*0.5) ) > tol){

        i = 2*i;

    }

    // stampo a terminale il risultato della regola insieme al numero 
    // delle iterazioni
    cout << "Trapezoidal: " << TrapezoidalRule(F, a, b, i) 
         << " iter: " << i << endl;

    return 0.;

}

// Funzione di convergenza per il metodo di Simpson
// gli do in input la funzione, gli estremi di integrazione e la tolleranza
double ConvergenceSimpson(double (*F)(double), double a, double b, double tol){

    // Controllo quanti intervalli servono
    int i = 2;

    // faccio il controllo (dimezzando l'intervallo) finché il risultato 
    // ottenuto non migliora meno della tolleranza
    while(fabs( ExtSimpsonRule(F,a,b,i) - ExtSimpsonRule(F,a,b,i*0.5) ) > tol){

        i = 2*i;

    }

    // stampo a terminale il risultato della regola insieme al numero 
    // delle iterazioni
    cout << "Simpson: " << ExtSimpsonRule(F, a, b, i) << " iter: " << i << endl;

    return 0.;

}



// ------------------ quadrature2.cpp ------------------- //

// definisco la funzione
double fSqrt(double x){

    return sqrt(1.+x);

}

// definisco la funzione
double func2(double x){

    return 1. - x + 2.*x*x + 0.5*x*x*x + x*x*x*x/4. - x*x*x*x*x/8.;

}

// il metodo di Simpson è identico a quello in "quadrature.cpp"



// implemento il metodo di Gauss in cui inserisco degli if per selezionare il 
// grado del polinomio di Legendre voluto. Puoi vedere il capitolo 5 in cui 
// calcoliamo zeri e pesi dei polinomi di Legendre di qualsiasi ordine.
// gli do in input la funzione, gli estremi di integrazione, il numero di 
// intervalli e il numero di punti
double Gauss(double (*F)(double), double a, double b, int N, int ng){

    // Implemento la funzione per funzionare fino a 5 punti

    // definisco gli array di pesi e zeri del polinomio di Legendre
    double w[NMAX_EQ], x[NMAX_EQ];
    // definisco la variabile che mi somma l'integrale
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

    // Calcolo l'integrale. Definisco l'ampiezza dell'intervallo e una somma
    // ausiliaria per trovare l'integrale di Gauss
    double h = fabs(b-a)/(double)N;
    double sumj = 0.0;

    for( int i=0 ; i<N ; i++ ){

            // Devo fare il cambio di variabili che poi viene iterato sui
            // vari intervalli successivi
            double x0 = a + i*h;
            double x1 = x0 + h;
            double s1 = (x1-x0)*0.5;
            double s2 = (x1+x0)*0.5;

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


// ------------------ integral_sin.cpp ------------------- //

// Funzioni
double funcSinx(double x){

    // Le singolarità vanno rimosse manualmente
    if( x == 0 ){

        return 1;

    }

    return sin(x)/x;
}

// funzioni Trapezoidal, ExtSimpson, Gauss identiche a quelle già fatte
// in "quadrature.cpp" e "quadrature2.cpp".



// ------------------ multid_quadrature.cpp ------------------- //

double func1(double x, double y){

    return x*x*x*x*y*y + 2*x*x*y*y - x*x*y + 2;

}

double funcCirc(double x, double y){

    // definisco la funzione definita a tratti
    if( sqrt(x*x + y*y) <= 1 ){

        return 1;

    }
    else{

        return 0;

    }

}

// gli do in input la funzione (con due variabili), gli estremi di integrazione,
// sia in x che in y, il numero di intervalli e di punti gaussiani
double Gauss2D(double (*F)(double, double), double x0, double x1, double y0, 
               double y1, int n, int ng){

    // Implemento la funzione per funzionare fino a 5 punti

    // definisco gli array di pesi e zeri del polinomio di Legendre
    double w[NMAX_EQ], x[NMAX_EQ];

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

    // Calcolo gli integrali

    // definisco l'ampiezza degli intervalli
    double dx = fabs(x0-x1)/(double)n;
    double dy = fabs(y0-y1)/(double)n;

    // definisco le grandezze ausiliarie per calcolare l'integrale
    double a , b , s1 , s2;
    double c , d , s3 , s4;

    // definisco una somma ausiliaria per trovare l'integrale di Gauss
    double sum = 0.0;

    for( int i = 0 ; i < n ; i++ ){  // Ciclo per le x

        // Definisco gli estremi dell'intervallo
        a = x0 + i*dx;
        b = a + dx;
        s1 = (b - a)*0.5;
        s2 = (b + a)*0.5;

        for( int j = 0 ; j < n ; j++ ){   // Ciclo per le y
            
            // Definisco gli estremi dell'intervallo
            c = y0 + j*dy;
            d = c + dy;
            s3 = (d - c)*0.5;
            s4 = (d + c)*0.5;

            for( int ik = 0 ; ik < ng ; ik ++ ){   // Ciclo per variare x

                for( int jk = 0 ; jk < ng ; jk ++){   // Ciclo per variare y

                    sum += w[ik]*s1*s3*w[jk]*F( s1*x[ik] + s2 , s3*x[jk] + s4 );

                }

            }

        }

    }

    return sum ;

}

// Funzione di convergenza per il metodo di Simpson
// gli do in input la funzione, gli estremi di integrazione e la tolleranza
void ConvergenceGauss(double (*F)(double, double), double x0, double x1, 
                      double y0, double y1, int ng, double tol){

    // Controllo quanti intervalli servono
    int i = 1;

    // faccio il controllo (incrementando l'intervallo) finché il risultato
    // ottenuto non migliora meno della tolleranza
    while ( fabs( Gauss2D(F, x0, x1, y0, y1, i, ng) - 3.1415926535897 ) > tol ){

        i ++ ;

    }

    cout << "Gauss = " << Gauss2D(F, x0, x1, y0, y1, i, ng) 
         << "\ncon " << i << " intervalli." << endl;

}




//---------------------------- 4-Random -----------------------------------//

// ------------------ guess.cpp ------------------- //
// no function...

// ------------------ prn_uniformity.cpp ------------------- //
// no function...

// ------------------ decay.cpp ------------------- //
// no function...

// ------------------ gauss_distrib.cpp ------------------- //

double gaussiana(const double &x, const double &sigma){

    return 1./(sigma*sqrt(2.*3.1415926)) * exp( -0.5 * (x*x) / (sigma*sigma) );

}

// ------------------ pi.cpp ------------------- //
// era definita solo la gaussiana già presente in "gauss_distrib.cpp"



//------------------------- 5-Root_finders --------------------------------//


// ------------ root finders function ------------- //

// metodo della bisezione
// gli do in input la funzione, gli estremi a e b, la tolleranza su x, 
// uno zero per riferimento e il numero di iterazioni
int bisection(double (*F)(double), double a, double b, double tol, 
              double &zero, int &l){

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
        cout << "(Bisection) # = " << n 
             << "    (l'estremo " << a << " è già lo zero)"<< endl;

        l = n;
        return 0;

    }
    else if( fb == 0.0 ){

        zero = b;

        // creo l'output delle iterazioni
        cout << "(Bisection) # = " << n 
             << "(l'estremo " << b << " è già lo zero)"<< endl;

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
            x = ( a+b ) * 0.5;

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
            //cout << "n = " << n << ";   [a,b] = [" << a << ", " << b 
            //       << "];    xm = " << x << ";   Deltax = " << fabs(a-b) 
            //       << ";   f(xm) = " << F(x) << endl;

        }

        // creo l'output delle iterazioni
        //cout << "(Bisection) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = x;
        return 0;

    }

}

// metodo della false position
// gli do in input la funzione, gli estremi a e b, la tolleranza su x, 
// uno zero per riferimento e il numero di iterazioni
int false_position(double (*F)(double), double a, double b, double tol, 
                   double &zero, int &l){

    double x = 3; // la guess di zero della funzione
    double xk = 0; // variabile che mi serve per la tolleranza
    double m, q; // parametri della retta
    int n = 0; // variabile per contare

    // definisco le variabili della funzione valutata
    double fa = F(a);
    double fb = F(b);
    double fx;

    // metto i controlli di non avere già uno zero
    if( fa == 0.0 ){

        zero = a;

        // creo l'output delle iterazioni
        cout << "(False position) # = " << n 
             << "    (l'estremo " << a << " è già lo zero)"<< endl;

        l = n;
        return 0;

    }
    else if( fb == 0.0 ){

        zero = b;

        // creo l'output delle iterazioni
        cout << "(False position) # = " << n 
             << "(l'estremo " << b << " è già lo zero)"<< endl;
        
        l = n;
        return 0;

    }
    else{
        // metto nel ciclo la condizione sia sulla tolleranza che sul 
        // numero di cicli (inserita dopo)
        while( fabs( x - xk ) > tol ){
        
            n++;
            fa = F(a);
            fb = F(b);

            if( n == 100 ){

                cout << "(False position) Troppe iterazioni." << endl;

                l = n;
                return 0;

            }

            xk = x;
        
            // trovo la retta
            m = ( fa - fb ) / (a - b);
            q = ( fb*a - fa*b ) / (a - b);
            // trovo lo zero della retta e lo chiamo x
            x = - q / m;

            fx = F(x);

            // controllo dove si trova x rispetto gli estremi a e b
            if( fa*fx < 0 ){

                b = x;

            }
            else if ( fa*fx > 0 ){

                a = x;

            }

            // creo l'output voluto (esercizio froot.cpp)
            //cout << "n = " << n << ";   [a,b] = [" << a << ", " << b 
            //       << "];    xm = " << x << ";   Deltax = " << fabs(a-b) 
            //       << ";   f(xm) = " << F(x) << endl;

        }

        // creo l'output delle iterazioni
        //cout << "(False position) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = x;
        return 0;

    }
}

// metodo della secante
// gli do in input la funzione, gli estremi a e b, la tolleranza su x, 
// uno zero per riferimento e il numero di iterazioni
int secant_method(double (*F)(double), double a, double b, double tol, 
                  double &zero, int &l){

    // definisco le variabili che mi servono per tenere traccia delle
    // varie iterazioni di x
    double xk1 = a, xk = b, xk2 = xk + 1; // dove uso xk come x_k, xk1 come
    // x_{k-1} e xk2 come x_{k+1} ; inizializzo gli zeri sugli estremi 
    // dell'intervallo in cui ricaviamo la retta
    int n = 0; // variabile per contare
    // una variabile di controllo per vedere di quanto miglioriamo la guess
    double xp = 0;

    // definisco le variabili della funzione valutata
    double fa = F(a);
    double fb = F(b);
    double fxk;
    double fxk1;

    // metto i controlli di non avere già uno zero
    if( fa == 0.0 ){

        zero = a;

        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n 
             << "    (l'estremo " << a << " è già lo zero)" << endl;

        l = n; // sono le iterazioni
        return 0;

    }
    else if( fb == 0.0 ){

        zero = b;

        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n 
             << "(l'estremo " << b << " è già lo zero)" << endl;

        l = n; // sono le iterazioni
        return 0;

    }
    else{
        
        // metto nel ciclo la condizione sia sulla tolleranza che sul numero 
        // di cicli (messo dopo)
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
            //cout << "n = " << n << ";   [a,b] = [" << xk1 << ", " << xk 
            //       << "];    x0 = " << xk2 << ";   Deltax = " << fabs(xk2-xk1) 
            //       << ";   f(x0) = " << F(xk2) << endl;

        }

        // creo l'output delle iterazioni
        //cout << "(Secant) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = xk2;
        return 0;

    }

}

// metodo di Newton
// gli do in input la funzione, la derivata, gli estremi a e b, la tolleranza
// su x, la tolleranza su y, uno zero per riferimento e il numero di iterazioni
int newton_method(double (*F)(double), double (*derF)(double), double a, 
                  double b, double xtol, double ytol, double &zero, int &l){

    int n = 0; // variabile per contare

    // controllo che gli estremi non siano già degli zeri
    double fa = F(a), fb = F(b);

    if( fabs(fa) < ytol ){

        zero = a;

        // creo l'output delle iterazioni
        cout << "(Newton) # = " << n 
             << "    (l'estremo " << a << " è già lo zero)" << endl;

        return 0;

    }

    if( fabs(fb) < ytol ){
        
        zero = b;
        // creo l'output delle iterazioni
        cout << "(Newton) # = " << n 
             << "    (l'estremo " << b << " è già lo zero)" << endl;

        return 0;

    }

    // definisco delle variabili utili all'algoritmo
    double Deltax = fabs( b-a ); // l'ampiezza dell'intervallo inziale
    double deltax = Deltax * 0.5; // la nuova ampiezza (inizializzata)
    // la mia prima guess, che prendo a metà dell'intervallo
    double x = ( a + b )*0.5;
    double fx = F(x); // funzione valutata nella guess

    // faccio il ciclo per il metodo di Newton, con i controlli su x, su y e 
    // sul numero di iterazioni (nel mezzo del while)
    while( fabs(deltax) > xtol && fabs(fx) > ytol ){

        n++;

        // blocco il ciclo se raggiungo troppe iterazioni
        if( n == 100 ){

            cout << "(Newton) Troppe iterazioni. (" << l << ")" << endl;
            l = n; // sono le iterazioni
            return 0;

        }

        // metto un controllo sugli intervalli, che mi blocca il ciclo se la 
        // nuova ampiezza deltax è maggiore di quella precedente. Mi blocca se 
        // cominciamo a non convergere
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

// bracketing function
// gli do in input la funzione, gli estremi a e b, il numero di intervalli, 
// un array per gli estremi (sinistro e destro) degli intervalli, e un 
// riferimento al numero di zeri
int Bracket(double (*F)(double), double a, double b, double n, double *xL, 
            double *xR, int &nroots){

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

        // metto la condizione in cui abbiamo un cambiamento di segno e quindi
        // con la quale ci ricordiamo il valore
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



// ------------------ froot.cpp ------------------- //

double funcFroot(double &x){

    return exp(-x) - x;

}

double derfuncFroot(double &x){

    return -exp(-x) - 1;

}


// ------------------ horner.cpp ------------------- //

// funzione polinomio
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

// funzione derivata del polinomio
double derpol(double x){

    // definisco l'array e il polinomio (derivata)
    double a[3] = {1, -6, 3}; // coefficienti della derivata
    //double a[4] = {5, 1, -3, 1}; // i coefficienti del polinomio
    // inizializzo der_p come a_n
    double dp = a[2];

    // valuto il polinomio con il metodo di horner
    for( int j = 2-1 ; j >= 0 ; j-- ){

        // moltiplico per x il termine in p
        dp = a[j] + dp*x;
        
    }

    return dp;
}


// ------------------ session4.cpp ------------------- //

double funcSes4(double x){

    return exp( 1./( x + 0.5 ) ) - ( 3 + 2*x )/( 1 + x );

}

double derfuncSes4(double x){

    return -1./((x + 0.5)*(x+0.5)) * exp( 1./(x + 0.5)) + 1./((1+x) * (1+x));

}

// ------------------ froot.cpp ------------------- //

double funcSin(double x){
    
    return sin(x) - ( x*x/100. + x/5. + 1./3. );

}

double derfuncSin(double x){
    
    return cos(x) - ( x/50. + 1./5. );

}

// ------------------ legendre.cpp ------------------- //

// polinomio di Legendre
double polLegendre(double x){

    if (g_LegendreN == 0){
        return 1.;
    }

    if (g_LegendreN == 1){
        return x;
    }

    // definisco i polinomi
    double P0 = 1., Pi = x, Pi1;  // P0 è P_0 ; Pi è P_n ; Pi1 è P_{n+1}

    for( int i = 1 ; i < g_LegendreN ; i++ ){

        Pi1 = ( (2.0 * i + 1.0) * x * Pi - i * P0 ) / ( i+1.0 );

        // Aggiorno i valori
        P0 = Pi;
        Pi = Pi1;
       
    }

    return Pi1;

}

// derivata polinomio di Legendre
double derpolLegendre(double x){

    if (g_LegendreN == 0.0){
        return 0.0;
    }

    if (g_LegendreN == 1.0){
        return 1.0;
    }

    double P0 = 1.0; // P0(x)
    double P1 = x; // P1(x)
    double Pi1; // P_{n+1}

    // stesso ciclo della funzione del polinomio
    for( int i = 1 ; i < g_LegendreN ; i++ ){

        Pi1 = ( (2.0 * i + 1.0) * x * P1 - i * P0 ) / ( i + 1.0 );

        // aggiorno i valori
        P0 = P1;
        P1 = Pi1;

    }

    // uso la formula per la derivata
    double dPi = ( g_LegendreN * (x * P1 - P0) ) / ( x * x - 1.0 );

    return dPi;
    
}

// pesi dei polinomi di Legendre
double wi(double x){

    return 2.0/( ( 1.0-x*x )*( derpolLegendre(x)*derpolLegendre(x) ) );

}



//-------------------------- 6-Derivative ---------------------------------//

// ------------------ derivative.cpp ------------------- //

double sinDer(double x){

    return sin(x);

}

double sinDerEx(double x){

    return cos(x);

}

// x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}
double derFD(double (*F)(double), double xi, double xp, double h){

    // definisco 
    double fiPrime = ( F(xp) - F(xi) )/h;

    return fiPrime;

}

// x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}
double derBD(double (*F)(double), double xi, double xm, double h){

    // definisco 
    double fiPrime = ( F(xi) - F(xm) )/h;

    return fiPrime;

}

// x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}
double derCD(double (*F)(double), double xm, double xp, double h){

    // definisco 
    double fiPrime = ( F(xp) - F(xm) )/( 2.0*h );

    return fiPrime;

}

// x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}
double der4th(double (*F)(double), double xmm, double xm, double xp, double xpp,
              double h){

    // definisco 
    double fiPrime = ( F(xmm) - 8.0*F(xm) + 8.0*F(xp) - F(xpp) )/( 12.0*h );

    return fiPrime;

}



// ------------------ trajectory.cpp ------------------- //

double position(double t){

    // elimino la singolarità
    if( t == 0 ){
        return 0.;
    }
    else{
        return alpha*t*t - t*t*t*( 1 - exp( -alpha*alpha/t ) );
    }

}

// definisco la velocità che mi serve solo per la derivata seconda in 0 avendo
// il problema di due limiti diversi
double velocity(double t){

    // elimino la singolarità
    if( t == 0 ){
        return 0.;
    }
    else{
        return 2*alpha*t - 3*t*t*( 1 - exp( -alpha*alpha/t ) ) + 
               ( alpha*alpha*t )*exp( -alpha*alpha/t ) ;
    }

}

// x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}
double SecondDerivative(double (*F)(double x), double xi, double xm, double xp, 
                        double h){

    // definisco
    double fprimeprime = ( F(xp) - 2.0*F(xi) + F(xm) )/( h*h );

    return fprimeprime;
    
}



//------------------------------ 7-ODE -------------------------------------//


// programmi utilizzati per tutte le esercitazioni

// implemento il metodo di Eulero. Gli do in input la variabile di integrazione,
// un puntatore alle soluzioni, un puntatore al Right-Hand-Side-Function, 
// lo step da utilizzare e la dimensionalità di Y, che non è altro che il numero
// di ODE di primo ordine che abbiamo.
void EulerStep(double t, double *Y, void (*RHSFunc)(double, double *, double *),
               double dt, int neq){
    
    // dt è lo step che utilizziamo per trovare la soluzione di dY/dt = rhs.
    // neq è il numero di ODE (dimensionalità di Y[])
    // *RHSFunc() punta al Right-Hand-Side-Function (in questo caso dYdt())
    
    int k; // variabile per visitare tutte le componenti di *Y
    double rhs[NMAX_EQ]; // per assicurarsi che rhs[] sia grande abbastanza
                         // usiamo NMAX_EQ (neq < 256)

    // calcolo il lato destro dell'equazione
    RHSFunc (t, Y, rhs);

    for( k = 0 ; k < neq ; k++ ){

        Y[k] += dt*rhs[k]; // Update solution array

    }

}

// implemento il metodo Runge-Kutta del secondo ordine (midpoint).
// gli do in input la variabile di integrazione, il puntatore alle soluzioni, 
// il puntatore alla funzione del Right-Hand-Side-Function, l'incremento e 
// l'ordine della ODE.
void RK2StepMid(double t, double *Y, void (*RHSFunc)(double t, double *Y, 
                double *R), double h, int neq){
    
    // definisco i vettori per gli step intermedi
    double Y1[NMAX_EQ], k1[NMAX_EQ], k2[NMAX_EQ];
    
    RHSFunc(t,Y,k1); // calcolo k1 con il RSH con t_n e Y_n

    // scrivo il ciclo per determinare Y_n + k1*h/2
    for( int i = 0 ; i < neq ; i++ ){
        
        Y1[i] = Y[i] + 0.5*h*k1[i];

    }

    RHSFunc(t+0.5*h,Y1,k2); // calcolo k2 con il RSH con t_n+h/2 e Y_n+k1*h/2
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + k2*h
    for (int i = 0 ; i < neq ; i++){
        
        Y[i] += h*k2[i];

    }

}

// implemento il metodo Runge-Kutta del secondo ordine (modified Eulero).
// gli do in input la variabile di integrazione, il puntatore alle soluzioni, 
// il puntatore alla funzione del Right-Hand-Side-Function, l'incremento e 
// l'ordine della ODE.
void RK2StepHeun(double t, double *Y, void (*RHSFunc)(double t, double *Y, 
                 double *R), double h, int neq){
    
    // definisco i vettori per gli step intermedi
    double Y1[NMAX_EQ], k1[NMAX_EQ], k2[NMAX_EQ];
    
    RHSFunc(t,Y,k1); // calcolo k1 con il RSH con t_n e Y_n

    // scrivo il ciclo per determinare Y_n + k1*h
    for( int i = 0 ; i < neq ; i++ ){
        
        Y1[i] = Y[i] + h*k1[i];

    }

    RHSFunc(t+h,Y1,k2); // calcolo k2 con il RSH con t_n+h e Y_n+k1*h
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + (k1+k2)*h/2
    for (int i = 0 ; i < neq ; i++){
        
        Y[i] += ( k1[i] + k2[i] )*h*0.5;

    }

}

// implemento il metodo Runge-Kutta del quarto ordine.
// gli do in input la variabile di integrazione, il puntatore alle soluzioni,
// il puntatore alla funzione del Right-Hand-Side-Function, l'incremento e
// l'ordine della ODE.
void RK4Step(double t, double *Y, void (*RHSFunc)(double t, double *Y, double *R),
             double h, int neq){
    
    // definisco i vettori per gli step intermedi
    double Y1[NMAX_EQ], k1[NMAX_EQ], k2[NMAX_EQ], k3[NMAX_EQ], k4[NMAX_EQ];
    
    RHSFunc(t,Y,k1); // calcolo k1 con il RSH con t_n e Y_n

    // scrivo il ciclo per determinare Y_n + k1*h/2
    for( int i = 0 ; i < neq ; i++ ){
        
        Y1[i] = Y[i] + 0.5*h*k1[i];

    }

    RHSFunc(t+0.5*h,Y1,k2); // calcolo k2 con il RSH con t_n+h/2 e Y_n+k1*h/2
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + k2*h/2
    for (int i = 0 ; i < neq ; i++){
        
        Y1[i] = Y[i] + h*k2[i]*0.5;

    }

    RHSFunc(t+0.5*h,Y1,k3); // calcolo k3 con il RSH con t_n+h/2 e Y_n+k2*h/2
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + k3*h
    for (int i = 0 ; i < neq ; i++){
        
        Y1[i] = Y[i] + h*k3[i];

    }

    RHSFunc(t+h,Y1,k4); // calcolo k4 con il RSH con t_n+h e Y_n+k3*h
    
    // scrivo il ciclo per calcolare 
    // Y_{n+1} = Y_n + h/6 * ( k1 + 2*k2 + 2*k3 + k4 )
    for (int i = 0 ; i < neq ; i++){
        
        Y[i] += h * ( k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i] ) / 6.0;

    }

}

// implemento il metodo di Verlet per la posizione
// gli do in input i puntatori ai vettori posizione e velocità, la dimensione
// degli array (ordine della EDO), l'incremento degli intervalli e la funzione 
// accelerazione
void PositionVerletStep(double *x, double *v, int neq, double dt,
                        void (*acceleration)(double *, double *, int)){

    double a[NMAX_EQ]; // creo il vettore per l'accelerazione

    // calcolo x per mezzi step
    for( int i = 0 ; i < neq ; i++ ){

        x[i] += 0.5 * dt * v[i];

    }

    // calcolo l'accelerazione a t = i*t + t/2
    acceleration( x, a, neq );

    // calcolo la velocità per step interi
    for( int i = 0 ; i < neq ; i++ ){

        v[i] += dt * a[i];

    }

    // calcolo la posizione per step intero
    for( int i = 0 ; i < neq ; i++ ){

        x[i] += 0.5 * dt * v[i];

    }

}

// implemento il metodo di Verlet per la posizione (alternativo)
// questo metodo è stato proposto dal professore come una versione che
// riutilizza l'approccio standard del RK-method, richiamando "RHSFunc" come
// "Acceleration".
// In questo caso neq è il numero totale di equazioni.
// gli do in input 
void PositionVerletStep2(double t, double *Y, void (*RHSFunc)(double,
                         double *, double *), double dt, int neq){
    
    int n;
    double dYdt[NMAX_EQ];
    // nota che nelle righe che seguono è corretto dividere per un intero poiché
    // in aritmetica dei puntatori si deve usare solo un intero.
    double *v = Y + neq/2; // simple pointer arithmetic
    double *x = Y;
    double *a = dYdt + neq/2;

    // 1. Evolve positions by half a step [drift]
    for( n = 0 ; n < neq/2 ; n++ ){
        
        x[n] += 0.5 * dt * v[n];
    
    }

    // 2. Compute acceleration at t= t(n+1/2)
    RHSFunc(t+dt/2.0, Y, dYdt);

    // 3. Evolve velocities for full step
    for( n = 0 ; n < neq/2 ; n++ ){
        
        v[n] += dt * a[n];
    
    }

    // 4. Evolve positions by half a step
    for( n = 0 ; n < neq/2.0 ; n++ ){
        
        x[n] += 0.5 * dt * v[n];

    }

}

// implemento il metodo di Verlet per la posizione
// gli do in input i puntatori ai vettori posizione e velocità, la dimensione
// degli array (ordine della EDO), l'incremento degli intervalli e la funzione 
// accelerazione
void VelocityVerletStep(double *x, double *v, int neq, double dt,
                        void (*acceleration)(double *, double *, int)){

    double a[NMAX_EQ]; // creo il vettore per l'accelerazione

    // calcolo l'accelerazione a t = i*t
    acceleration( x, a, neq );

    // calcolo v per mezzi step
    for( int i = 0 ; i < neq ; i++ ){

        v[i] += 0.5 * dt * a[i];

    }

    // calcolo la posizione  per step interi
    for( int i = 0 ; i < neq ; i++ ){

        x[i] += dt * v[i];

    }

    // calcolo l'accelerazione per x^{n+1}
    acceleration( x, a, neq );

    // calcolo v per step intero
    for( int i = 0 ; i < neq ; i++ ){

        v[i] += 0.5 * dt * a[i];

    }
    
}


// ------------------ ode1.cpp ------------------- //

// funzione problem dependent
void RHSFuncOde1(double t, double *Y, double *R){

    // Compute the right-hand side of the ODE dy/dt = -t*y
    double y = Y[0];
    R[0] = -t*y;

}

// soluzione esatta
double ode1Sol(double t){

    return exp( -t*t/2 );

}

// ------------------ ode2.cpp ------------------- //

double ode2Sol(double t){

    return cos(t);

}
// definisco il Right-Hand-Side-Function (è problem dependent). 
// Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore) 
// mi faccio dare R
void RHSFuncOde2(double t, double *Y, double *R){

    // Compute the right-hand side of the ODE (2 equation)
    double x = Y[0];
    double y = Y[1];
    R[0] = y;
    R[1] = -x;

}

// ------------------ kepler.cpp ------------------- //

// definisco il Right-Hand-Side-Function (è problem dependent). 
// Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore) mi 
// faccio dare R
void RHSFuncOde3(double t, double *Y, double *R){

    // setto le condizioni iniziali
    double x = Y[0];
    double y = Y[1];
    double vx = Y[2];
    double vy = Y[3];

    double r = sqrt( x*x + y*y );

    R[0] = vx;
    R[1] = vy;
    R[2] = -x/( r*r*r );
    R[3] = -y/( r*r*r );

}

// ------------------ harmonic.cpp ------------------- //

// definisco il Right-Hand-Side-Function (è problem dependent). 
// Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore)
// mi faccio dare R
void RHSFuncOde4(double t, double *Y, double *R){

    // Compute the right-hand side of the ODE (2 equation)
    double x = Y[0];
    double y = Y[1];

    R[0] = y;
    R[1] = - OMEGA * OMEGA * x;

}

// implemento la funzione accelerazione che mi serve per il metodo di Verlet
// gli do in input gli array di posizione e accelerazione, oltre ad il numero di
// punti
void acceleration(double *x, double *a, int n){

    for( int i = 0 ; i < n ; i++ ){

        a[i] = - OMEGA * OMEGA * x[i];
    }

}