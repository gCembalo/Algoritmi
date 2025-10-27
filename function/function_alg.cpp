//
// Vedi le slide ch05 per vedere come organizzare il tuo codice.
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// includo l'header
#include "function_alg.h"

// Dichiarazione delle funzioni in base alla lezione. Potrebbero esserci funzioni definite più volte, quindi è sempre meglio controllare, non solo il funzionamento, ma anche la sintassi con cui le si chiama


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




//-------------------------------------------------- 1-Introduzione ---------------------------------------------------------//

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



//-------------------------------------------------- 2-Precision ---------------------------------------------------------//


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
void normQuad(const double& x, const double& y, const double& z, double& s1, double& s2){

    s1 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);
    s2 = ( -y - sqrt(y*y - 4*x*z) )/(2*x);

    cout << "La soluzione con (+) con il metodo standard: " << s1 << endl;
    cout << "La soluzione con (-) con il metodo standard: " << s2 << endl;

}

// funzione con formula alternativa per eq. secondo grado
void altQuad(const double& x, const double& y, const double& z, double& s3, double& s4){

    s3 = -(2*z)/( y + sqrt(y*y - 4*x*z) );
    s4 = -(2*z)/( y - sqrt(y*y - 4*x*z) );

    cout << "La soluzione con (+) con il metodo alternativa: " << s3 << endl;
    cout << "La soluzione con (-) con il metodo alternativa: " << s4 << endl;

}

// funzione che in base al segno del coefficiente b seleziona un metodo risolutivo piuttosto che l'altro
void signbQuad(const double& x, const double& y, const double& z, double& s1, double& s2){

    if(y>=0){

        s1  = ( -y - sqrt(y*y - 4*x*z) )/(2*x);
        s2 = -(2*z)/( y + sqrt(y*y - 4*x*z) );

    } else if (y<0){

        s1 = -(2*z)/( y - sqrt(y*y - 4*x*z) );
        s2 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);

    }

    cout << "+--------------------------------------------------------+" << endl;
    cout << "Si vede che a seconda del segno di b bisogna utilizzare una forma o l'altra. Le soluzioni sono:" << endl;
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

// funzioni "norm" , "alt" e "singb" erano identiche a "quadratic.cpp" ma senza l'output



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
// gli do in input il numero di cui fare la radice, la guess, l'errore per il controllo e il riferimento al risultato
double heron(const double &S, double &x, double &err, double &y){

    // sviluppo la funzione mettendo come controllo l'errore (il cui limite è dato in input) sulla stima del metodo
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



//-------------------------------------------------- 3-Quadrature ---------------------------------------------------------//


// ------------------ quadrature1.cpp ------------------- //

// Funzione che vogliamo integrare
double fExp(double x){

    return exp(-x);

}

// Regola rettangolo
// gli do in input la funzione, gli estremi di integrazione e il numero di intervalli
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
// gli do in input la funzione, gli estremi di integrazione e il numero di intervalli
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
// gli do in input la funzione, gli estremi di integrazione e il numero di intervalli
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
// gli do in input la funzione, gli estremi di integrazione e il numero di intervalli
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

    // faccio il controllo (dimezzando l'intervallo) finché il risultato ottenuto non migliora meno della tolleranza
    while( fabs( RectangularRule(F, a, b, i) - RectangularRule(F, a, b, i/2) ) > tol ){

        i = 2*i;

    }

    // stampo a terminale il risultato della regola insieme al numero delle iterazioni
    cout << "Rectangular: " << RectangularRule(F, a, b, i) << " iter: " << i << endl;

    return 0.;

}

// Funzione di convergenza per il metodo del punto medio
// gli do in input la funzione, gli estremi di integrazione e la tolleranza
double ConvergenceMidPoint(double (*F)(double), double a, double b, double tol){

    // Controllo quanti intervalli servono
    int i = 2;

    // faccio il controllo (dimezzando l'intervallo) finché il risultato ottenuto non migliora meno della tolleranza
    while ( fabs( MidPointRule(F, a, b, i) - MidPointRule(F, a, b, i/2) ) > tol ){

        i = 2*i;

    }

    // stampo a terminale il risultato della regola insieme al numero delle iterazioni
    cout << "Mid Point: " << MidPointRule(F, a, b, i) << " iter: " << i << endl;

    return 0.;

}

// Funzione di convergenza per il metodo del trapezio
// gli do in input la funzione, gli estremi di integrazione e la tolleranza
double ConvergenceTrapezoidal(double (*F)(double), double a, double b, double tol){

    // Controllo quanti intervalli servono
    int i = 2;

    // faccio il controllo (dimezzando l'intervallo) finché il risultato ottenuto non migliora meno della tolleranza
    while ( fabs( TrapezoidalRule(F, a, b, i) - TrapezoidalRule(F, a, b, i/2) ) > tol ){

        i = 2*i;

    }

    // stampo a terminale il risultato della regola insieme al numero delle iterazioni
    cout << "Trapezoidal: " << TrapezoidalRule(F, a, b, i) << " iter: " << i << endl;

    return 0.;

}

// Funzione di convergenza per il metodo di Simpson
// gli do in input la funzione, gli estremi di integrazione e la tolleranza
double ConvergenceSimpson(double (*F)(double), double a, double b, double tol){

    // Controllo quanti intervalli servono
    int i = 2;

    // faccio il controllo (dimezzando l'intervallo) finché il risultato ottenuto non migliora meno della tolleranza
    while ( fabs( ExtSimpsonRule(F, a, b, i) - ExtSimpsonRule(F, a, b, i/2) ) > tol ){

        i = 2*i;

    }

    // stampo a terminale il risultato della regola insieme al numero delle iterazioni
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



// implemento il metodo di Gauss in cui inserisco degli if per selezionare il grado del polinomio di Legendre voluto. Puoi vedere il capitolo 5 in cui calcoliamo zeri e pesi dei polinomi di Legendre di qualsiasi ordine.
// gli do in input la funzione, gli estremi di integrazione, il numero di intervalli e il numero di punti
double Gauss(double (*F)(double), double a, double b, int N, int ng){

    // Implemento la funzione per funzionare fino a 5 punti

    // definisco gli array di pesi e zeri del polinomio di Legendre
    double w[ng], x[ng];
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

    // Calcolo l'integrale. Definisco l'ampiezza dell'intervallo e una somma ausiliaria per trovare l'integrale di Gauss
    double h = fabs(b-a)/(double)N;
    double sumj = 0.0;

    for( int i=0 ; i<N ; i++ ){

            // Devo fare il cambio di variabili che poi viene iterato sui vari intervalli successivi
            double x0 = a + i*h;
            double x1 = x0 + h;
            double s1 = (x1-x0)/2;
            double s2 = (x1+x0)/2;

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

// funzioni Trapezoidal, ExtSimpson, Gauss identiche a quelle già fatte in "quadrature.cpp" e "quadrature2.cpp".



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

// gli do in input la funzione (con due variabili), gli estremi di integrazione, sia in x che in y, il numero di intervalli e di punti gaussiani
double Gauss2D(double (*F)(double, double), double x0, double x1, double y0, double y1, int n, int ng){

    // Implemento la funzione per funzionare fino a 5 punti

    // definisco gli array di pesi e zeri del polinomio di Legendre
    double w[ng], x[ng];

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
        s1 = (b - a)/2;
        s2 = (b + a)/2;

        for( int j = 0 ; j < n ; j++ ){   // Ciclo per le y
            
            // Definisco gli estremi dell'intervallo
            c = y0 + j*dy;
            d = c + dy;
            s3 = (d - c)/2;
            s4 = (d + c)/2;

            for( int ik = 0 ; ik < ng ; ik ++ ){   // Ciclo per variare x

                for( int jk = 0 ; jk < ng ; jk ++){   // Ciclo per variare y

                    sum += w[ik] * s1 * s3 * w[jk] * F( s1*x[ik] + s2 , s3*x[jk] + s4 );

                }

            }

        }

    }

    return sum ;

}

// Funzione di convergenza per il metodo di Simpson
// gli do in input la funzione, gli estremi di integrazione e la tolleranza
void ConvergenceGauss(double (*F)(double, double), double x0, double x1, double y0, double y1, int ng, double tol){

    // Controllo quanti intervalli servono
    int i = 1;

    // faccio il controllo (incrementando l'intervallo) finché il risultato ottenuto non migliora meno della tolleranza
    while ( fabs( Gauss2D(F, x0, x1, y0, y1, i, ng) - 3.1415926535897932384626433 ) > tol ){

        i ++ ;

    }

    cout << "Gauss = " << Gauss2D(F, x0, x1, y0, y1, i, ng) << "\ncon " << i << " intervalli." << endl;

}




//-------------------------------------------------- 4-Random ---------------------------------------------------------//

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



//-------------------------------------------------- 5-Root_finders ---------------------------------------------------------//

// ------------------ froot.cpp ------------------- //

// ------------------ horner.cpp ------------------- //

// ------------------ session4.cpp ------------------- //

// ------------------ froot.cpp ------------------- //

// ------------------ legendre.cpp ------------------- //