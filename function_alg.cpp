#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// includo l'header
#include "function_alg.h"

// Dichiarazione delle funzioni in base alla lezione. Potrebbero esserci funzioni definite più volte,
// per questo motivo questo headre non è fatto per essere utilizzato, ma solo per avere una raccolta 
// di tutte le funzioni scritte nelle varie esercitazioni.


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

// ------------------ heron.cpp ------------------- //




// ------------------ quadratic.cpp ------------------- //

void dati(double& x, double& y, double& z){
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

void norm(const double& x, const double& y, const double& z, double& s1, double& s2){
    s1 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);
    s2 = ( -y - sqrt(y*y - 4*x*z) )/(2*x);
    cout << "La soluzione con (+) con il metodo standard: " << s1 << endl;
    cout << "La soluzione con (-) con il metodo standard: " << s2 << endl;
}

void alt(const double& x, const double& y, const double& z, double& s3, double& s4){
    s3 = -(2*z)/( y + sqrt(y*y - 4*x*z) );
    s4 = -(2*z)/( y - sqrt(y*y - 4*x*z) );
    cout << "La soluzione con (+) con il metodo alternativa: " << s3 << endl;
    cout << "La soluzione con (-) con il metodo alternativa: " << s4 << endl;

}

void signb(const double& x, const double& y, const double& z, double& s1, double& s2){
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

void sol(double& x1, double& x2){
    // Prendo i coefficienti
    cout << "Abbiamo l'equazione:" << endl;
    cout << "    ax^2 + bx + c = 0      " << endl;
    cout << "Dammi la soluzione 1: ";
    cin >> x1;
    cout << "Dammi la soluzione 2: ";
    cin >> x2;
}

void ordsol(double& x1, double& x2){
    // Ordino le soluzioni
    double xtemp;
    if(x2 >= x1){
        xtemp = x1;
        x1 = x2;
        x2 = xtemp;
    }
}

void norm(const double& x, const double& y, const double& z, double& s1, double& s2){
    s1 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);
    s2 = ( -y - sqrt(y*y - 4*x*z) )/(2*x);
    //cout << "La soluzione con (+) con il metodo standard: " << s1 << endl;
    //cout << "La soluzione con (-) con il metodo standard: " << s2 << endl;
}

void alt(const double& x, const double& y, const double& z, double& s3, double& s4){
    s3 = -(2*z)/( y + sqrt(y*y - 4*x*z) );
    s4 = -(2*z)/( y - sqrt(y*y - 4*x*z) );
    //cout << "La soluzione con (+) con il metodo alternativa: " << s3 << endl;
    //cout << "La soluzione con (-) con il metodo alternativa: " << s4 << endl;

}

void signb(const double& x, const double& y, const double& z, double& s1, double& s2){
    if(y>=0){
        s1  = ( -y - sqrt(y*y - 4*x*z) )/(2*x);
        s2 = -(2*z)/( y + sqrt(y*y - 4*x*z) );
    } else if (y<0){
        s1 = -(2*z)/( y - sqrt(y*y - 4*x*z) );
        s2 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);
    }
}


// ------------------ roundoff.cpp ------------------- //

// Definizione funzione sqrt
double fx1r(const double& x){
    double a = sqrt(x*x + 1.) - x;
    return a;
}

// Definizione funzione cos
double fx1c(const double& x){
    double b = 1. - cos(x);
    return b;
}

// Definizione funzione sqrt razionalizzata
// moltiplico sopra e sotto per la funzione cambiata di segno
double fx2r(const double& x){
    double c = 1./(sqrt(x*x + 1.) + x);
    return c;
}

// Definizione funzione cos razionalizzata
// moltiplico sopra e sotto per la funzione cambiata di segno
double fx2c(const double& x){
    double d = sin(x)/(1.+cos(x));
    return d;
}

// Definizione sviluppo di Taylor della sqrt
double frTaylor(const double& x){
    double c = x + 1./(2.*x) - 1./(8.*x*x*x) - x;
    return c;
}

// Definizione sviluppo di Taylor del cos
double fcTaylor(const double& x){
    double d;
    d = 1. - (1. - (x*x)/2. + (x*x*x*x)/(24.));
    return d;
}



//-------------------------------------------------- 3-Quadrature ---------------------------------------------------------//

// ------------------ quadrature1.cpp ------------------- //

// Funzione che vogliamo integrare
double func(double x){
    return exp(-x);
}

// Regola quadratura
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

double ExtSimpsonRule(double (*F)(double), double a, double b, int N){
    int w = 4; // Definizione del peso
    double h, sum;
    // Definisco l'ampiezza di ogni intervallo
    h = fabs(a-b)/(double)N;
    sum = (F(a)+F(b))*h/3;
    for (int i=1; i<N; i++ )
    {
        sum += w*F(a + i*h)*h/3;
        w = 6-w; // Calcolo il peso ogni volta
    }
    return sum;
}

double ConvergenceRectangular(double (*F)(double), double a, double b, double tol){
    // Controllo quanti intervalli servono
    int i = 2;
    while ( fabs( RectangularRule(func, a, b, i) - RectangularRule(func, a, b, i/2) ) > tol ){
        i = 2*i;
    }
    cout << "Rectangular: " << RectangularRule(func, a, b, i) << " iter: " << i << endl;
    return 0.;
}

double ConvergenceMidPoint(double (*F)(double), double a, double b, double tol){
    // Controllo quanti intervalli servono
    int i = 2;
    while ( fabs( MidPointRule(func, a, b, i) - MidPointRule(func, a, b, i/2) ) > tol ){
        i = 2*i;
    }
    cout << "Mid Point: " << MidPointRule(func, a, b, i) << " iter: " << i << endl;
    return 0.;
}

double ConvergenceTrapezoidal(double (*F)(double), double a, double b, double tol){
    // Controllo quanti intervalli servono
    int i = 2;
    while ( fabs( TrapezoidalRule(func, a, b, i) - TrapezoidalRule(func, a, b, i/2) ) > tol ){
        i = 2*i;
    }
    cout << "Trapezoidal: " << TrapezoidalRule(func, a, b, i) << " iter: " << i << endl;
    return 0.;
}

double ConvergenceSimpson(double (*F)(double), double a, double b, double tol){
    // Controllo quanti intervalli servono
    int i = 2;
    while ( fabs( ExtSimpsonRule(func, a, b, i) - ExtSimpsonRule(func, a, b, i/2) ) > tol ){
        i = 2*i;
    }
    cout << "Simpson: " << ExtSimpsonRule(func, a, b, i) << " iter: " << i << endl;
    return 0.;
}

// ------------------ quadrature2.cpp ------------------- //

double func1(double x){
    return sqrt(1.+x);
}

double func2(double x){
    return 1. - x + 2.*x*x + 0.5*x*x*x + x*x*x*x/4. - x*x*x*x*x/8.;
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
            for(int j=0; j<ng; j++){
                sumj += s1 * w[j] * F( s1*x[j] + s2 );
            }
            sum += sumj;
        } 
    
    // Restituisco il risultato
    return sum;
}


// ------------------ integral_sin.cpp ------------------- //

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
            for(int j=0; j<ng; j++){
                sumj += s1 * w[j] * F( s1*x[j] + s2 );
            }
            sum += sumj;
        } 
    
    // Restituisco il risultato
    return sum;
}


// ------------------ multid_quadrature.cpp ------------------- //

double func(double x, double y){
    return x*x*x*x*y*y + 2*x*x*y*y - x*x*y + 2;
}

double func2(double x, double y){
    if( sqrt(x*x + y*y) <= 1 ){
        return 1;
    }
    else{
        return 0;
    }
}

double Gauss2D(double (*F)(double, double), double x0, double x1, double y0, double y1, int n, int ng){

    // Implemento la funzione per funzionare fino a 5 punti
    double w[ng], x[ng]; // array di pesi e zeri del polinomio di Legendre

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
    double dx = fabs(x0-x1)/(double)n;
    double dy = fabs(y0-y1)/(double)n;

    double a , b , s1 , s2;
    double c , d , s3 , s4;

    double sum = 0.0;
    for( int i = 0 ; i < n ; i++ ){  // Ciclo per le x

        // Definisco gli estremi dell'intervallo
        a = x0 + i*dx;
        b = a + dx;
        s1 = (b - a)/2;
        s2 = (b + a)/2;

        for( int j = 0 ; j < n ; j++ ){   // Ciclo per le y
            
            // Definisco le variabili dell'intervallo
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

void ConvergenceGauss(double (*F)(double, double), double x0, double x1, double y0, double y1, int ng, double tol){
    // Controllo quanti intervalli servono
    int i = 1;
    while ( fabs( Gauss2D(func2, x0, x1, y0, y1, i, ng) - 3.1415926535897932384626433 ) > tol ){
        i ++ ;
    }
    cout << "Gauss = " << Gauss2D(func2, x0, x1, y0, y1, i, ng) << "\ncon " << i << " intervalli." << endl;
}