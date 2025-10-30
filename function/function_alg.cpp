//
// Vedi le slide ch05 per vedere come organizzare il tuo codice.
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// includo l'header
#include "function_alg.h"

int g_LegendreN; // variabile che serve nel capitolo 5 per i polinomi di Legendre

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

// --------------------------------------- //
// definisco come variabile globale i punti gaussiani
//int Ng = 3;
// definisco come variabile globale il raggio 0
//int r0 = 1; // r0 = GMm
/*int main(){
    // definisco la tolleranza
    double tol = 1.e-6;
    // definisco i sottointervalli
    int ng = 2 , nt = 2;
    // calcolo gli integrali con i due metodi mettendo la condizione sull'errore e sul numero di sotto intervalli
    double epsilong = 10.0; // errore; lo inizializzo così per fare partire il ciclo
    double epsilont = 10.0; // errore; lo inizializzo così per fare partire il ciclo
    double Wng, Wn2g; // sono i lavori W_n e W_{n/2} che calcolo con gauss
    double Wnt, Wn2t; // sono i lavori W_n e W_{n/2} che calcolo con il trapezio
    // inizializzo i due integrali con 2 sotto intervalli
    Wn2g = Gauss(integrand, 0, 4*M_PI, ng, Ng);
    Wn2t = TrapezoidalRule(integrand, 0, 4*M_PI, ng);
    // implemento il ciclo di Gauss
    while( epsilong > tol && ng < 1000 ){
        // calcolo W_{n} raddoppiando i sottointervalli
        ng *= 2;
        Wng = Gauss(integrand, 0.0, 4*M_PI, ng, Ng);
        // calcolo l'errore
        epsilong = fabs( Wng - Wn2g );
        // sovrascrivo il valore di W_{n/2}
        Wn2g = Wng; }
    // implemento il ciclo di Trapezio (volendo puoi mettere nel while la condizione && nt < 1000 per bloccare i cicli a nt = 1000 sottointervalli )
    while( epsilont > tol  ){
        // calcolo W_{n} raddoppiando i sottointervalli
        nt *= 2;
        Wnt = TrapezoidalRule(integrand, 0.0, 4*M_PI, nt);
        // calcolo l'errore
        epsilont = fabs( Wnt - Wn2t );
        // sovrascrivo il valore di W_{n/2}
        Wn2t = Wnt; }
    cout << "Trapezoidal:   n = " << nt << ";       W = " << Wnt << endl;
    cout << "Gaussian   :   n = " << ng << ";       W = " << Wng << endl;
    // usando il teorema delle forze conservative:
    //      W = - Delta U       con U = -1/r
    double work = 9.9196;
    cout << "Il risultato esatto è: " << work << endl;
    cout << "Discostamento (trapezoidal): " << fabs(Wnt - work) << endl;
    cout << "Discostamento (gauss)      : " << fabs(Wng - work) << endl;
    return 0; }
double integrand(double theta){ // l'integrale da fare è:
    //      \int \vec{F}_{grav} * ds
    // con ds = \hat{r} dr , dunque l'integranda è:
    //      F_{grav} dr
    // ma non possiamo integrare in dr e dobbiamo cambiare in dtheta:
    //      dr = (dr/dtheta) dtheta
    // in questo modo l'integrale da fare diventa:
    //      \int_0^{4\pi} F_{grav} * (dr/dtheta) * dtheta
    // per cui dobbiamo anche calcolare lo Jacobiano (dr/dtheta), che è:
    //      r0 * [ e^(-b/theta)/b - e^(-b/theta)*(1+b/theta)/b ]
    // riscritto come:
    //      -( r0 * theta )/( b*b )*e^( -theta/b )
    double b = M_PI;
    double r_theta = (double)r0 * ( 1.0 + theta/b ) * exp( -theta/b ); // orbita
    double force = - (double)r0/( r_theta * r_theta ); // forza
    double derForce = - exp( -theta/b )*( (double)r0 * theta )/( b*b );
    // restituisco la funzione da integrare
    return force*derForce;//*r_theta*r_theta; }*/
// --------------------------------------- //

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

// --------------------------------------- //
// potenziali e quanto variabili globali
// definisci tolleranza, gli zeri e le iterazioni dei metodi. Definisci anche due variabili che ti fanno da estremi
// stampo a video gli zeri trovati con i due metodi
/*  cout << "Bisection, results:" << endl;
    for( s = 1 ; s <= 3 ; s++ ){ // intervalli che definisco per ridurre leggermente 0<E<80; da una prima compilazione vedo che il secondo 0 è a 28, per cui il primo sarà sicuramente minore.
        double a = ( (s - 1) * 20.0 + 1.0 );
        double b = s * 25.0;
        // richiamo la funzione
        bisection(Energy, a, b, tol, x1, l1);
        // stampo l'output
        cout << "s = " << s << "; Root: " << x1 << "; ntry = " << l1 << endl;}
    cout << "Secant, results:" << endl;
    for( s = 1 ; s <= 3 ; s++ ){ // intervalli che definisco per ridurre leggermente 0<E<80; da una prima compilazione vedo che il secondo 0 è a 28, per cui il primo sarà sicuramente minore.
        double a = ( (s - 1) * 20.0 + 1.0 );
        double b = s * 25.0;
        // richiamo la funzione
        secant_method(Energy, a, b, tol, x2, l2);
        // stampo l'output
        cout << "s = " << s << "; Root: " << x2 << "; ntry = " << l2 << endl;}
// definisco la funzione dei livelli energetici
double Energy(double E){ return sqrt(E) - (double)s*M_PI + asin( sqrt( E/V1 ) ) + asin( sqrt( E/V2 ) ); }
// --------------------------------------- // */


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

// --------------------------------------- //
/*// definisco come variabile globale
int ng = 5; // punti gaussiani
double tlo = 0.0; // tempo iniziale
int main(){ // definisco le costanti
    double tol = 1.e-6; // tolleranza sia x che y
    // definisco le iterazioni
    int l;
    // definisco il punto di inversione
    double inv;
    // richiamo Newton per cercare lo zero
    // int newton_method(double (*F)(double), double (*derF)(double), double a, double b, double xtol, double ytol, double &zero, int &l)
    newton_method(vel, acc, tlo, 10, tol, tol, inv, l);
    // stampo a terminale lo zero
    cout << "Inversion time (zero) = " << inv << endl;
    return 0;}
double vel(double t){ static int nfv = 0; // Cumulative number of function evaluations
    // definisco le costanti che mi servono
    double v0 = -1.0; // velocità iniziale
    double thi = t;
    // questo è il numero di intervalli
    int nint = ceil( 2.0*fabs( tlo - thi ) );
    double integral = Gauss(acc, tlo, thi, nint, ng);
    nfv += ng*nint;
    // stampo a terminale
    cout << "tlo = " << tlo << ";    thi = " << thi << ";    nint = " << nint << ";     Func eval = " << nfv << endl;
    // ritorno il valore di v(t)
    return v0 + integral;}
double acc(double t){ return ( 1.0 - exp(-t) )/( sqrt( 1.0 + t*t*t*t ) ); } */
// --------------------------------------- //


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


// ------------ root finders function ------------- //

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

// metodo della false position
// gli do in input la funzione, gli estremi a e b, la tolleranza su x, uno zero per riferimento e il numero di iterazioni
int false_position(double (*F)(double), double a, double b, double tol, double &zero, int &l){

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
        cout << "(False position) # = " << n << "    (l'estremo " << a << " è già lo zero)"<< endl;

        l = n;
        return 0;

    }
    else if( fb == 0.0 ){

        zero = b;

        // creo l'output delle iterazioni
        cout << "(False position) # = " << n << "(l'estremo " << b << " è già lo zero)"<< endl;
        
        l = n;
        return 0;

    }
    else{
        // metto nel ciclo la condizione sia sulla tolleranza che sul numero di cicli (inserita dopo)
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
            //cout << "n = " << n << ";   [a,b] = [" << a << ", " << b << "];    xm = " << x << ";   Deltax = " << fabs(a-b) << ";   f(xm) = " << F(x) << endl;

        }

        // creo l'output delle iterazioni
        //cout << "(False position) # = " << n << endl;

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

// bracketing function
// gli do in input la funzione, gli estremi a e b, il numero di intervalli, un array per gli estremi (sinistro e destro) degli intervalli, e un riferimento al numero di zeri
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

    return exp( 1/( x + 0.5 ) ) - ( 3 + 2*x )/( 1 + x );

}

double derfuncSes4(double x){

    return - 1 / ((x + 0.5)*(x+0.5)) * exp( 1 / (x + 0.5)) + 1 / ((1+x) * (1+x));

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



//-------------------------------------------------- 6- ---------------------------------------------------------//