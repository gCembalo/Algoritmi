// Name: Gabriele Cembalo
// Date: 28/10/2025
//
// Code output:
// ************************************************************
// Trapezoidal:   n = 8192;       W = 9.919630e+00
// Gaussian   :   n = 32;       W = 9.919630e+00
// Il risultato esatto è: 9.919600e+00
// Discostamento (trapezoidal): 3.013429e-05
// Discostamento (gauss)      : 3.000711e-05
// ************************************************************

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// definisco come variabile globale i punti gaussiani
int Ng = 3;
// definisco come variabile globale il raggio 0
int r0 = 1; // r0 = GMm

double integrand(double theta);
double TrapezoidalRule(double (*)(double), double, double, int);
double Gauss(double (*)(double), double, double, int, int);

int main(){

    cout << setprecision(6);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

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
        Wn2g = Wng;

    }

    // implemento il ciclo di Trapezio (volendo puoi mettere nel while la condizione && nt < 1000 per bloccare i cicli a nt = 1000 sottointervalli )
    while( epsilont > tol  ){

        // calcolo W_{n} raddoppiando i sottointervalli
        nt *= 2;
        Wnt = TrapezoidalRule(integrand, 0.0, 4*M_PI, nt);

        // calcolo l'errore
        epsilont = fabs( Wnt - Wn2t );

        // sovrascrivo il valore di W_{n/2}
        Wn2t = Wnt;

    }

    cout << "Trapezoidal:   n = " << nt << ";       W = " << Wnt << endl;
    cout << "Gaussian   :   n = " << ng << ";       W = " << Wng << endl;

    // usando il teorema delle forze conservative:
    //      W = - Delta U       con U = -1/r
    double work = 9.9196;

    cout << "Il risultato esatto è: " << work << endl;
    cout << "Discostamento (trapezoidal): " << fabs(Wnt - work) << endl;
    cout << "Discostamento (gauss)      : " << fabs(Wng - work) << endl;


    return 0;
}

double integrand(double theta){

    // l'integrale da fare è:
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
    return force*derForce;//*r_theta*r_theta;

}

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