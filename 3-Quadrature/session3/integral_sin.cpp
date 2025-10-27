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

int main(){
    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    fdata << setprecision(7);  //set the significant digits
    fdata << setiosflags ( ios::scientific );  //set the scientific notation

    double x0 = 0.; // L'estremo di integrazione
    double x = 0.8; // Il nostro punto
    // I diversi intervalli li seleziono dando diversi valori di N, così che si divida (x-x0) in N parti e si ottengano effettivamente 0.8, 0.4, 0.2, 0.1

    cout << "\n+-----------------------------------+" << endl;
    cout << "\nValore esatto Si(0.8) = 0.77209 57855\n" << endl;
    cout << "Valutazione con i diversi metodi:\n" << endl;
    cout << "Gauss (n = 1, ng = 3) = " << Gauss(func,x0,x,1,3) << endl;
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
        gauss += Gauss(func,x0,x,1,3);
        fdata << x << "    " << gauss << endl;
        x0 = x;
    }
    fdata.close();  // Chiudo il file

    cout << "\n+-----------------------------------+\n" << endl;
    double asint = 0.0;
    double x1 = 0.0;
    for( double i = 0.0 ; i<=1000 ; i += 0.1 ){
        // Scrivo il file
        asint += Gauss(func,x1,i,1,3);
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