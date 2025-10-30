// compute the numerical derivative f(x)=sin(x) in x=1 using FD, BD and CD (or higher) using different increments h=0.5,0.25,0.125, …
//
// Plot the error as a function of h using a log-log scaling.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double sinDer(double);
double sinDerEx(double);
double derFD(double (*)(double), double, double, double);
double derBD(double (*)(double), double, double, double);
double derCD(double (*)(double), double, double, double);
double der4th(double (*)(double), double, double, double, double, double);

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata1, fdata2, fdata3, fdata4;
    fdata1.open("FD.dat");
    fdata2.open("BD.dat");
    fdata3.open("CD.dat");
    fdata4.open("4thorder.dat");

    //set the significant digits and the scientific notation
    fdata1 << setprecision(7);
    fdata1 << setiosflags ( ios::scientific );
    fdata2 << setprecision(7);
    fdata2 << setiosflags ( ios::scientific );
    fdata3 << setprecision(7);
    fdata3 << setiosflags ( ios::scientific );
    fdata4 << setprecision(7);
    fdata4 << setiosflags ( ios::scientific );

    // definisco la spaziatura
    double h = 0.5; // poi lo dimezzerò ad ogni ciclo
    // definisco le mie x_i , x_{i-1} e x_{i+1}
    double xi = 1.0 , xm , xp , xmm , xpp; // x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}

    // definisco le variabili in cui metto i risultati delle derivate
    double ValderFD , ValderBD, ValderCD, Valder4th;

    // definisco le variabili degli errori
    double errDerFD , errDerBD , errDerCD, errDer4th;

    // definisco la variabile per la derivata esatta
    double DerEx = sinDerEx(xi);

    for( int i = 0 ; i <= 10 ; i++ ){

        // modifico i punti delle x
        xp = xi + h;
        xm = xi - h;
        xmm = xi - 2*h;
        xpp = xi + 2*h;

        // calcolo le derivate con i diversi metodi
        ValderFD = derFD(sinDer, xi, xp, h);
        ValderBD = derBD(sinDer, xi, xm, h);
        ValderCD = derCD(sinDer, xm, xp, h);
        Valder4th = der4th(sinDer, xmm, xm, xp, xpp, h);

        // calcolo il valore degli errori
        errDerFD = fabs( ValderFD - DerEx );
        errDerBD = fabs( ValderBD - DerEx );
        errDerCD = fabs( ValderCD - DerEx );
        errDer4th = fabs( Valder4th - DerEx );

        // inserico i valori nel file .dat
        fdata1 << 1.0/h << "     " << errDerFD << endl;
        fdata2 << 1.0/h << "     " << errDerBD << endl;
        fdata3 << 1.0/h << "     " << errDerCD << endl;
        fdata4 << 1.0/h << "     " << errDer4th << endl;

        // modifico il valore di h
        h /= 2.0;

    }

    fdata1.close();
    fdata2.close();
    fdata3.close();
    fdata4.close();

    return 0;

}


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
double der4th(double (*F)(double), double xmm, double xm, double xp, double xpp, double h){

    // definisco 
    double fiPrime = ( F(xmm) - 8.0*F(xm) + 8.0*F(xp) - F(xpp) )/( 12.0*h );

    return fiPrime;

}