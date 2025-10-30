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
    double xi = 1.0 , xi1 , x1i , xi2 , x2i; // x_i è x_{i} ; xi1 è x_{i+1} ; x1i è x_{i-1} ; x2i è x_{i-2} ; xi2 è x_{i+2}

    // definisco le variabili in cui metto i risultati delle derivate
    double ValderFD , ValderBD, ValderCD, Valder4th;

    // definisco le variabili degli errori
    double errDerFD , errDerBD , errDerCD, errDer4th;

    // definisco la variabile per la derivata esatta
    double DerEx = sinDerEx(xi);

    for( int i = 0 ; i <= 10 ; i++ ){

        // modifico i punti delle x
        xi1 = xi + h;
        x1i = xi - h;
        x2i = xi - 2*h;
        xi2 = xi + 2*h;

        // calcolo le derivate con i diversi metodi
        ValderFD = derFD(sinDer, xi, xi1, h);
        ValderBD = derBD(sinDer, xi, x1i, h);
        ValderCD = derCD(sinDer, x1i, xi1, h);
        Valder4th = der4th(sinDer, x2i, x1i, xi1, xi2, h);

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

    return 0;

}


double sinDer(double x){

    return sin(x);

}

double sinDerEx(double x){

    return cos(x);

}

// x_i è x_{i} ; xi1 è x_{i+1} ; x1i è x_{i-1}
double derFD(double (*F)(double), double xi, double xi1, double h){

    // definisco 
    double fiPrime = ( F(xi1) - F(xi) )/h;

    return fiPrime;

}

// x_i è x_{i} ; xi1 è x_{i+1} ; x1i è x_{i-1}
double derBD(double (*F)(double), double xi, double x1i, double h){

    // definisco 
    double fiPrime = ( F(xi) - F(x1i) )/h;

    return fiPrime;

}

// x_i è x_{i} ; xi1 è x_{i+1} ; x1i è x_{i-1}
double derCD(double (*F)(double), double x1i, double xi1, double h){

    // definisco 
    double fiPrime = ( F(xi1) - F(x1i) )/h;

    return fiPrime;

}

// x_i è x_{i} ; xi1 è x_{i+1} ; x1i è x_{i-1} ; x2i è x_{i-2} ; xi2 è x_{i+2}
double der4th(double (*F)(double), double x2i, double x1i, double xi1, double xi2, double h){

    // definisco 
    double fiPrime = ( F(x2i) - 8.0*F(x1i) + F(xi1) - F(xi2) )/( 12.0*h );

    return fiPrime;

}