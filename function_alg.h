#ifndef function_alg
// La parte da qua in poi sara` considerato
// solo se funcion_alg non e` stato ancora definito.
// Serve per non includere piu` volte function_alg.h.

#define function_alg   // definiamo VET_H

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// Dichiarazione delle funzioni in base alla lezione. Potrebbero esserci funzioni definite più volte,
// per questo motivo questo headre non è fatto per essere utilizzato, ma solo per avere una raccolta 
// di tutte le funzioni scritte nelle varie esercitazioni.


//-------------------------------------------------- 1-Introduzione ---------------------------------------------------------//

// practice2.cpp
double sum(double, double);
double addone(double);
int Quotient(int, int, int&, int&);


//-------------------------------------------------- 2-Precision ---------------------------------------------------------//

// heron.cpp
double heron(const double &, double &, double &, double &);

// quadratic.cpp
void dati(double&, double&, double&);
void norm(const double&, const double&, const double&, double&, double&);
void alt(const double&, const double&, const double&, double&, double&);
void signb(const double&, const double&, const double&, double&, double&);

// quadratic2.cpp
void sol(double&, double&);
void ordsol(double&, double&);
void norm(const double&, const double&, const double&, double&, double&);
void alt(const double&, const double&, const double&, double&, double&);
void signb(const double&, const double&, const double&, double&, double&);

// roundoff.cpp
double fx1r(const double&);
double fx1c(const double&);
double fx2r(const double&);
double fx2c(const double&);
double frTaylor(const double&);
double fcTaylor(const double&);


//-------------------------------------------------- 3-Quadrature ---------------------------------------------------------//

// quadrature1.cpp
double func(double);
// Diversi metodi di quadratura
double RectangularRule(double (*)(double), double, double, int);
double MidPointRule(double (*)(double), double, double, int);
double TrapezoidalRule(double (*)(double), double, double, int);
double ExtSimpsonRule(double (*)(double), double, double, int);
// Funzioni per la convergenza
double ConvergenceRectangular(double (*)(double), double, double, double);
double ConvergenceMidPoint(double (*)(double), double, double, double);
double ConvergenceTrapezoidal(double (*)(double), double, double, double);
double ConvergenceSimpson(double (*)(double), double, double, double);

// quadrature2.cpp
double func1(double);
double func2(double);
// Diversi metodi
double Simpson(double (*)(double), double, double, int);
double Gauss(double (*)(double), double, double, int, int);

// integral_sin.cpp
double func(double); // funzione
// Diversi metodi
double Trapezoidal(double (*)(double), double, double,int);
double Simpson(double (*)(double), double, double, int);
double Gauss(double (*)(double), double, double, int, int);

// multid_quadrature.cpp
double func(double, double);
double func2(double, double);
double Gauss2D(double (*)(double, double), double, double, double, double, int, int);
void ConvergenceGauss(double (*F)(double, double), double , double , double , double , int , double);


//-------------------------------------------------- 4- ---------------------------------------------------------//






#endif // fine del if della prima riga