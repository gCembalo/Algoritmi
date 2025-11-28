//
// Vedi le slide ch05 per vedere come organizzare il tuo codice.
//

#ifndef function_alg
// La parte da qua in poi sara` considerato
// solo se funcion_alg non e` stato ancora definito.
// Serve per non includere piu` volte function_alg.h.

#define function_alg   // definiamo VET_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

// Dichiarazione delle funzioni in base alla lezione. Potrebbero esserci 
// funzioni definite più volte, quindi è sempre meglio controllare, non solo 
// il funzionamento, ma anche la sintassi con cui le si chiama.


//-------------------------- 1-Introduzione ---------------------------------//

// practice2.cpp
double sum(double, double);
double addone(double);
int Quotient(int, int, int&, int&);


//--------------------------- 2-Precision ----------------------------------//

// quadratic.cpp
void datiQuad(double&, double&, double&);
void normQuad(const double&, const double&, const double&, double&, double&);
void altQuad(const double&, const double&, const double&, double&, double&);
void signbQuad(const double&, const double&, const double&, double&, double&);

// quadratic2.cpp
void solQuad(double&, double&);
void ordsolQuad(double&, double&);
// alcune funzioni erano identiche a "quadratic.cpp"

// roundoff.cpp
double fx1r(const double&);
double fx1c(const double&);
double fx2r(const double&);
double fx2c(const double&);
double frTaylor(const double&);
double fcTaylor(const double&);

// heron.cpp
double heron(const double &, double &, double &, double &);


//-------------------------- 3-Quadrature ---------------------------------//

// quadrature1.cpp
double fExp(double);
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
double fSqrt(double);
double func2(double);
// Diversi metodi
// il metodo di Simpson è identico a quello in "quadrature.cpp"
double Gauss(double (*)(double), double, double, int, int);

// integral_sin.cpp
double funcSinx(double); // funzione
// Diversi metodi implementati in altre esercitazioni

// multid_quadrature.cpp
double func1(double, double);
double funcCirc(double, double);
double Gauss2D(double (*)(double, double), double, double, double,
               double, int, int);
void ConvergenceGauss(double (*F)(double, double), double , double , 
               double , double , int , double);


//------------------------ 4-Random_numbers -------------------------------//

// gauss_distrib.cpp
double gaussiana(const double&, const double&);



//------------------------- 5-Root_finders --------------------------------//

// funzioni uguali per tutti i file del capitolo.
// metodi di root finding
int bisection(double (*F)(double), double, double, double, double&, int &);
int false_position(double (*F)(double), double, double, double, double&, int &);
int secant_method(double (*F)(double), double, double, double, double &, int &);
int newton_method(double (*F)(double), double (*derF)(double), double, double,
                  double, double, double &, int &);

// metodo di bracketing
int Bracket(double (*F)(double), double a, double b, double n, double *xL, 
            double *xR, int &nroots);


// froot.cpp
double funcFroot(double &);
double derfuncFroot(double &);

// horner.cpp
double pol(double);
double derpol(double);

// session4.cpp
double funcSes4(double);
double derfuncSes4(double);

// froot.cpp
double funcSin(double);
double derfuncSin(double);

// legendre.cpp
double polLegendre(double);
double derpolLegendre(double);
double wi(double);


//-------------------------- 6-Derivative ---------------------------------//

// derivative.cpp
double sinDer(double);
double sinDerEx(double);
double derFD(double (*)(double), double, double, double);
double derBD(double (*)(double), double, double, double);
double derCD(double (*)(double), double, double, double);
double der4th(double (*)(double), double, double, double, double, double);

// trajectory.cpp

// le funzioni derivate (derFD , derBD , derCD , der4th) sono uguali a "derivative.cpp"
double position(double);
double velocity(double);
double SecondDerivative(double (*)(double), double, double, double, double);


//------------------------------ 7-ODE -------------------------------------//

// funzioni generiche
void EulerStep(double, double *, void (*)(double, double *, double *), 
               double, int);
void RK2StepMid(double, double *, void (*)(double, double *, double *), 
                double, int);
void RK2StepHeun(double, double *, void (*)(double, double *, double *), 
                 double, int);
void RK4Step(double, double *, void (*)(double, double *, double *), 
             double, int);
void PositionVerletStep(double *, double *, int, double, void (*)(double *,
                        double *, int));
void PositionVerletStep2(double, double *, void (*)(double,
                         double *, double *), double, int);
void VelocityVerletStep(double *, double *, int, double, void (*)(double *,
                        double *, int));

// ode1.cpp
void RHSFuncOde1(double, double *,double *);
double ode1Sol(double);

// ode2.cpp
double ode2Sol(double);
void RHSFuncOde2(double, double *,double *);

// kepler.cpp
void RHSFuncOde3(double, double *,double *);

// harmonic.cpp
void RHSFuncOde4(double, double *,double *);
void acceleration(double *, double *, int);



//------------------------------ 8-ODE -------------------------------------//

// poisson.cpp
void RHSFuncPoisson(double, double *, double *);
double ResidualPoisson(double);

// wave.cpp
void RHSFuncWave(double, double *, double *);
double ResidualWave(double);

// qho.cpp
void RHSFuncQho(double, double *, double *);
double ResidualQho(double);



//--------------------------- 9-MultiArray ----------------------------------//

// matrix.cpp
// gauss_elim.cpp
// gauss_elim2.cpp
// tridiag.cpp
// bvp.cpp



//--------------------------- 10-PDE ----------------------------------//

// elliptic.cpp
// example.cpp
// example.cpp


#endif // fine del if della prima riga