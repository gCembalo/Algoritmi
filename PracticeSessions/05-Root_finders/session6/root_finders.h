#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

// metodi di root finding
int bisection(double (*F)(double), double, double, double, double&, int &);
int false_position(double (*F)(double), double, double, double, double&, int &);
int secant_method(double (*F)(double), double, double, double, double &, int &);
int newton_method(double (*F)(double), double (*derF)(double), double, double, double, double, double &, int &);

// metodo di bracketing
int Bracket(double (*F)(double), double a, double b, double n, double *xL, double *xR, int &nroots);