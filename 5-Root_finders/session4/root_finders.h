#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

// metodi di root finding
int bisection(double (*F)(double), double, double, double, double&);
int false_position(double (*F)(double), double, double, double, double&);
int secant_method(double (*F)(double), double, double, double, double &);
int newton_method(double (*F)(double), double (*derF)(double), double, double, double &);