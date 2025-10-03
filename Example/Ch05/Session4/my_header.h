#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

int Bisection ( double (*)(double), double, double, double, double, double&, int&, int );
int FalsePos ( double (*)(double), double, double, double, double, double&, int&, int );
int Secant ( double (*)(double), double, double, double, double, double&, int&, int );
int Newton ( double (*)(double), double (*)(double), double, double, double, double&, int&, int );