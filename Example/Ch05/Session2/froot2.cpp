#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

int Secant ( double (*)(double), double, double, double, double, double&, int&, int );
int Newton ( double (*)(double), double (*)(double), double, double, double, double&, int&, int );

double Funzione ( double );
double Derivata ( double );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 6 );
    
    int ks, kn, kmax = 128;
    int err_s, err_n;
    
    double a = -1., b = 1., c = (a+b)*0.5;
    double x_tol = 1.e-7, y_tol = 1.e-7;
    double zero_s, zero_n;

    err_s = Secant ( Funzione, a, b, x_tol, y_tol, zero_s, ks, kmax );
    err_n = Newton ( Funzione, Derivata, c, x_tol, y_tol, zero_n, kn, kmax );                        

    cout << endl;
    cout << "Calcolo dello zero di e^(-x)-x" << endl;
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << "secant: \t\t err = " << err_s << "\t k = " << ks << "\t\t x0 = " << zero_s << "\t F(x0) = " << Funzione(zero_s) << endl;
    cout << "newton: \t\t err = " << err_n << "\t k = " << kn << "\t\t x0 = " << zero_n << "\t F(x0) = " << Funzione(zero_n) << endl;
    cout << endl;
 
    return 0;
}

double Funzione ( double x )
{
    return exp(-x)-x;
}

double Derivata ( double x )
{
    return -exp(-x)-1;
}

int Secant ( double (*F)(double), double a, double b, double x_tol, double y_tol, double& zero_s, int& ks, int kmax )
{
    int err;
    
    double fa = F(a), fb = F(b);
    double dx = 0.;
 
    do
    {
        dx = fb*(b-a)/(fb-fa);
        ks++;
        
        // cout << "Secant:  k = " << ks << "\t a = " << a << "\t b = " << b << "\t dx = " << dx << endl;    // debugging

        a = b;
        fa = fb;
        b = b-dx;
        zero_s = b;
        fb = F(b);

        if ( ks>kmax) break;
    }
    while ( fabs(dx)>x_tol );

    if ( dx<=x_tol && fabs(fb)<=y_tol )
    {
        err = 0;                                                                             
    }
    else
    {
        err = 1;                                                                                
    }

    return err;  
}

int Newton ( double (*F)(double), double (*D)(double), double c, double x_tol, double y_tol, double& zero_n, int& kn, int kmax )
{
    int err;
    
    double f = F(c), d = D(c);
    double dx = 0.;
  
    do
    {
        dx = f/d;
        kn++;
        
        // cout << "Newton:  k = " << kn << "\t c = " << c << "\t dx = " << dx << endl;    // debugging

        c = c-dx;
        zero_n = c;
        f = F(c);
        d = D(c);

        if ( kn>kmax) break;
    }
    while ( fabs(dx)>x_tol );

    if ( dx<=x_tol && fabs(f)<=y_tol )
    {
        err = 0;                                                                             
    }
    else
    {
        err = 1;                                                                                
    }

    return err;
}