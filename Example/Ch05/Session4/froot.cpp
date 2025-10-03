#include "my_header.h"

double Funzione ( double );
double Derivata ( double );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 6 );
    
    int err_b, err_f, err_s, err_n;
    int kb, kf, ks, kn, kmax = 99;

    double a = -1., b = 1., c = (a+b)*0.5;
    double x_tol = 1.e-7, y_tol = 1.e-7;
    double zero_b, zero_f, zero_s, zero_n;
    
    err_b = Bisection ( Funzione, a, b, x_tol, y_tol, zero_b, kb, kmax );
    err_f = FalsePos ( Funzione, a, b, x_tol, y_tol, zero_f, kf, kmax );
    err_s = Secant ( Funzione, a, b, x_tol, y_tol, zero_s, ks, kmax );
    err_n = Newton ( Funzione, Derivata, c, x_tol, y_tol, zero_n, kn, kmax );                       

    cout << endl;
    cout << "Calcolo dello zero di e^(-x)-x" << endl;
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << "bisection: \t\t err = " << err_b << "\t k = " << kb << "\t\t x0 = " << zero_b << "\t F(x0) = " << Funzione(zero_b) << endl;
    cout << "false position: \t err = " << err_f << "\t k = " << kf << "\t\t x0 = " << zero_f << "\t F(x0) = " << Funzione(zero_f) << endl;
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