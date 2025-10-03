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

double Funzione ( double );
double Derivata ( double );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 6 );
    
    int kb, kf, ks, kn, kmax = 128;
    int err_b, err_f, err_s, err_n;

    double a = -5., b = 0., c = (a+b)*0.5;
    double x_tol = 1.e-8, y_tol = 1.e-8;
    double zero_b, zero_f, zero_s, zero_n;
    
    err_b = Bisection ( Funzione, a, b, x_tol, y_tol, zero_b, kb, kmax ); 
    err_f = FalsePos ( Funzione, a, b, x_tol, y_tol, zero_f, kf, kmax ); 
    err_s = Secant ( Funzione, a, b, x_tol, y_tol, zero_s, ks, kmax );
    err_n = Newton ( Funzione, Derivata, c, x_tol, y_tol, zero_n, kn, kmax );                     

    cout << endl;
    cout << "Calcolo dello zero di x^3-3x^2+x+5" << endl;
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
    int i, n = 3;                                                        // n è il grado del polinomio

    double p = 0., a[n+1] = {5,1,-3,1};

    for ( i=n; i>=0; i-- )
    {
        p = p*x + a[i];
    }

    return p;
}

double Derivata ( double x )
{
    int i, n = 3;                                                        // n è il grado del polinomio

    double p = 0., dp = 0., a[n+1] = {5,1,-3,1};
    
    for ( i=n; i>=0; i-- )
    {
        dp = dp*x + p;
        p = p*x + a[i];
    }

    return dp;
}

int Bisection ( double (*F)(double), double a, double b, double x_tol, double y_tol, double& zero_b, int& kb, int kmax )
{
    int err;

    double xm = 0., dx = 0.;
    double fa = F(a), fb = F(b), fxm = 0.;
    
    do
    {
        xm = (b+a)*0.5;
        dx = fabs(b-a);
        fxm = F(xm);

        if ( fa==0 )
        {
            dx = 0.;
            zero_b = a;
        }
        else if ( fb==0 )
        {
            dx = 0.;
            zero_b = b;
        }
        else if ( fa*fxm<0 )
        {
            b = xm;
            fb = fxm;
            zero_b = xm;
        }
        else if ( fa*fxm>0 )
        {
            a = xm;
            fa = fxm;
            zero_b = xm;
        }
        else
        {
            dx = 0.;
            zero_b = xm;
        }

        kb++;

        if ( kb>kmax) break;
    } 
    while ( dx>x_tol );

    if ( dx<=x_tol && fabs(fxm)<=y_tol )
    {
        err = 0;                                                                             
    }
    else
    {
        err = 1;                                                                                
    }

    return err;    
}

int FalsePos ( double (*F)(double), double a, double b, double x_tol, double y_tol, double& zero_f, int& kf, int kmax )
{
    int err;

    double del = 0., xm = 0.;
    double fa = F(a), fb = F(b), fxm = 0.;

    do
    {
        xm = F(a)*(a-b)/(fb-fa)+a;                                                        
        del = fabs(xm-zero_f);
        fxm = F(xm);

        if ( fa==0 )
        {
            del = 0.;
            zero_f = a;
        }
        else if ( fb==0 )
        {
            del = 0.;
            zero_f = b;
        }
        else if ( fa*fxm<0 )
        {
            b = xm;
            fb = fxm;
            zero_f = xm;
        }
        else if ( fa*fxm>0 )
        {
            a = xm;
            fa = fxm;
            zero_f = xm;
        }
        else
        {
            del = 0.;
            zero_f = xm;
        }

        kf++;

        if ( kf>kmax) break;
    } 
    while ( del>x_tol );

    if ( del<=x_tol && fabs(fxm)<=y_tol )
    {
        err = 0;                                                                             
    }
    else
    {
        err = 1;                                                                                
    }

    return err;    
}

int Secant ( double (*F)(double), double a, double b, double x_tol, double y_tol, double& zero_s, int& ks, int kmax )
{
    int err;
    
    double fa = F(a), fb = F(b);
    double dx = 0.;
 
    do
    {
        dx = fb*(b-a)/(fb-fa);
        a = b;
        fa = fb;
        b = b-dx;
        zero_s = b;
        fb = F(b);
        ks++;

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
    
    double dx = 0.;
    double f = F(c), d = D(c);
  
    do
    {
        dx = f/d;
        c = c-dx;
        zero_n = c;
        f = F(c);
        d = D(c);
        kn++;

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