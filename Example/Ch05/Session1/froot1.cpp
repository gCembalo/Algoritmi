#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

int Bisection ( double (*)(double), double , double , double, double, double&, int&, int );
int FalsePos ( double (*)(double), double , double , double, double, double&, int&, int );

double Funzione ( double );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 6 );
    
    int kb, kf, kmax = 128;
    int err_b, err_f;
    
    double a = -1., b = 1.;
    double x_tol = 1.e-7, y_tol = 1.e-7;
    double zero_b, zero_f;

    err_b = Bisection ( Funzione, a, b, x_tol, y_tol, zero_b, kb, kmax );                                      
    err_f = FalsePos ( Funzione, a, b, x_tol, y_tol, zero_f, kf, kmax );                                        

    cout << endl;
    cout << "Calcolo dello zero di e^(-x)-x" << endl;
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << "bisection: \t\t err = " << err_b << "\t k = " << kb << "\t\t x0 = " << zero_b << "\t F(x0) = " << Funzione(zero_b) << endl;
    cout << "false position: \t err = " << err_f << "\t k = " << kf << "\t\t x0 = " << zero_f << "\t F(x0) = " << Funzione(zero_f) << endl;
    cout << endl;
    
    return 0;
}

double Funzione ( double x )
{
    return exp(-x)-x;
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
        kb++;

        // cout << "Bisection: k = " << kb << "\t [a,b] = [" << a << "," << b << "] \t xm = " << xm << "\t dx = " << dx << "\t fm = " << fxm << endl;    // debugging

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
        else                                                      // cioè fm==0
        {
            dx = 0.;
            zero_b = xm;
        }

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

    double del = 0., xm = 0., fa = F(a), fb = F(b), fxm = 0.;

    do
    {
        xm = F(a)*(a-b)/(fb-fa)+a;                               // retta per due punti                                                    
        del = fabs(xm-zero_f);
        fxm = F(xm);
        kf++;

        // cout << "FalsePos:  k = " << kf << "\t [a,b] = [" << a << "," << b << "]" << "\t xm = " << xm << "\t fm = " << fxm << "\t |del| = " << del << endl;    // debugging

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
        else                                                     // cioè fm==0
        {
            del = 0.;
            zero_f = xm;
        }

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