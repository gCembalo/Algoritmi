#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

const int g_leg = 3;                      // grado del polinomio di legendre

void Bracket ( double (*F)(double), double, double, double*, double*, int, int& );

int Bisection ( double (*F)(double), double, double, double, double, double&, int&, int );
int FalsePos ( double (*F)(double), double, double, double, double, double&, int&, int );
int Secant ( double (*F)(double), double, double, double, double, double&, int&, int );
int Newton ( double (*F)(double), double (*D)(double), double, double, double, double&, int&, int );

double Funzione ( double );
double Derivata ( double );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 12 );
    
    int err_b, err_f, err_s, err_n;
    int kb, kf, ks, kn, kmax = 99;
    int i, n = 11, nr = 32;

    double a = -1., b = 1.;
    double x_tol = 1.e-7, y_tol = 1.e-7;
    double zero_b, zero_f, zero_s, zero_n;
    double xl[nr], xr[nr];
    double wi_b, wi_f, wi_s, wi_n;

    cout << endl;
    cout << "Calcolo degli zeri del polinomio di Legendre di grado " << g_leg << endl;

    Bracket ( Funzione, a, b, xl, xr, n, nr );

    if ( nr==0 )
    {
        cout << "No roots" << endl;
    }
    
    for ( i=0; i<nr; i++ )
    {
        err_b = Bisection ( Funzione, xl[i], xr[i], x_tol, y_tol, zero_b, kb, kmax );
        err_f = FalsePos ( Funzione, xl[i], xr[i], x_tol, y_tol, zero_f, kf, kmax );
        err_s = Secant ( Funzione, xl[i], xr[i], x_tol, y_tol, zero_s, ks, kmax );
        err_n = Newton ( Funzione, Derivata, (xl[i]+xr[i])*0.5, x_tol, y_tol, zero_n, kn, kmax );

        wi_b = 2/((1-zero_b*zero_b)*Derivata(zero_b)*Derivata(zero_b));
        wi_f = 2/((1-zero_f*zero_f)*Derivata(zero_f)*Derivata(zero_f));
        wi_s = 2/((1-zero_s*zero_s)*Derivata(zero_s)*Derivata(zero_s));
        wi_n = 2/((1-zero_n*zero_n)*Derivata(zero_n)*Derivata(zero_n));
        
        cout << "root # " << i+1 << " --------------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << "bisection: \t\t err = " << err_b << "\t k = " << kb << "\t\t x0 = " << zero_b << "\t F(x0) = " << Funzione(zero_b) << "\t w = " << wi_b << endl;
        cout << "false position: \t err = " << err_f << "\t k = " << kf << "\t\t x0 = " << zero_f << "\t F(x0) = " << Funzione(zero_f) << "\t w = " << wi_f << endl;
        cout << "secant: \t\t err = " << err_s << "\t k = " << ks << "\t\t x0 = " << zero_s << "\t F(x0) = " << Funzione(zero_s) << "\t w = " << wi_s << endl;
        cout << "newton: \t\t err = " << err_n << "\t k = " << kn << "\t\t x0 = " << zero_n << "\t F(x0) = " << Funzione(zero_n) << "\t w = " << wi_n << endl;

        zero_b = 0., kb = 0;
        zero_f = 0., kf = 0;
        zero_s = 0., ks = 0;
        zero_n = 0., kn = 0;
    }
    
    cout << endl;
    
    return 0;
}

double Funzione ( double x )
{
    int i;

    double pn = 0., p0 = 1., p1 = x;
    
    if ( g_leg==0 ) 
    {
        return p0;
    }
    else if ( g_leg==1 )
    {
        return p1;
    }
    else 
    {
        for ( i=2; i<=g_leg; i++ )
        {
            pn = ((2*i-1)*x*p1-(i-1)*p0)/(double)i;
            p0 = p1;
            p1 = pn;   
        }
        
        return pn;
    }
}

double Derivata ( double x )
{
    int i;

    double d_pn = 0., d_p0 = 0., pn = 0., p0 = 1., p1 = x;

    if ( g_leg==0 ) 
    {
        return d_p0;
    }
    else if ( g_leg==1 )
    {
        return p0;
    }
    else 
    {
        for ( i=2; i<=g_leg; i++ )
        {
            pn = ((2*i-1)*x*p1-(i-1)*p0)/(double)i;
            d_pn = (i*(x*pn-p1))/(x*x-1);
            p0 = p1;
            p1 = pn;   
        }
        
        return d_pn;
    }
}

void Bracket ( double (*F)(double), double a, double b, double* xl, double* xr, int n, int& nr )
{
    int i, j = 0;

    double fa = F(a), fb = 0.;
    double h, xa = 0., xb = 0.;

    h = (b-a)/(double)n;
    
    for ( i=0; i<n; i++ )
    {
        xa = a+i*h;
        xb = a+(i+1)*h;
        fb = F(xb);
        
        if ( fa*fb<0 )
        {
            xl[j] = xa;
            xr[j] = xb;
            j++;
        }
        else if ( fa==0 || ( fb==0 && xb==b ))
        {
            xl[j] = xa;
            xr[j] = xb;
            j++;
        }

        fa = fb;
    }

    nr = j;
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

    double f = F(c), d = D(c);
    double dx = 0.;
    
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