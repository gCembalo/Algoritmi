#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

static double g_E;

void dYdx ( double, double*, double* );
void RK4 ( double, double*, void (*)(double,double*,double*), double, int );
void Bracket ( double (*)(double), double, double, double*, double*, int, int& );

int Bisection ( double (*)(double), double, double, double, double, double&, int&, int );

double Residual ( double );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 8 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 8 );

    int i, neq = 2;
    int nst = 800, n = 1000;
    int nint = 20, nr = 10;
    int err, k, kmax = 99;

    double Y[neq], x, dx;
    double xl = -10., xr = 10.;
    double E, res, zero;
    double Ei = 0., Ef = 5.;
    double El[nr], Er[nr];
    double xtol = 1.e-10, ytol = 1.e-8;
    double psi[nst], xx[nst];

    fdata.open ( "qho.dat" );

    g_E = 0.5;
    dx = (xr-xl)/(double)nst;
    x = xl;

    Y[0] = exp(-x*x*0.5);
    Y[1] = -x*Y[0];
    
    for ( i=1; i<nst; i++ )
    {
        RK4 ( x, Y, dYdx, dx, neq );

        x += dx;

        fdata << x << " " << Y[0] << endl;
    }

    fdata << endl << endl;

    dx *= -1.;
    x = xr;

    Y[0] = exp(-x*x*0.5);
    Y[1] = -x*Y[0];
    
    for ( i=1; i<nst; i++ )
    {
        RK4 ( x, Y, dYdx, dx, neq );

        x += dx;

        fdata << x << " " << Y[0] << endl;
    }

    fdata << endl << endl;

    for ( i=0; i<n; i++ )
    {
        E = 5*i/(double)(n-1)+1.e-3;
        
        res = Residual(E);
        
        fdata << E << " " << res << endl;
    }

    fdata.close ();
    
    cout << endl;
    
    Bracket ( Residual, Ei, Ef, El, Er, nint, nr );

    for ( i=0; i<nr; i++ )
    {
        err = Bisection ( Residual, El[i], Er[i], xtol, ytol, zero, k, kmax );
 
        if ( err==0 )
        {
            cout << "E" << i+1 << " = " << zero << "\t k = " << k << endl;
        }
        else 
        {
            cout << "Nessuna soluzione" << endl;
        }
        
        k = 0;
    }

    cout << endl;
    
    return 0;
}

void dYdx ( double x, double* Y, double* R )
{   
    R[0] = Y[1];
    R[1] = x*x*Y[0]-2*g_E*Y[0];
}

void RK4 ( double t, double* Y, void (*RHSFunc)(double,double*,double*), double dt, int neq )
{
    int i;

    double Y1[neq], k1[neq], k2[neq], k3[neq], k4[neq];

    RHSFunc ( t, Y, k1 );

    for ( i=0; i<neq; i++ )
    {
        Y1[i] = Y[i]+dt*0.5*k1[i];
    }

    RHSFunc ( t+0.5*dt, Y1, k2 );

    for ( i=0; i<neq; i++ )
    {
        Y1[i] = Y[i]+0.5*dt*k2[i];
    }

    RHSFunc ( t+0.5*dt, Y1, k3 );
    
    for ( i=0; i<neq; i++ )
    {
        Y1[i] = Y[i]+dt*k3[i];
    }

    RHSFunc ( t+dt, Y1, k4 );

    for ( i=0; i<neq; i++ )
    {
        Y[i] = Y[i]+dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.;
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
        else if ( fxm==0 )
        {
            dx = 0.;
            zero_b = xm;
        }
        else
        {
            cout << "errore" << endl;
            dx = 0.;
            zero_b = xm;
        }
        
        kb++;

        cout.flush ();

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

double Residual ( double E )
{   
    int controllo;
    int n = 1200, neq = 2;
    
    double xl = -10., xr = 10., xm = 0.1;
    double x = 0., dx0, dx;
    double xml = 0., xmr = 0.;
    double Yf[neq], Yb[neq];                                           // f = forward e b = backward
    double A, B, R;

    dx0 = (xr-xl)/(double)n;
    g_E = E;
    x = xl;
    controllo = 1;

    Yf[0] = exp(-x*x*0.5);
    Yf[1] = -x*Yf[0];

    while ( controllo )
    {
        dx = dx0;

        if ( x+dx>xm )
        {
            dx = xm-x;
            controllo = 0;
        }

        RK4 ( x, Yf, dYdx, dx, neq );
        
        x += dx;
    }

    xml = x;
    x = xr;
    controllo = 1;

    Yb[0] = exp(-x*x*0.5);
    Yb[1] = -x*Yb[0];

    while ( controllo )
    {
        dx = -dx0;

        if ( dx<xm-x )
        {
            dx = xm-x;
            controllo = 0;
        }

        RK4 ( x, Yb, dYdx, dx, neq );

        x += dx;
    }

    xmr = x;

    if ( fabs(xmr-xml)>1.e-10 )
    {
        cout << "error at matching point" << endl;
    }
    
    A = Yf[1]*Yb[0];
    B = Yf[0]*Yb[1];
    R = (A-B)/sqrt(A*A+B*B)+1.e-7;
    
    return R;   
}