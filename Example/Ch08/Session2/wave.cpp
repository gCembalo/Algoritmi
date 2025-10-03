#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

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

    int i, k, err;
    int n = 50, nr = 10;
    int nit, nit_max = 99;
    int neq = 3, nst = 100;

    double xi = 0., xf = 1., x, dx;
    double Y[neq], s = 1., zero;
    double a = 1., b = 20.;
    double kl[nr], kr[nr];
    double xtol = 1.e-9, ytol = 1.e-7;

    dx = (xf-xi)/(double)nst;

    fdata.open ( "wave.dat" );

    for ( k=1; k<=5; k++ )
    {
        x = xi;

        Y[0] = 0.;
        Y[1] = s;
        Y[2] = (double)k;

        fdata << x << " " << Y[0] << " " << Y[1] << endl;
        
        for ( i=0; i<nst; i++ )
        {
            RK4 ( x, Y, dYdx, dx, neq );

            x += dx;

            fdata << x << " " << Y[0] << " " << Y[1] << endl;
        }

        fdata << endl << endl;
    }

    fdata.close ();

    cout << endl;

    Bracket ( Residual, a, b, kl, kr, n, nr );

    for ( i=0; i<nr; i++ )
    {
        err = Bisection ( Residual, kl[i], kr[i], xtol, ytol, zero, nit, nit_max );
 
        if ( err==0 )
        {
            cout << "x" << i+1 << " = " << zero << "\t nit = " << nit << endl;
        }
        else 
        {
            cout << "Nessuna soluzione" << endl;
        }
        
        nit = 0;
    }
    
    cout << endl;
    
    return 0;
}

void dYdx ( double r, double* Y, double* R )
{
    R[0] = Y[1];
    R[1] = -Y[2]*Y[2]*Y[0];   
    R[2] = 0.;
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

double Residual ( double k )
{
    int i, nst = 100, neq = 3;

    double xi = 0., xf = 1., x, dx;
    double Y[neq], res;

    x = xi;
    dx = (xf-xi)/(double)nst;
    
    Y[0] = 0.;
    Y[1] = 1.;
    Y[2] = k;
    
    for ( i=0; i<nst; i++ )
    {
        RK4 ( x, Y, dYdx, dx, neq );

        x += dx;
    }
    
    res = Y[0]-0.;

    return res;  
}