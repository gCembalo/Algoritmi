#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

void dYdr ( double, double*, double* );
void RK4 ( double, double*, void (*)(double,double*,double*), double, int );

int Bisection ( double (*)(double), double, double, double, double, double&, int&, int );

double Residual ( double );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 8 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 8 );

    int i, j;
    int err, k = 0, kmax = 99;
    int neq = 2, nst = 1000;

    double s = 0., si = 0., sf = 5., ds; 
    double ra = 0., rb = 20., r = ra, dr;
    double Y[neq], Yex, res;
    double a = 0.1, b = 1.;
    double x_tol = 1.e-9, y_tol = 1.e-6; 

    dr = (rb-ra)/(double)nst;

    fdata.open ( "poisson.dat" );

    for ( j=0; j<=10; j+=2 )
    {
        s = j*0.1;
        r = ra;

        Y[0] = 0.;                                                              // phi(r)
        Y[1] = s;                                                               // d[phi(r)]/dr

        fdata << r << " " << Y[0] << " " << Y[1] << endl;
        
        for ( i=0; i<nst; i++ )
        {
            RK4 ( r, Y, dYdr, dr, neq );

            r += dr;

            fdata << r << " " << Y[0] << " " << Y[1] << endl;
        }

        fdata << endl << endl;
    }

    ds = (sf-si)/(double)nst;

    for ( i=0; i<nst; i++ )
    {
        res = Residual(s);

        s = si+i*ds;

        fdata << s << " " << res << endl;
    }

    fdata << endl << endl;

    err = Bisection ( Residual, a, b, x_tol, y_tol, s, k, kmax );

    cout << endl;
    
    if ( err==0 )
    {
        cout << "s = " << s << "\t k = " << k << endl;
    }
    else 
    {
        cout << "Nessuna soluzione" << endl;
    }
    
    cout << endl;

    r = 0., rb = 50.;

    Y[1] = 0.5;

    for ( i=1; i<=nst; i++ )
    {
        RK4 ( r, Y, dYdr, dr, neq );

        r += dr;
        Yex = 1-0.5*(r+2)*exp(-r);
        
        fdata << r << " " << Y[0]/r << " " << Yex/r << endl;
    }

    fdata.close ();

    return 0;
}

void dYdr ( double r, double* Y, double* R )
{
    double v, rho;
    
    v = Y[1];
    rho = exp(-r)/(8.*M_PI);

    R[0] = v;
    R[1] = -4.*M_PI*r*rho;
}

void RK4 ( double t, double* Y, void (*RHSFunc)(double,double*,double*), double dt, int neq )
{
    int i;

    double Y1[neq], k1[neq], k2[neq], k3[neq], k4[neq];

    RHSFunc ( t, Y, k1 );

    for ( i=0; i<neq; i++ )
    {
        Y1[i] = Y[i]+0.5*dt*k1[i];
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

        cout.flush ();                                                          // svuota il buffer intasato

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

double Residual ( double s )
{   
    int i, n = 1000, neq = 2;

    double ri = 0., rf = 20., dr, r;
    double res, Y[neq];

    dr = (rf-ri)/(double)n;
    
    Y[0] = 0.;
    Y[1] = s;
    
    for ( i=0; i<n; i++ )
    {
        RK4 ( r, Y, dYdr, dr, neq );

        r += dr;
    }

    res = Y[0]-1;

    return res;
}