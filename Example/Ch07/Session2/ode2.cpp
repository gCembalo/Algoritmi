#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

void EulerStep ( double, double*, void (*)(double,double*,double*), double, int );
void RK2_midpoint ( double, double*, void (*)(double,double*,double*), double, int );
void RK2_heun ( double, double*, void (*)(double,double*,double*), double, int );
void RK4 ( double, double*, void (*)(double,double*,double*), double, int );
void dYdt ( double, double*, double* );

int main ()
{
    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 6 );
    
    int i, neq = 2, nst = 200;

    double ti = 0., tf = 20*M_PI, dt, t, Y[neq];

    fdata.open ( "ode2.dat" );

    dt = (tf-ti)/(double)nst;
    t = ti;

    Y[0] = 1.;
    Y[1] = 0.;

    fdata << t << " " << Y[0] << " " << Y[1] << endl;

    for ( i=0; i<nst; i++ )
    {
        EulerStep ( t, Y, dYdt, dt, neq );

        t += dt;
        
        fdata << t << " " << Y[0] << " " << Y[1] << endl;
    }

    fdata << endl << endl;

    t = ti;
    
    Y[0] = 1.;
    Y[1] = 0.;

    fdata << t << " " << Y[0] << " " << Y[1] << endl;

    for ( i=0; i<nst; i++ )
    {
        RK2_midpoint ( t, Y, dYdt, dt, neq );

        t += dt;
        
        fdata << t << " " << Y[0] << " " << Y[1] << endl;
    }

    fdata << endl << endl;

    t = ti;

    Y[0] = 1.;
    Y[1] = 0.;

    fdata << t << " " << Y[0] << " " << Y[1] << endl;

    for ( i=0; i<nst; i++ )
    {
        RK2_heun ( t, Y, dYdt, dt, neq );

        t += dt;
        
        fdata << t << " " << Y[0] << " " << Y[1] << endl;
    }

    fdata << endl << endl;

    t = ti;

    Y[0] = 1.;
    Y[1] = 0.;

    fdata << t << " " << Y[0] << " " << Y[1] << endl;

    for ( i=0; i<nst; i++ )
    {
        RK4 ( t, Y, dYdt, dt, neq );

        t += dt;
        
        fdata << t << " " << Y[0] << " " << Y[1] << endl;
    }

    fdata.close ();

    cout << endl;
    cout << " (vedi file associato) " << endl;
    cout << endl;

    return 0;
}

void dYdt ( double t, double* Y, double* R )
{
    double x, y;
    
    x = Y[0]; 
    y = Y[1];

    R[0] = y; 
    R[1] = -x;
}

void EulerStep ( double t, double* Y, void (*RHSFunc)(double,double*,double*), double dt, int neq )
{
    int k;

    double rhs[neq];

    RHSFunc ( t, Y, rhs );

    for ( k=0; k<neq; k++ )
    {
        Y[k] += dt*rhs[k];
    }
}

void RK2_midpoint ( double t ,double* Y, void (*RHSFunc)(double,double*,double*), double dt, int neq )
{
    int i;

    double Y1[neq], k1[neq], k2[neq];

    RHSFunc ( t, Y, k1 );

    for ( i=0; i<neq; i++ )
    {
        Y1[i] = Y[i]+dt*0.5*k1[i];
    }

    RHSFunc ( t+0.5*dt, Y1, k2 );

    for ( i=0; i<neq; i++ )
    {
        Y[i] = Y[i]+dt*k2[i];
    }
}

void RK2_heun ( double t ,double* Y, void (*RHSFunc)(double,double*,double*), double dt, int neq )
{
    int i;

    double Y1[neq], k1[neq], k2[neq];

    RHSFunc ( t, Y, k1 );

    for ( i=0; i<neq; i++ )
    {
        Y1[i] = Y[i]+dt*k1[i];
    }

    RHSFunc ( t+dt, Y1, k2 );

    for ( i=0; i<neq; i++ )
    {
        Y[i] = Y[i]+dt*0.5*(k1[i]+k2[i]);
    }
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