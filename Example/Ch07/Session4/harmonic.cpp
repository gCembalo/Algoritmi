#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

void Acceleration ( double*, double*, int );
void Velocity_verlet ( double, double*, double*, int, void (*)(double*,double*,int) );
void Position_verlet ( double, double*, double*, int, void (*)(double*,double*,int) );
void RK2_midpoint ( double, double*, void (*)(double,double*,double*), double, int );
void RK4 ( double, double*, void (*)(double,double*,double*), double, int );
void dYdt ( double, double*, double* );

int main ()
{
    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 6 );
    
    int i, neq = 2, nstep, np = 1;                                                                                        // np sono le particelle

    double omega = 1., E, Y[neq];
    double t, T = 2.*M_PI/omega, dt = 0.02*T, ti = 0., tf = 10*T;
    double x0 = 1., v0 = 0., x[1], v[1];

    nstep = (tf-ti)/dt;
    t = ti;

    Y[0] = x0;
    Y[1] = v0;

    fdata.open ( "harmonic.dat" );

    for ( i=0; i<nstep; i++ )
    {
        RK2_midpoint ( t, Y, dYdt, dt, neq );
        
        t = t+dt;
        
        E = 0.5*Y[0]*Y[0]+0.5*Y[1]*Y[1]; 
        
        fdata << t << " " << Y[0] << " " << Y[1] << " " << E << endl;
    }

    fdata << endl << endl;

    t = ti;

    Y[0] = x0;
    Y[1] = v0;

    for ( i=0; i<nstep; i++ )
    {
        RK4 ( t, Y, dYdt, dt, neq );
        
        t = t+dt;
        
        E = 0.5*Y[0]*Y[0]+0.5*Y[1]*Y[1]; 
        
        fdata << t << " " << Y[0] << " " << Y[1] << " " << E << endl;
    }

    fdata << endl << endl;

    t = ti;

    x[0] = x0;
    v[0] = v0;

    for ( i=0; i<nstep; i++ )
    {
        Position_verlet ( dt, x, v, np, Acceleration );
        
        t = t+dt;
        
        E = 0.5*x[0]*x[0]+0.5*v[0]*v[0]; 
        
        fdata << t << " " << x[0] << " " << v[0] << " " << E << endl;
    }

    fdata << endl << endl;

    t = ti;

    x[0] = x0;
    v[0] = v0;

    for ( i=0; i<nstep; i++ )
    {
        Velocity_verlet ( dt, x, v, np, Acceleration );
        
        t = t+dt;
        
        E = 0.5*x[0]*x[0]+0.5*v[0]*v[0]; 
        
        fdata << t << " " << x[0] << " " << v[0] << " " << E << endl;
    }

    fdata.close ();

    cout << endl;
    cout << " (vedi file associato) " << endl;
    cout << endl;

    return 0;
}

void Acceleration ( double* x, double* a, int n )
{
    int i;

    for ( i=0; i<n; i++ )
    {
        a[i] = -x[i];
    }
}

void Velocity_verlet ( double h, double* x, double* v, int np, void (*acceleration)(double*,double*,int) )  
{
    int k;

    double a[np];

    acceleration ( x, a, np );

    for ( k=0; k<np; k++ )
    {
        v[k] = v[k]+0.5*h*a[k];
    }

    for ( k=0; k<np; k++ )
    {
        x[k] = x[k]+h*v[k];
    }

    acceleration ( x, a, np );

    for ( k=0; k<np; k++ )
    {
        v[k] = v[k]+0.5*h*a[k];
    }
}

void Position_verlet ( double h, double* x, double* v, int np, void (*acceleration)(double*,double*,int) )    
{
    int k;

    double a[np];

    for ( k=0; k<np; k++ )
    {
        x[k] = x[k]+0.5*h*v[k];
    }

    acceleration ( x, a, np );

    for ( k=0; k<np; k++ )
    {
        v[k] = v[k]+h*a[k];
    }

    for ( k=0; k<np; k++ )
    {
        x[k] = x[k]+0.5*h*v[k];
    }
}

void dYdt ( double t, double* Y, double* R )
{
    double y1 = Y[0], y2 = Y[1];

    R[0] = y2;
    R[1] = -y1;
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