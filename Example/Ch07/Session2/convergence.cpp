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

    int i, j, neq = 2, nst = 4;

    double ti = 0., tf = 3., dt, t;
    double Y1[neq], Y2[neq], Y3[neq], Y4[neq];
    double err_eul, err_rk2m, err_rk2h, err_rk4;

    fdata.open ( "convergence.dat" );

    for ( j=0; j<9; j++ )
    {
        dt = (tf-ti)/(double)nst;
        t = ti;

        Y1[0] = Y2[0] = Y3[0] = Y4[0] = 1.;
        Y1[1] = Y2[1] = Y3[1] = Y4[1] = 0.;
        
        for ( i=0; i<nst; i++ )
        {
            EulerStep ( t, Y1, dYdt, dt, neq );
            RK2_midpoint ( t, Y2, dYdt, dt, neq );
            RK2_heun ( t, Y3, dYdt, dt, neq );
            RK4 ( t, Y4, dYdt, dt, neq );

            t += dt;
        }

        err_eul = fabs(Y1[0] - cos(tf));
        err_rk2m = fabs(Y2[0] - cos(tf));
        err_rk2h = fabs(Y3[0] - cos(tf));
        err_rk4 = fabs(Y4[0] - cos(tf));

        fdata << dt << " " << err_eul << " " << err_rk2m << " " << err_rk2h << " " << err_rk4 << endl;
        
        nst = 2*nst;
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