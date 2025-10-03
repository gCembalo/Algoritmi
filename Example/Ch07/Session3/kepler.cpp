#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

void RK4 ( double, double*, void (*)(double,double*,double*), double, int );
void dYdt ( double, double*, double* );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 6 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 6 );

    int i, count = 0;
    int neq = 4, n = 10;                                                                    // n sono le semiorbite

    double dt, t = 0., alpha = 1., E;
    double r, x0 = 4., y0 = 0., Y[neq];
    double v, vx0 = 0., vy0, vy_old;
    
    r = sqrt(x0*x0+y0*y0);
    vy0 = sqrt(alpha/r);
    v = sqrt(vx0*vx0+vy0*vy0);
    E = 0.5*v*v-1./r;
    
    cout << endl;
    cout << "E = " << E << endl;
    cout << "v = " << v << endl;
    cout << endl;

    if ( v==sqrt(2./r) )
    {
        cout << "orbita parabolica" << endl;
    }
    else if ( v>sqrt(2./r) )
    {
        cout << "orbita iperbolica" << endl;
    }
    else if ( v==sqrt(1./r) )
    {
        cout << "orbita circolare" << endl;
    }
    else
    {
        cout << "orbita ellittica" << endl;
    }

    cout << endl;

    Y[0] = x0;
    Y[1] = y0;
    Y[2] = vx0;
    Y[3] = vy0;
    
    fdata.open ( "kepler.dat" );

    for ( i=0; i<10000; i++ )
    {
        dt = 0.1*r/v;                                                                      // dt adattivo: frazione di r/v 
        vy_old = Y[3];
        t = t+dt;
            
        RK4 ( t, Y, dYdt, dt, neq );

        if ( Y[3]*vy_old<0 )                                                               // conteggio delle orbite con i turning points
        {       
            count++;

            if ( count==2*n+1 ) break;
        }
    
        r = sqrt(Y[0]*Y[0]+Y[1]*Y[1]);
        v = sqrt(Y[2]*Y[2]+Y[3]*Y[3]);
        
        fdata << Y[0] << " " << Y[1] << " " << Y[2] << " " << Y[3] << endl;
    }

    fdata.close();

    return 0;
}

void dYdt ( double t, double* Y, double* R )
{
    double x, y, vx, vy, r;
    
    x = Y[0];
    y = Y[1];
    vx = Y[2];
    vy = Y[3]; 

    r = sqrt(x*x+y*y);

    R[0] = vx;
    R[1] = vy;
    R[2] = -x/(r*r*r); 
    R[3] = -y/(r*r*r);
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