#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

void EulerStep ( double, double*, void (*)(double,double*,double*), double, int );
void dYdt ( double, double*, double* );

double Funcy ( double );

int main ()
{
    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 6 );

    int i, k;
    int neq = 1, nst;

    double t, dt, ti = 0., tf = 3.;
    double f, Y[neq], st[10] = {1.,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001};
    double abs_err, rel_err;

    fdata.open ( "ode1.dat" );

    for ( k=0; k<10; k++ )
    {
        Y[0] = 1.;
        
        dt = st[k];
        t = ti;
        nst = (tf-ti)/dt;
        f = Funcy(t);
        abs_err = fabs(Y[0]-f);
        rel_err = abs_err/f;

        fdata << t << " " << Y[0] << " " << abs_err << " " << rel_err << endl;

        for ( i=0; i<nst; i++ )
        {
            EulerStep ( t, Y, dYdt, dt, neq );

            t += dt;
            f = Funcy(t);
            abs_err = fabs(Y[0]-f);
            rel_err = abs_err/f;

            fdata << t << " " << Y[0] << " " << abs_err << " " << rel_err << endl;
        }

        fdata << endl << endl;
    }

    cout << endl;
    cout << " (vedi file associato) " << endl;
    cout << endl;

    return 0;
}

double Funcy ( double t )
{
  return exp(-t*t*0.5);
}

void dYdt ( double t, double* Y, double* R )
{
    double y;
    
    y = Y[0];
    
    R[0] = -t*y;
}

void EulerStep ( double t, double* Y, void (*RHSFunc)(double,double*,double*), double dt, int neq )    // prende un array Y e lo avanza al tempo t
{
    int k;
    
    double rhs[neq];

    RHSFunc ( t, Y, rhs );

    for ( k=0; k<neq; k++ )
    {
        Y[k] += dt*rhs[k];
    }
}