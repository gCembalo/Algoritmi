#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

const double alpha = 10;

double Funcx ( double );

int main ()
{
    int i, n = 10000;
    
    double t, ti = 0., tf = alpha;
    double h, v, a;
    double f, fp, fm, fpp, fmm;

    h = (tf-ti)/(double)n;
    
    fdata.open ( "trajectory.dat" );

    for ( i=0; i<=n; i++ )
    {
        t = ti+i*h;
        f = Funcx(t);
        fp = Funcx(t+h);
        fm = Funcx(t-h);
        fpp = Funcx(t+2*h);
        fmm = Funcx(t-2*h);

        if ( t==0 )
        {
            v = (fp-f)/h;
            a = (fpp-2*fp+f)/(h*h);
        }
        else
        {
            v = (fp-fm)/(2*h);
            a = (fp-2*f+fm)/(h*h);
        }
        
        fdata << t << " " << Funcx(t) << " " << v << " " << a << endl;
    }

    fdata.close();

    cout << endl;
    cout << " (vedi file associato) " << endl;
    cout << endl;
    
    return 0;
}

double Funcx ( double t )
{
    if ( t==0 ) 
    {
        return 0.;
    }
    else 
    {
        return alpha*t*t-t*t*t*(1-exp(-alpha*alpha/t));
    }
}