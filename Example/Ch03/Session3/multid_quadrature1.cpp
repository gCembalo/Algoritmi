#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

double Func2D_p ( double, double );
double Func2D_g ( double, double );
double Gauss2D ( double (*)(double,double), double, double, double, double, int, int );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 12 );  

    int n = 1;
    
    double d = 1., tol = 10e-5;

    cout << endl;
    cout << "Funzione polinomiale" << endl;
    cout << "integrale calcolato: \t" << Gauss2D ( Func2D_p, -1, 1, -1, 1, 1, 3 ) << endl;
    cout << "integrale vero: \t" << 412./45. << endl;
    cout << endl;
    cout << "Funzione a gradino" << endl;
    cout << "integrale calcolato: \t" << Gauss2D ( Func2D_g, -1, 1, -1, 1, 1, 4 ) << endl;
    cout << "integrale vero: \t" << M_PI << endl;

    while ( d>tol )
    {
        d = abs ( Gauss2D ( Func2D_g, -1, 1, -1, 1, n, 4 ) - M_PI );
        n++;
    }

    cout << "precisione: \t\t" << tol << "\t N = " << n << endl;
    cout << endl;

    return 0;
}

double Func2D_p ( double x, double y )                                                    
{
    return x*x*x*x*y*y+2*x*x*y*y-y*x*x+2;
}

double Func2D_g ( double x, double y )                                                    
{
    if ( sqrt(x*x+y*y)<=1 )
    {
        return 1.;
    }
    else
    {
        return 0.;
    }
}

double Gauss2D ( double (*F)(double,double), double xa, double xb, double ya, double yb, int n, int ng )
{
    double w[ng], x[ng];

    if ( ng==1 ) 
    {
        x[0] = 0.;

        w[0] = 2.;
    }
    else if ( ng==2 ) 
    {
        x[0] = -sqrt(1./3.);
        x[1] = sqrt(1./3.);

        w[0] = w[1] = 1.;
    }
    else if ( ng==3 ) 
    {
        x[0] = -sqrt(3./5.);
        x[1] = 0.;
        x[2] = sqrt(3./5.);

        w[1] = 8./9.;
        w[0] = w[2] = 5./9.;
    }
    else if ( ng==4 ) 
    {
        x[0] = -sqrt(3./7.+2./7.*sqrt(6./5.));
        x[1] = -sqrt(3./7.-2./7.*sqrt(6./5.));
        x[2] = sqrt(3./7.-2./7.*sqrt(6./5.));
        x[3] = sqrt(3./7.+2./7.*sqrt(6./5.));

        w[0] = w[3] = (18-sqrt(30))/36.;
        w[1] = w[2] = (18+sqrt(30))/36.;
    }
    else if ( ng==5 ) 
    {
        x[0] = -1./3.*sqrt(5+2*sqrt(10./7.));
        x[1] = -1./3.*sqrt(5-2*sqrt(10./7.));
        x[2] = 0.;
        x[3] = 1./3.*sqrt(5-2*sqrt(10./7.));
        x[4] = 1./3.*sqrt(5+2*sqrt(10./7.));
        
        w[0] = w[4] = (322-13*sqrt(70))/900.;
        w[1] = w[3] = (322+13*sqrt(70))/900.;
        w[2] = 128./225.;
    }
    else
    {
        cout << "codice non implementato per Ng>5" << endl;
    }

    int i, ik, j, jk;
    
    double sum = 0., sum_x, sum_y;
    double xc, yc, dx, dy;

    dx = (xb-xa)/(double)n;
    dy = (yb-ya)/(double)n;
    
    for ( j=0; j<n; j++ )                                                                                           // sottointervalli y
    {
        for ( i=0; i<n; i++ )                                                                                       // sottointervalli x
        {
            xc = xa+(i+0.5)*dx;                                                                                     // centro intervallo x
            yc = ya+(j+0.5)*dy;                                                                                     // centro intervallo y

            sum_y = 0.;
            
            for ( jk=0; jk<ng; jk++ )
            {
                sum_x = 0.;

                for ( ik=0; ik<ng; ik++ ) 
                {
                    sum_x += w[ik]*F ( xc+0.5*dx*x[ik], yc+0.5*dy*x[jk] );
                }
                
                sum_x *= 0.5*dx;
                sum_y += w[jk]*sum_x;
            }
            
            sum_y *= 0.5*dy;
            sum += sum_y;
        }
    }

    return sum;
}