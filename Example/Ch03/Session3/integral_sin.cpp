#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

double Function ( double );
double Trapezoidal ( double (*)(double), double, double, int );
double Simpson ( double (*)(double), double, double, int );
double Gauss ( double (*)(double), double, double, int, int );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 6 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 6 );

    int i, n, ng;
    
    double a = 0., b = 0.8;

    cout << endl;
    cout << "Funzione: \t sin(x)/x" << endl;
    cout << "exact: \t\t" << 0.7720957855 << endl;
    cout << "----------------------------------------------" << endl;
    cout << "trapezoidal: \t" << Trapezoidal ( Function, a, b, 1 ) << "\t N = 1" << endl;
    cout << "trapezoidal: \t" << Trapezoidal ( Function, a, b, 2 ) << "\t N = 2" << endl;
    cout << "trapezoidal: \t" << Trapezoidal ( Function, a, b, 4 ) << "\t N = 4" << endl;
    cout << "trapezoidal: \t" << Trapezoidal ( Function, a, b, 8 ) << "\t N = 8" << endl;
    cout << "----------------------------------------------" << endl;
    cout << "simpson: \t" << Simpson ( Function, a, b, 2 ) << "\t N = 2" << endl;
    cout << "simpson: \t" << Simpson ( Function, a, b, 4 ) << "\t N = 4" << endl;
    cout << "simpson: \t" << Simpson ( Function, a, b, 8 ) << "\t N = 8" << endl;
    cout << "----------------------------------------------" << endl;
    cout << "gauss: \t\t" << Gauss ( Function, a, b, 1, 3 ) << "\t N = 1  Ng = 3" << endl;
    cout << "----------------------------------------------" << endl;
    cout << "asymptote: \t" << Gauss ( Function, a, 1000, 10000, 3 ) << endl;
    cout << endl;

    fdata.open ( "integral_sin.dat" );

    for ( i=0; i<=250; i++ )
    {
        fdata << i*0.1 << " " << Gauss ( Function, a, i*0.1, i, 3 ) << endl;
    }

    fdata.close ();

    return 0;
}

double Function ( double x )                                                          // tutte le singolarità vanno rimosse manualmente
{
    if ( x==0 ) 
    {
        return 1.;
    }
    else
    {
        return sin(x)/x;
    }
}

double Trapezoidal ( double (*F)(double), double a, double b, int n )
{
    int i;
    
    double h, sum;

    h = (b-a)/(double)n;
    sum = 0.5*(F(a)+F(b))*h;
  
    for ( i=1; i<n; i++ )
    {
        sum += F(a+i*h)*h;
    }

    return sum;
}

double Simpson ( double (*F)(double), double a, double b, int n )
{
    int i, w = 4;
    
    double h, sum;

    h = (b-a)/(double)n;
    sum = (F(a)+F(b))*h/3;

    for ( i=1; i<n; i++ )
    {
        sum += w*F(a+i*h)*h/3;
        w = 6-w;
    }

    return sum;
}

double Gauss ( double (*F)(double), double a, double b, int n, int ng )
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

    int i, k;
    
    double sum = 0., sumk;
    double x0, x1, c1, c2, h;

    h = (b-a)/(double)n;
    
    for ( i=0; i<n; i++ )                                                          
    {
        x0 = a+i*h; 
        x1 = x0+h;
        c1 = (x1-x0)*0.5;
        c2 = (x0+x1)*0.5;
        sumk = 0.;

        for ( k=0; k<ng; k++ )                                                    
        {
            sumk += c1*w[k]*F(c1*x[k]+c2);
        }

        sum += sumk;
    }
    
    return sum;
}