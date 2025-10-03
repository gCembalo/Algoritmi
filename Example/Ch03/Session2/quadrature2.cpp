#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

double Function1 ( double );
double Function2 ( double );
double Simpson ( double (*)(double), double, double, int );
double Gauss ( double (*)(double), double, double, int, int );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 12 );

    int n, ng;
    
    double a, b;

    cout << endl;
    cout << "Funzione: sqrt(1+x)" << endl;
    cout << "---------------------------------------------" << endl;
    cout << "exact: \t\t" << 4.666666666667 << endl;
    cout << "simpson: \t" << Simpson ( Function1, 0., 3., 2 ) << endl;
    cout << "gauss: \t\t" << Gauss ( Function1, 0., 3., 1, 3 ) << endl;
    cout << endl;
    cout << "Funzione: 1-x+2x^2+(1/2)x^3+(1/4)x^4-(1/8)x^5" << endl;
    cout << "---------------------------------------------" << endl;
    cout << "exact: \t\t" << -66./5. << endl;
    cout << "simpson: \t" << Simpson ( Function2, -1., 5., 2 ) << endl;
    cout << "gauss: \t\t" << Gauss ( Function2, -1., 5., 1, 3 ) << endl;
    cout << endl;

    return 0;
}

double Function1 ( double x )
{
    return sqrt(1+x);
}

double Function2 ( double x )
{
    return 1-x+2*x*x+0.5*x*x*x+0.25*x*x*x*x-0.125*x*x*x*x*x;
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
    double w[ng], x[ng];                                                                // l'array richiede l'utilizzo di un ng specifico a priori

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
    
    for ( i=0; i<n; i++ )                                                               // somma sugli intervalli
    {
        x0 = a+i*h; 
        x1 = x0+h;
        c1 = (x1-x0)*0.5;
        c2 = (x0+x1)*0.5;
        sumk = 0.;

        for ( k=0; k<ng; k++ )                                                          // somma sugli intervalli gaussiani
        {
            sumk += c1*w[k]*F(c1*x[k]+c2);
        }

        sum += sumk;
    }
    
    return sum;
}