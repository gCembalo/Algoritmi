#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

double Function ( double );
double Rectangular ( double (*)(double), double, double, int );
double Midpoint ( double (*)(double), double, double, int );
double Trapezoidal ( double (*)(double), double, double, int );
double Simpson ( double (*)(double), double, double, int );

int main ()
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 12 );

    int n, nr = 2, nm = 2, nt = 2, ns = 2;
    
    double a = 0., b = 1., tol = 1.e-5;
    double i2_nr = 0., i2_nm = 0., i2_nt = 0., i2_ns = 0., i_nr, i_nm, i_nt, i_ns;                    

    cout << endl;
    cout << "Funzione: e^-x" << endl;
    cout << "--------------------------------------" << endl;
    cout << "inserisci il numero di sottointervalli" << endl;

    cin >> n;

    cout << "----------------------------------------------" << endl;
    cout << "exact: \t\t " << 1-exp(-1) << endl;                                                                                 // risultato del calcolo analitico
    cout << "rectangular: \t " << Rectangular ( Function, a, b, n ) << "\t N = " << n << endl;
    cout << "midpoint: \t " << Midpoint ( Function, a, b, n ) << "\t N = " << n << endl;
    cout << "trapezoidal: \t " << Trapezoidal ( Function, a, b, n ) << "\t N = " << n << endl;
    cout << "simpson: \t " << Simpson ( Function, a, b, n ) << "\t N = " << n << "\t\t (il metodo richiede N pari)" << endl;               
    cout << endl;

    do
    {
        nr = 2*nr;
        i2_nr = Rectangular ( Function, a, b, 2*nr );
        i_nr = Rectangular ( Function, a, b, nr );
    }
    while ( abs(i2_nr-i_nr)>=tol );

    do
    {
        nm = 2*nm;
        i2_nm = Midpoint ( Function, a, b, 2*nm );
        i_nm = Midpoint ( Function, a, b, nm );
    }
    while ( abs(i2_nm-i_nm)>=tol );

    do
    {
        nt = 2*nt;
        i2_nt = Trapezoidal ( Function, a, b, 2*nt );
        i_nt = Trapezoidal ( Function, a, b, nt );
    }
    while ( abs(i2_nt-i_nt)>=tol );

    do
    {
        ns = 2*ns;
        i2_ns = Simpson ( Function, a, b, 2*ns );
        i_ns = Simpson ( Function, a, b, ns );
    }
    while ( abs(i2_ns-i_ns)>=tol );

    cout << "Tolleranza: 10^-5" << endl;
    cout << "--------------------------------------------------" << endl;
    cout << "rectangular: \t " << i2_nr << "\t N = " << 2*nr << endl;
    cout << "midpoint: \t " << i2_nm << "\t N = " << 2*nm << endl;
    cout << "trapezoidal: \t " << i2_nt << "\t N = " << 2*nt << endl;
    cout << "simpson: \t " << i2_ns << "\t N = " << 2*ns << endl;
    cout << endl;

    return 0;
}

double Function ( double x )
{
    return exp(-x);
}

double Rectangular ( double (*F)(double), double a, double b, int n )                // funzione approssimata come costante con l'estremo sinistro
{  
    int i;
    
    double h, sum = 0.;

    h = (b-a)/(double)n;
  
    for ( i=0; i<n; i++ )
    {
        sum += F(a+h*i)*h;
    }

    return sum;
}

double Midpoint ( double (*F)(double), double a, double b, int n )                   // funzione approssimata come costante con il punto centrale
{  
    int i;
    
    double h, sum = 0.;

    h = (b-a)/(double)n;
  
    for ( i=0; i<n; i++ )
    {
        sum += F(a+h*(i+0.5))*h;
    }

    return sum;
}

double Trapezoidal ( double (*F)(double), double a, double b, int n )                // funzione approssimata come retta
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

double Simpson ( double (*F)(double), double a, double b, int n )                    // funzione approssimata come parabola
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