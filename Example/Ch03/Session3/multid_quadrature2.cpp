#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>                                                                                 // fare operazioni con dei file (ora non usato)

using namespace std;

static int stat_nig = 1;                                                                           // numero di intervalli gaussiani 
static double stat_yp, stat_yg;                                                                    // static per tutte le variabili globali                                                          

double GFunc1D_p ( double );                                                                       // funzione 1D in y con somma sulle x
double GFunc1D_g ( double );
double FFunc1D_p ( double );                                                                       // funzione 1D in x con y fissato
double FFunc1D_g ( double );
double Func2D_p ( double, double );                                                                // funzione 2D da integrare
double Func2D_g ( double, double ); 
double Gauss ( double (*)(double), double, double, int, int );                                     // codice preesistente per l' integrale 1D 

int main ()                                                                                        // somma sulle y
{
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 12 );

    double yi = -1., yf = 1.;                                      

    cout << endl;
    cout << "Funzione polinomiale" << endl;
    cout << "integrale calcolato: \t" << Gauss ( GFunc1D_p, yi, yf, stat_nig, 3 ) << endl;
    cout << "integrale vero: \t" << 412./45. << endl;
    cout << endl;
    cout << "Funzione a gradino" << endl;
    cout << "integrale calcolato: \t" << Gauss ( GFunc1D_g, yi, yf, stat_nig, 4 ) << endl;
    cout << "integrale vero: \t" << M_PI << endl;

    double d = 1., tol = 10e-5;

    while ( d>tol )
    {
        d = abs ( Gauss ( GFunc1D_g, yi, yf, stat_nig, 4 ) - M_PI );
        stat_nig++;
    }

    cout << "precisione \t\t" << tol << "\t N = " << stat_nig << endl;
    cout << endl;

    return 0;
}

double GFunc1D_p ( double y )                                                                                                                                
{
    double sum_x = 0., xi = -1., xf = 1.;                                                   

    stat_yp = y;                                                                           

    sum_x = Gauss ( FFunc1D_p, xi, xf, stat_nig, 3 );                                       
    
    return sum_x;
}

double GFunc1D_g ( double y )                                                                                                                                
{
    double sum_x = 0., xi = -1., xf = 1.;                                                   

    stat_yg = y;                                                                           

    sum_x = Gauss ( FFunc1D_g, xi, xf, stat_nig, 4 );                                       
    
    return sum_x;
}

double FFunc1D_p ( double x )                                                              
{
    return Func2D_p ( x, stat_yp );
}

double FFunc1D_g ( double x )                                                              
{
    return Func2D_g ( x, stat_yg );
}

double Func2D_p ( double x, double y )                                                             // p = polinomio                                    
{
    return x*x*x*x*y*y+2*x*x*y*y-y*x*x+2;
}

double Func2D_g ( double x, double y )                                                             // g = gradino
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