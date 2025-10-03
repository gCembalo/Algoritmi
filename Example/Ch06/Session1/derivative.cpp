#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream fdata;

double Func ( double );
double Der ( double );

int main ()
{
    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 10 );

    int i;
    
    double x = 1., h = 0.5;
    double fd, bd, cd, ho;
    double err_fd, err_bd, err_cd, err_ho;
    double f = Func(x), fex = Der(x);

    fdata.open ( "derivative.dat" );

    for ( i=0; i<10; i++ )
    {
        double  fp, fm, fpp, fmm;
        
        fp = Func(x+h);
        fm = Func(x-h);
        fpp = Func(x+2*h);
        fmm = Func(x-2*h);

        fd = (fp-f)/h;
        bd = (f-fm)/h;
        cd = (fp-fm)/(2*h);
        ho = (fmm-8*fm+8*fp-fpp)/(12*h);

        err_fd = fabs(fd-fex);
        err_bd = fabs(bd-fex);
        err_cd = fabs(cd-fex);
        err_ho = fabs(ho-fex);
        
        fdata << 1/h << " " << err_fd << " "  << err_bd << " " << err_cd << " " << err_ho << endl;
        
        h = h*0.5;
    }

    fdata.close ();

    cout << endl;
    cout << " (vedi file associato) " << endl;
    cout << endl;
    
    return 0;
}

double Func ( double x )
{
    return sin(x);
}

double Der ( double x )
{
    return cos(x);
}