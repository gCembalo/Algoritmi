// Consider the function:
// exp( 1/( x + 0.5 ) ) - ( 3 + 2*x )/( 1 + x )
// Search the zero over the interval [0,2] with tolerance 10-7. The solution is x0 ≈ 5.235934e-01 Any problem ? Explain.
// Bisection: 26
// False position: 55
// Secant: ??
// Newton: 8
//

// includo le funzioni di root finding
#include "root_finders.h"
// definisco le funzioni che studio
double func(double);
double derfunc(double);

int main(){

    // definisco gli estremi dell'intervallo
    double a = 0, b = 2;
    // defiisco una tolleranza
    double tol = 1.e-8;
    //double facc = 1e-6;
    // definisco gli zeri
    double x0, x1, x2, x3;

    cout << "\n+-----------------------------------------------------------------------------+\nLo zero è a x = " << 5.235934e-01 << endl;

    cout << "prendendo l'intervallo [ " << a << " , " << b << " ]\n" << endl;

    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    cout << "\nCon # iterazioni:\n" << endl;
     // stampo gli algoritmi
    bisection(func, a, b, tol, x0);
    false_position(func, a, b, tol, x1);
    secant_method(func, a, b, tol, x2);
    newton_method(func, derfunc, a, tol, x3);

    cout << "\nI diversi metodi restituiscono:\n" << endl;

    cout << "(Bisezione): " << x0 << "\n(False position): " << x1 << "\n(Secant): " << x2 << "\n(Newton): " << x3 << endl;

    cout << "\n+-----------------------------------------------------------------------------+\n" << endl;

    return 0;
}




double func(double x){

    return exp( 1/( x + 0.5 ) ) - ( 3 + 2*x )/( 1 + x );

}

double derfunc(double x){

    return - 1 / ((x + 0.5)*(x+0.5)) * exp( 1 / (x + 0.5)) + 1 / ((1+x) * (1+x));

}