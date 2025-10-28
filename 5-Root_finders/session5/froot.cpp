// Modify your program to find all of the zeros of the function:
// sin(x) - ( x*x/100. + x/5. + 1./3. )
// by first using a bracketing function and then any of the root finder methods to refine the search. This function has five roots in the interval [-10,10], corresponding to:
// -8.716925e+00
// -6.889594e+00
// -2.968485e+00
// 4.361680e-01
// 2.183971e+00
// 
// Are you able to find them all ? Which root-finder performs the best ?
// 

#include "root_finders.h"

double func(double);
double derfunc(double);

int main(){

    // definisco gli estremi dell'intervallo
    double a = -10, b = 10;
    // defiisco una tolleranza
    double tol = 1.e-8;
    //double facc = 1e-6;

    cout << "\n+-----------------------------------------------------------------------------+\nLa funzione che vogliamo studiare è:\n        sin(x) - ( x*x/100 + x/5 + 1/3 )" << endl;
    cout << "\nGli zeri sono: \n       " << -8.716925e+00 << "      " << -6.889594e+00 << "      " << -2.968485e+00 << "      " << 4.361680e-01 << "      " << 2.183971e+00 << "      \n" << endl;
    
    int n;
    cout << "In quanti intervalli vuoi dividere [ " << a << " , " << b << " ] ?" << endl;
    cin >> n;

    // definisco il numero di roots che avremo
    int nroots = 0;
    // definisco gli arrai degli estremi
    double xL[n*n], xR[n*n];
    // definisco un array dove mettere gli zeri
    double x0[n*n];

    // invoco il bracketing
    Bracket(func,a,b,n,xL,xR,nroots);

    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    // stampo a terminale gli estremi di ogni intervallo
    cout << "\n+-----------------------------------------------------------------------------+" << endl;
    cout << "Le radici si trovano negli intervalli:\n" << endl;
    for( int j = 0 ; j < nroots ; j++){

        cout << "n = " << j << ":       [ " << xL[j] << "   ,   " << xR[j] << " ]" << endl;

    }

    // definisco la variabile che mi da le iterazioni
    int l = 0;
    // stampo a terminale i risultati dei diversi metodi per trovare le radici

    cout << "\n+-----------------------------------------------------------------------------+" << endl;
    cout << "Bisezione:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        bisection(func, xL[k], xR[k], tol, x0[k], l);
        cout << "x" << k << " = " << x0[k] << "     trovata con " << l << " iterazioni" << endl;
    }

    cout << "\n+-----------------------------------------------------------------------------+" << endl;
    cout << "False position:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        false_position(func, xL[k], xR[k], tol, x0[k], l);
        cout << "x" << k << " = " << x0[k] << "     trovata con " << l << " iterazioni" << endl;
    }

    cout << "\n+-----------------------------------------------------------------------------+" << endl;
    cout << "Secante:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        secant_method(func, xL[k], xR[k], tol, x0[k], l);
        cout << "x" << k << " = " << x0[k] << "     trovata con " << l << " iterazioni" << endl;
    }

    cout << "\n+-----------------------------------------------------------------------------+" << endl;
    cout << "\nNewton:" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        newton_method(func, derfunc, xL[k], tol, x0[k], l);
        cout << "x" << k << " = " << x0[k] << "     trovata con " << l << " iterazioni" << endl;
    }

    return 0;
}



double func(double x){
    
    return sin(x) - ( x*x/100. + x/5. + 1./3. );

}

double derfunc(double x){
    
    return cos(x) - ( x/50. + 1./5. );

}