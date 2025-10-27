// using single precision arithmetic, obtain a numerical approximation to sqrt(x^2 + 1) – x (valid for large x) and 1 - cos(x) (valid for x ≈ 0). Write your code such that the output looks like the slide. with the first column with the function, in the second one the rationalized one and in the third the Taylor expansion.
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

double fx1r(const double&);
double fx1c(const double&);
double fx2r(const double&);
double fx2c(const double&);
double frTaylor(const double&);
double fcTaylor(const double&);

int main(){
    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    double xr, xc; // Definisco le x da usare per le funzioni

    cout << "+----------------------------------------+" << endl;
    cout << "Esempio funzione: sqrt(x^2 + 1) - x per grandi x" << endl;
    for(int i=3; i<10; i++){
        xr = pow(10,i);
        cout << "x = " << xr << "; fx1 = " << fx1r(xr) << "; fx2 = " << fx2r(xr) << "; f(Taylor) = " << frTaylor(xr) << endl;
    }

    cout << endl;
    cout << endl;
    cout << "Esempio funzione: 1 - cos(x) per piccoli x" << endl;
    for(int i=1; i<8; i++){
        xc = pow(10,-i);
        cout << "x = " << xc << "; fx1 = " << fx1c(xc) << "; fx2 = " << fx2c(xc) << "; f(Taylor) = " << fcTaylor(xc) << endl;
    }
    cout << "+----------------------------------------+" << endl;

    return 0;
}




// Definizione funzione sqrt
double fx1r(const double& x){
    double a = sqrt(x*x + 1.) - x;
    return a;
}

// Definizione funzione cos
double fx1c(const double& x){
    double b = 1. - cos(x);
    return b;
}

// Definizione funzione sqrt razionalizzata
// moltiplico sopra e sotto per la funzione cambiata di segno
double fx2r(const double& x){
    double c = 1./(sqrt(x*x + 1.) + x);
    return c;
}

// Definizione funzione cos razionalizzata
// moltiplico sopra e sotto per la funzione cambiata di segno
double fx2c(const double& x){
    double d = sin(x)/(1.+cos(x));
    return d;
}

// Definizione sviluppo di Taylor della sqrt
double frTaylor(const double& x){
    //double c = x + 1./(2.*x) - 1./(8.*x*x*x) - x;
    double c = 1./(2.*x);
    return c;
}

// Definizione sviluppo di Taylor del cos
double fcTaylor(const double& x){
    double d;
    //d = 1. - (1. - (x*x)/2. + (x*x*x*x)/(24.));
    d = (x*x)/2.;
    return d;
}