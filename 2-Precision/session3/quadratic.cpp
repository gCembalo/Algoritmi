// using double precision arithmetic, write a computer program to solver the quadratic quadratic equation. Using, at first, the standard formula. In order to avoid catastrophic cancellation, implement the selective expressions depending on the sign of the b coefficient.
//
// Test your solver on the following cases:
// a = 1    b = -(x1+x2)    c = x1*x2
//
// with:
// x1 = 2   x2 = -3
// x1 = 1.e-5    x2 = 1.e8
// x1 = 1.e-12      x2 = 1.e12
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

void dati(double&,double&,double&);
void norm(const double&,const double&,const double&,double&,double&);
void alt(const double&,const double&,const double&,double&,double&);
void signb(const double&,const double&,const double&,double&,double&);

int main(){
    
    // Prendo i dati
    double a, b, c;
    dati(a, b, c);

    cout << "+--------------------------------------------------------+" << endl;
    cout << setprecision(10);

    // Formula standard
    double x1, x2;      // con 1 sol + e 2 sol -
    norm(a, b, c, x1, x2);

    // Formula alternativa
    double x3, x4;      // con 3 sol + e 4 sol -
    alt(a, b, c, x3, x4);
    cout << "+--------------------------------------------------------+" << endl;

    cout << "La differenza tra le soluzioni positive è: " << x1-x3 << endl;
    cout << "La differenza tra le soluzioni negative è: " << x2-x4 << endl;
    
    // Implemento la selezione della formula a seconda di b
    double y1, y2;
    signb(a, b, c, y1, y2);

    return 0;
}


void dati(double& x, double& y, double& z){
    // Prendo i coefficienti
    cout << "Abbiamo l'equazione:" << endl;
    cout << "    ax^2 + bx + c = 0      " << endl;
    cout << "Dammi a: ";
    cin >> x;
    cout << "Dammi b: ";
    cin >> y;
    cout << "Dammi c: ";
    cin >> z;
}

void norm(const double& x, const double& y, const double& z, double& s1, double& s2){
    s1 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);
    s2 = ( -y - sqrt(y*y - 4*x*z) )/(2*x);
    cout << "La soluzione con (+) con il metodo standard: " << s1 << endl;
    cout << "La soluzione con (-) con il metodo standard: " << s2 << endl;
}

void alt(const double& x, const double& y, const double& z, double& s3, double& s4){
    s3 = -(2*z)/( y + sqrt(y*y - 4*x*z) );
    s4 = -(2*z)/( y - sqrt(y*y - 4*x*z) );
    cout << "La soluzione con (+) con il metodo alternativa: " << s3 << endl;
    cout << "La soluzione con (-) con il metodo alternativa: " << s4 << endl;

}

void signb(const double& x, const double& y, const double& z, double& s1, double& s2){
    if(y>=0){
        s1  = ( -y - sqrt(y*y - 4*x*z) )/(2*x);
        s2 = -(2*z)/( y + sqrt(y*y - 4*x*z) );
    } else if (y<0){
        s1 = -(2*z)/( y - sqrt(y*y - 4*x*z) );
        s2 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);
    }

    cout << "+--------------------------------------------------------+" << endl;
    cout << "Si vede che a seconda del segno di b bisogna utilizzare una forma o l'altra. Le soluzioni sono:" << endl;
    cout << s1 << endl;
    cout << s2 << endl;
    cout << "+--------------------------------------------------------+" << endl;
}