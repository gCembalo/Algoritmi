// Esercizio analogo a "quadratic.cpp" ma con maggior dettaglio
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

void sol(double&,double&);
void ordsol(double&,double&);
void norm(const double&,const double&,const double&,double&,double&);
void alt(const double&,const double&,const double&,double&,double&);
void signb(const double&,const double&,const double&,double&,double&);

int main(){
    cout << setprecision(5);
    cout << setiosflags ( ios::scientific );

    // Prendo le soluzioni
    double x1, x2;
    sol(x1, x2);
    ordsol(x1, x2);

    cout << "+--------------------------------------------------------+" << endl;

    // Calcolo i coefficienti
    double a = 1.0;
    double b = -(x1 + x2);
    double c = x1*x2;
    cout << "I coefficienti sono:\na = " << a << "\nb = " << b << "\nc = " << c << endl; 

    // Calcolo il delta
    double delta = b*b - 4*a*c;
    if(delta < 0){
        cout << "Non esistono soluzioni." << endl;
        return 0;
    }

    // Formula standard
    double y1, y2;      // con 1 sol + e 2 sol -
    norm(a, b, c, y1, y2);
    ordsol(y1, y2);
    cout << "+--------------------------------------------------------+" << endl;
    cout << "Con il metodo tradizionale:" << endl;
    cout << "La soluzione x1 = " << y1 << " err = " << fabs(y1/x1 - 1.0) << endl;
    cout << "La soluzione x2 = " << y2 << " err = " << fabs(y2/x2 - 1.0) << endl;

    // Formula alternativa
    double y3, y4;      // con 3 sol + e 4 sol -
    alt(a, b, c, y3, y4);
    ordsol(y3, y4);
    cout << "+--------------------------------------------------------+" << endl;
    cout << "Con il metodo alternativo (fattorizzazione):" << endl;
    cout << "La soluzione x1 = " << y3 << " err = " << fabs(y3/x1 - 1.0) << endl;
    cout << "La soluzione x2 = " << y4 << " err = " << fabs(y4/x2 - 1.0) << endl;
    
    // Implemento la selezione della formula a seconda di b
    double y5, y6;
    signb(a, b, c, y5, y6);
    ordsol(y5, y6);
    cout << "+--------------------------------------------------------+" << endl;
    cout << "Si vede che a seconda del segno di b bisogna utilizzare una forma o l'altra. Le soluzioni sono:" << endl;
    cout << y5 << " err = " << fabs(y5/x1 - 1) << endl;
    cout << y6 << " err = " << fabs(y6/x2 - 1) << endl;
    cout << "+--------------------------------------------------------+" << endl;

    return 0;
}


void sol(double& x1, double& x2){
    // Prendo i coefficienti
    cout << "Abbiamo l'equazione:" << endl;
    cout << "    ax^2 + bx + c = 0      " << endl;
    cout << "Dammi la soluzione 1: ";
    cin >> x1;
    cout << "Dammi la soluzione 2: ";
    cin >> x2;
}

void ordsol(double& x1, double& x2){
    // Ordino le soluzioni
    double xtemp;
    if(x2 >= x1){
        xtemp = x1;
        x1 = x2;
        x2 = xtemp;
    }
}

void norm(const double& x, const double& y, const double& z, double& s1, double& s2){
    s1 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);
    s2 = ( -y - sqrt(y*y - 4*x*z) )/(2*x);
    //cout << "La soluzione con (+) con il metodo standard: " << s1 << endl;
    //cout << "La soluzione con (-) con il metodo standard: " << s2 << endl;
}

void alt(const double& x, const double& y, const double& z, double& s3, double& s4){
    s3 = -(2*z)/( y + sqrt(y*y - 4*x*z) );
    s4 = -(2*z)/( y - sqrt(y*y - 4*x*z) );
    //cout << "La soluzione con (+) con il metodo alternativa: " << s3 << endl;
    //cout << "La soluzione con (-) con il metodo alternativa: " << s4 << endl;

}

void signb(const double& x, const double& y, const double& z, double& s1, double& s2){
    if(y>=0){
        s1  = ( -y - sqrt(y*y - 4*x*z) )/(2*x);
        s2 = -(2*z)/( y + sqrt(y*y - 4*x*z) );
    } else if (y<0){
        s1 = -(2*z)/( y - sqrt(y*y - 4*x*z) );
        s2 = ( -y + sqrt(y*y - 4*x*z) )/(2*x);
    }
}