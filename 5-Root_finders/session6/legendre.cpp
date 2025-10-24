#include "root_finders.h"

double polLegendre(double);
double derpolLegendre(double);
double wi(double);

int g_LegendreN;

int main(){

    // definisco gli estremi dell'intervallo
    double a = -1, b = 1;
    // defiisco una tolleranza
    double tol = 1.e-8;
    //double facc = 1e-6;

    cout << "\n+-----------------------------------------------------------------------------+\nVogliamo studiare i polinomi di Legendre, che grado vuoi?" << endl;
    cin >> g_LegendreN;

    // definisco il numero di roots che avremo
    int nroots = 0;
    // definisco gli arrai degli estremi
    double xL[g_LegendreN], xR[g_LegendreN];
    // definisco un array dove mettere gli zeri
    double x0[g_LegendreN], x1[g_LegendreN], x2[g_LegendreN], x3[g_LegendreN];

    cout << setprecision(12);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation


    // stampo nel file il polinomio
    ofstream fdata; // declare Output stream class to operate on files
    fdata.open("Legendre.dat"); // open output file
    //set the significant digits and the scientific notation
    fdata << setprecision(7);
    fdata << setiosflags ( ios::scientific );

    for( double i = -1.0 ; i < 1.0 ; i += 1.e-3 ){

        double pol = polLegendre(i);
        fdata << i << "     " << pol << endl;

    }
    fdata.close();


    // invoco il bracketing
    Bracket(polLegendre,a,b,2*g_LegendreN+1,xL,xR,nroots);

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
    cout << "Troviamo gli zeri:\n" << endl;

    cout << "Bisezione:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        bisection(polLegendre, xL[k], xR[k], tol, x0[k], l);
        cout << "x" << k << " = " << x0[k] << "     trovata con " << l << " iterazioni" << endl;
    }

    cout << "\nFalseposition:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        false_position(polLegendre, xL[k], xR[k], tol, x1[k], l);
        cout << "x" << k << " = " << x1[k] << "     trovata con " << l << " iterazioni" << endl;
    }

    cout << "\nSecante:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        secant_method(polLegendre, xL[k], xR[k], tol, x2[k], l);
        cout << "x" << k << " = " << x2[k] << "     trovata con " << l << " iterazioni" << endl;
    }

    cout << "\nNewton:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        newton_method(polLegendre, derpolLegendre, xR[k], tol, x3[k], l);
        cout << "x" << k << " = " << x3[k] << "     trovata con " << l << " iterazioni" << endl;
    }

    cout << "\n+-----------------------------------------------------------------------------+" << endl;
    cout << "Troviamo gli zeri:\n" << endl;

    cout << "Bisezione:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        //bisection(polLegendre, xL[k], xR[k], tol, x0[k], l);
        cout << "w" << k << " = " << wi(x0[k]) << endl;
    }

    cout << "\nFalseposition:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        //false_position(polLegendre, xL[k], xR[k], tol, x1[k], l);
        cout << "w" << k << " = " << wi(x1[k]) << endl;
    }

    cout << "\nSecante:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        //secant_method(polLegendre, xL[k], xR[k], tol, x2[k], l);
        cout << "w" << k << " = " << wi(x2[k]) << endl;
    }

    cout << "\nNewton:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        //newton_method(polLegendre, derpolLegendre, xR[k], tol, x3[k], l);
        cout << "w" << k << " = " << wi(x3[k]) << endl;
    }
    cout << "\n+-----------------------------------------------------------------------------+" << endl;

    return 0;
}


double polLegendre(double x){

if (g_LegendreN == 0) return 1.;
    if (g_LegendreN == 1) return x;

    double Pi_1 = 1, Pi = x, Pi1;

    for (int i=1; i<g_LegendreN; i++){

        Pi1 = ((2 * i + 1) * x * Pi - i* Pi_1) / (i+1.);

        // Aggiorno i valori
        Pi_1 = Pi;
        Pi = Pi1;
       
    }

    return Pi1;

}

double derpolLegendre(double x){

if (g_LegendreN == 0) return 0.0;
    if (g_LegendreN == 1) return 1.0;

    double Pnm1 = 1.0;   // P0(x)
    double Pn = x;       // P1(x)
    double Pnp1;

    // Calcolo Pn e Pn-1
    for (int i = 1; i < g_LegendreN; i++) {
        Pnp1 = ((2.0 * i + 1.0) * x * Pn - i * Pnm1) / (i + 1.0);
        Pnm1 = Pn;
        Pn = Pnp1;
    }

    // Uso la formula per la derivata
    double dPn = (g_LegendreN * (x * Pn - Pnm1)) / (x * x - 1.0);

    return dPn;
    
}

double wi(double x){

    return 2.0/( (1.0-x*x)*( derpolLegendre(x)*derpolLegendre(x) ) );

}