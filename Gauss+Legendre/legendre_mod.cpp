// write a code that can found the zero of Legendre polynomial of arbitrary order n. In order to evaluate Pn(x) use Bonnet’s recursion formula written in the slide.
//
// with P0(x) = 1, P1 (x)= x. For Newton-Raphson method, the derivative may be computed in terms of Pn(x) and Pn-1(x) through (the formula in the slide).
//
// Now that you have found the roots with high precision, you may use them to produce your Gaussian weights (written in the slide).
//

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

    // definisco il valore del polinomio e della derivata che uso solo per stamparlo nel file
    double pol, derPol;

    for( double i = -1.0 ; i < 1.0 ; i += 1.e-3 ){

        Legendre(i, pol, derPol);
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
    cout << "(Iterazioni)           Zeri            Pesi\n" << endl;

    cout << "Bisezione:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        bisection(polLegendre, xL[k], xR[k], tol, x0[k], l);
        cout << "(" << l << ")      " << "x" << k << " = " << x0[k] << "        w" << k << " = " << wi(x0[k]) << endl;
    }

    cout << "\nFalseposition:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        false_position(polLegendre, xL[k], xR[k], tol, x1[k], l);
        cout << "(" << l << ")      " << "x" << k << " = " << x1[k] << "        w" << k << " = " << wi(x1[k]) << endl;
    }

    cout << "\nSecante:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        secant_method(polLegendre, xL[k], xR[k], tol, x2[k], l);
        cout << "(" << l << ")      " << "x" << k << " = " << x2[k] << "        w" << k << " = " << wi(x2[k]) << endl;
    }

    cout << "\nNewton:\n" << endl;
    for( int k = 0 ; k < nroots ; k++ ){
        newton_method(polLegendre, derpolLegendre, xL[k], xR[k], tol, 1.e-8, x3[k], l);
        cout << "(" << l << ")      " << "x" << k << " = " << x3[k] << "        w" << k << " = " << wi(x3[k]) << endl;
    }

    cout << "\n+-----------------------------------------------------------------------------+" << endl;
    /*cout << "I pesi sono:\n" << endl;

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
    cout << "\n+-----------------------------------------------------------------------------+" << endl;*/

    return 0;
}


// funzione unica per polinomi, derivata e pesi di Legendre. Prendo in input la x, un riferimento al polinomio e uno alla derivata; ho una variabile globale g_LegendreN definita
void Legendre(double x, double &pn, double &dpn){

    // mi faccio restituire i polinomi di grado 0 o 1
    if (g_LegendreN == 0){

        pn = 1.0;
        dpn = 0.0;
        
    }

    if (g_LegendreN == 1){

        pn = x;
        dpn = 1.0;
        
    }

    // definisco i polinomi
    double P0 = 1., Pi = x, Pi1;  // P0 è P_0 ; Pi è P_n ; Pi1 è P_{n+1}


    for( int i = 1 ; i < g_LegendreN ; i++ ){

        Pi1 = ( (2.0 * i + 1.0) * x * Pi - i * P0 ) / ( i+1.0 );

        // Aggiorno i valori
        P0 = Pi;
        Pi = Pi1;
       
    }

    // ora Pi1 è il polinomio di grado g_LegendreN mentre P0 è il polinomio di grado g_LegendreN-1
    pn = Pi1;

    // calcolo la derivata del polinomio di grado g_LegendreN usando la formula nota
    double dPi = ( g_LegendreN * (x * Pi1 - P0) ) / ( x * x - 1.0 );
    dpn = dPi;

}


double polLegendre(double x){

    if (g_LegendreN == 0){
        return 1.;
    }

    if (g_LegendreN == 1){
        return x;
    }

    // definisco i polinomi 0 e 1
    double P0 = 1., Pi = x, Pi1;  // P0 è P_0 ; Pi è P_n ; Pi1 è P_{n+1}

    for( int i = 1 ; i < g_LegendreN ; i++ ){

        Pi1 = ( (2.0 * i + 1.0) * x * Pi - i * P0 ) / ( i+1.0 );

        // Aggiorno i valori
        P0 = Pi;
        Pi = Pi1;
       
    }

    return Pi1;

}


double derpolLegendre(double x){

    if (g_LegendreN == 0.0){
        return 0.0;
    }

    if (g_LegendreN == 1.0){
        return 1.0;
    }

    double P0 = 1.0; // P0(x)
    double P1 = x; // P1(x)
    double Pi1; // P_{n+1}

    // stesso ciclo della funzione del polinomio
    for( int i = 1 ; i < g_LegendreN ; i++ ){

        Pi1 = ( (2.0 * i + 1.0) * x * P1 - i * P0 ) / ( i + 1.0 );

        // aggiorno i valori
        P0 = P1;
        P1 = Pi1;

    }

    // uso la formula per la derivata
    double dPi = ( g_LegendreN * (x * P1 - P0) ) / ( x * x - 1.0 );

    return dPi;
    
}

double wi(double x){

    return 2.0/( ( 1.0-x*x )*( derpolLegendre(x)*derpolLegendre(x) ) );

}