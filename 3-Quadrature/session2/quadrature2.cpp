// Using Simpson rule with 2 intervals (3 points in total) and Gauss-Legendre (1 interval, 3 Gaussian points), compute the integral:
// int_{0}^{3} \sqrt{ 1+t } dt = 4.66666667
//
// Now, using the same code, test the two algorithms in computing the integral:
// int_{-1}^{5} f(x) dt
// with f(x) = 1 - x + 2x^2 + x^3 / 2 + x^4 / 4 - x^5 / 8
// again using 2 intervals for Simpson’s rule and 1 interval for Gauss-Legendre with 3 Gaussian points. The exact result is –66/5. Which of the two is better ?
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// Funzione
double func1(double);
double func2(double);
// Diversi metodi
double Simpson(double (*)(double),double,double,int);
double Gauss(double (*)(double),double,double,int,int);

int main(){
    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    // Definisco il numero di intervalli
    int Ns, Ng;
    // Definisco il numero di punti
    int ng;

    // Definisco gli estremi di integrazione per i due integrali
    double a = 0., b = 3., c = -1., d = 5.;

    cout << "\n+-----------------------------------+" << endl;
    cout << "Funzione: sqrt(1+x)\n" << endl;
    //cout << "Gli esercizi sono da fare con 2 intervalli per Simpson\n e 1 intervallo e 3 punti per Gauss.\n" << endl;
    //cout << "Dammi il numero di intervalli per il metodo di Simpson: ";
    //cin >> Ns;
    //cout << "Dammi il numero di intervalli per il metodo di Gauss: ";
    //cin >> Ng;
    //cout << "Dammi il numero di punti per il metodo di Gauss: ";
    //cin >> ng;
    Ns = 2 , Ng = 1 , ng = 3;
    double simp = Simpson(func1,a,b,Ns);
    double gauss = Gauss(func1,a,b,Ng,ng);
    cout << "Valore esatto: 4.666666677" << endl;
    cout << "Simpson:       " << simp << endl;
    cout << "Gauss:         " << gauss << endl;

    cout << "\n+-----------------------------------+" << endl;
    double Ns2, Ng2, ng2;
    cout << "Funzione: 1 - x + 2x^2 + x^3/2 + x^4/4 - x^5/8 \n" << endl;
    //cout << "Gli esercizi sono da fare con 2 intervalli per Simpson\n e 1 intervallo e 3 punti per Gauss.\n" << endl;
    //cout << "Dammi il numero di intervalli per il metodo di Simpson: ";
    //cin >> Ns2;
    //cout << "Dammi il numero di intervalli per il metodo di Gauss: ";
    //cin >> Ng2;
    //cout << "Dammi il numero di punti per il metodo di Gauss: ";
    //cin >> ng2;
    Ns2 = 2 , Ng2 = 1 , ng2 = 3;
    double simp2 = Simpson(func2,c,d,Ns2);
    double gauss2 = Gauss(func2,c,d,Ng2,ng2);
    cout << "Valore esatto: " << -66./5. << endl;
    cout << "Simpson: " << simp2 << endl;
    cout << "Gauss: " << gauss2 << endl;

    return 0;
}



// Funzioni
double func1(double x){
    return sqrt(1.+x);
}

double func2(double x){
    return 1. - x + 2.*x*x + 0.5*x*x*x + x*x*x*x/4. - x*x*x*x*x/8.;
}

double Simpson(double (*F)(double), double a, double b, int N){
    int w = 4; // Definizione del peso
    double h, sum;
    // Definisco l'ampiezza di ogni intervallo
    h = fabs(a-b)/(double)N;
    sum = (F(a)+F(b))*h/3.;
    for (int i=1; i<N; i++ )
    {
        sum += w*F(a + i*h)*h/3.;
        w = 6.-w; // Calcolo il peso ogni volta
    }
    return sum;
}

double Gauss(double (*F)(double), double a, double b, int N, int ng){
    // Implemento la funzione per funzionare fino a 5 punti
    double w[ng], x[ng]; // array di pesi e zeri del polinomio di Legendre
    double sum = 0.0;

    // Riempio gli array in base al numero di punti che abbiamo
    if(ng == 1){
        x[0] = 0.;
        w[0] = 2.;
    }
    else if(ng == 2){
        x[1] = sqrt(1./3.);
        x[0] = -sqrt(1./3.);
        w[0] = 1.;
        w[1] = 1.;
    }
    else if(ng == 3){
        x[0] = -sqrt(3./5.);
        x[1] = 0.;
        x[2] = sqrt(3./5.);
        w[0] = 5./9.;
        w[1] = 8./9.;
        w[2] = 5./9.;
    }
    else if(ng == 4){
        x[0] = -sqrt( 3./7. - 2./7.*sqrt(6./5.) );
        x[1] = -sqrt( 3./7. + 2./7.*sqrt(6./5.) );
        x[2] = sqrt( 3./7. - 2./7.*sqrt(6./5.) );
        x[3] = sqrt( 3./7. + 2./7.*sqrt(6./5.) );
        w[0] = ( 18. + sqrt(30.) )/( 36. );
        w[1] = ( 18. - sqrt(30.) )/( 36. );
        w[2] = ( 18. + sqrt(30.) )/( 36. );
        w[3] = ( 18. - sqrt(30.) )/( 36. );
    }
    else if(ng == 5){
        x[0] = -1./3.*sqrt( 5. - 2.*sqrt(10./7.) );
        x[1] = -1./3.*sqrt( 5. + 2.*sqrt(10./7.) );
        x[2] = 0.;
        x[3] = 1./3.*sqrt( 5. - 2.*sqrt(10./7.) );
        x[4] = 1./3.*sqrt( 5. + 2.*sqrt(10./7.) );

        w[0] = ( 322. + 13.*sqrt(70.) )/( 900. );
        w[1] = ( 322. - 13.*sqrt(70.) )/( 900. ); 
        w[2] = 128./225.;
        w[3] = ( 322. + 13.*sqrt(70.) )/( 900. );
        w[4] = ( 322. - 13.*sqrt(70.) )/( 900. );
    }
    else{
        cout << "codice non implementato per Ng>5." << endl;
        return 0.;
    }

    // Calcolo l'integrale
    double h = fabs(b-a)/(double)N;
    double sumj = 0.0;

    for(int i=0; i<N; i++){
            // Devo fare il cambio di variabili che poi viene iterato sui vari intervalli successivi
            double x0 = a + i*h;
            double x1 = x0 + h;
            double s1 = (x1-x0)/2;
            double s2 = (x1+x0)/2;
            for(int j=0; j<ng; j++){
                sumj += s1 * w[j] * F( s1*x[j] + s2 );
            }
            sum += sumj;
        } 
    
    // Restituisco il risultato
    return sum;
}