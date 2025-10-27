// given an interval [a,b] to compute the integral of f(x). Divide the interval [a,b] into N equally spaced sub-interval separated by N+1 points {x0, x1, x2, … xN}. In each sub-interval apply the rectangular, trapezoidal and Simpson rules.
//
// Consider f(x) = exp(-x) and use [a, b] = [0,1].
//
// Next, iterate by doubling the value of intervals N ( = 4, 8, 16, …) until convergence is achieved. tol is a prescribed tolerance (e.g. 10-5).
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// Funzione
double func(double);
// Diversi metodi di quadratura
double RectangularRule(double (*)(double),double,double,int);
double MidPointRule(double (*)(double),double,double,int);
double TrapezoidalRule(double (*)(double),double,double,int);
double ExtSimpsonRule(double (*)(double),double,double,int);
// Funzioni per la convergenza
double ConvergenceRectangular(double (*)(double),double,double,double);
double ConvergenceMidPoint(double (*)(double),double,double,double);
double ConvergenceTrapezoidal(double (*)(double),double,double,double);
double ConvergenceSimpson(double (*)(double),double,double,double);

int main(){
    cout << setprecision(11);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    // Definisco gli estremi dell'intervallo e il numero di intervalli
    double a = 0, b = 1;
    int N = 4;

    // Con un numero N fissato
    cout << "\n+------------------------------------------+" << endl;
    cout << "Con N = " << N << "\n" << endl;
    // Rectangular Rule
    double rec = RectangularRule(func, a, b, N);
    // Mid Point Rule
    double mid = MidPointRule(func, a, b, N);
    // Trapezoidal Rule
    double trap = TrapezoidalRule(func, a, b, N);
    // Simpson Rule
    double simp = ExtSimpsonRule(func, a, b, N);

    cout << "Rectangular: " << rec << ";\nMid Point: " << mid << ";\nTrapezoidal: " << trap << "\nSimpson: " << simp << ";" << endl;
    cout << "+------------------------------------------+" << endl;
    cout << "Convergenza: \n" << endl;

    // Controllo quanti step prima della convergenza
    double tol = 1e-5;
    ConvergenceRectangular(func,a,b,tol);
    ConvergenceMidPoint(func,a,b,tol);
    ConvergenceTrapezoidal(func,a,b,tol);
    ConvergenceSimpson(func,a,b,tol);

    return 0;
}





// Funzione che vogliamo integrare
double func(double x){
    return exp(-x);
}

// Regola quadratura
double RectangularRule(double (*F)(double), double a, double b, int N){
    // Definisco l'ampiezza di ogni intervallo
    double h = fabs(a-b)/(double)N;
    // Definisco la variabile che somma i rettangoli
    double sum = 0.0;
    // Calcolo dell'integrale
    for(int i=0; i<=N-1; i++){
        sum += F(a + i*h) * h; // Integrale
    }
    return sum;
}

double MidPointRule(double (*F)(double), double a, double b, int N){
    // Definisco l'ampiezza di ogni intervallo
    double h = fabs(a-b)/(double)N;
    // Definisco la variabile che somma i rettangoli
    double sum = 0.0;
    // Calcolo dell'integrale
    for(int i=0; i<=N-1; i++){
        sum += F(a + i*h + h/2) * h;
    }
    return sum;
}

double TrapezoidalRule(double (*F)(double), double a, double b, int N){
    // Definisco l'ampiezza di ogni intervallo
    double h = fabs(a-b)/(double)N;
    // Definisco la variabile che somma i rettangoli
    double sum = 0.0;
    // Calcolo dell'integrale
    for(int i=0; i<=N-1; i++){
        sum += 0.5 * F(a + i*h) * h + 0.5 * F(a + (i+1)*h) * h;
    }
    return sum;
}

double ExtSimpsonRule(double (*F)(double), double a, double b, int N){
    int w = 4; // Definizione del peso
    double h, sum;
    // Definisco l'ampiezza di ogni intervallo
    h = fabs(a-b)/(double)N;
    sum = (F(a)+F(b))*h/3;
    for (int i=1; i<N; i++ )
    {
        sum += w*F(a + i*h)*h/3;
        w = 6-w; // Calcolo il peso ogni volta
    }
    return sum;
}

double ConvergenceRectangular(double (*F)(double), double a, double b, double tol){
    // Controllo quanti intervalli servono
    int i = 2;
    while ( fabs( RectangularRule(func, a, b, i) - RectangularRule(func, a, b, i/2) ) > tol ){
        i = 2*i;
    }
    cout << "Rectangular: " << RectangularRule(func, a, b, i) << " iter: " << i << endl;
    return 0.;
}

double ConvergenceMidPoint(double (*F)(double), double a, double b, double tol){
    // Controllo quanti intervalli servono
    int i = 2;
    while ( fabs( MidPointRule(func, a, b, i) - MidPointRule(func, a, b, i/2) ) > tol ){
        i = 2*i;
    }
    cout << "Mid Point: " << MidPointRule(func, a, b, i) << " iter: " << i << endl;
    return 0.;
}

double ConvergenceTrapezoidal(double (*F)(double), double a, double b, double tol){
    // Controllo quanti intervalli servono
    int i = 2;
    while ( fabs( TrapezoidalRule(func, a, b, i) - TrapezoidalRule(func, a, b, i/2) ) > tol ){
        i = 2*i;
    }
    cout << "Trapezoidal: " << TrapezoidalRule(func, a, b, i) << " iter: " << i << endl;
    return 0.;
}

double ConvergenceSimpson(double (*F)(double), double a, double b, double tol){
    // Controllo quanti intervalli servono
    int i = 2;
    while ( fabs( ExtSimpsonRule(func, a, b, i) - ExtSimpsonRule(func, a, b, i/2) ) > tol ){
        i = 2*i;
    }
    cout << "Simpson: " << ExtSimpsonRule(func, a, b, i) << " iter: " << i << endl;
    return 0.;
}