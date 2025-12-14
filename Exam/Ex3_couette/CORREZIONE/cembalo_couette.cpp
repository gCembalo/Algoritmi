// Last name: Cembalo
// u'(r_1)               = -8.9899007e-01
// L1err(Shooting)       = 2.9317587e-08
// L1err(Finite Diff)    = 3.5447246e-05
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#define NMAX_EQ 64 // numero massimo di eq (sicurezza)

using namespace std;

double ExSolution(double);
void RHSFunc(double, double *, double *);
double Residual(double);
void RK4Step(double, double *, void (*)(double, double *, double *), double, int);
int secant_method(double (*F)(double), double, double, double, double &, int &);
void TridiagSolver(double *, double *, double *, double *, double *, int);

int main(){

    ofstream fdata;
    fdata.open("couette.dat"); // file per le soluzioni

    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 7 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 10 );

    int neq = 2; // numero equazioni
    int nstep = 200; // numero punti
    double tol = 1.e-8;

    // errori L1 normati per entrambi i metodi
    double errL1RK4 , errL1FD;

    double r; // variabile raggio
    double r1 = 1.0;
    double r2 = 10.0;
    double dr = fabs(r2 - r1) / nstep; // step fisso

    // condizioni iniziali velocità angolare
    double u1 = 1.0;
    double u2 = 0.6;

    double res; // variabile residuo
    double szero; // variabile per lo zero del residuo
    int l = 0; // variaible per le iterazioni 
               // (inutile per l'esercizio)

    // variabile derivata di u
    double s; // guess, derivata di u in r
    double s0 = -1.0 , sf = -5.0;
    double ds = fabs( s0 - sf )/(double)nstep;

    // plotto il residuo ( non utilizzato )
    for( int i = 0 ; i < nstep ; i++ ){

        s = s0 + i*ds;
        res = Residual(s);

        fdata << s << "  " << res << endl;

    }

    fdata << endl << endl;


    // uso secante per trovare lo zero del residuo e quindi il valore di s
    double status;
    status = secant_method(Residual, s0, sf, tol, szero, l);

    // aggiungo un controllo sullo zero del residuo
    if (status == 0){
        cout << "\nu'(r_1) = " << szero << endl;
    }else{
        cout << "No solution!" << endl;
    }

    // definisco l'array delle soluzioni e imposto le condizioni iniziali
    double Y[neq];
    Y[0] = u1;
    Y[1] = szero;

    r = r1;

    // plotto la condizione iniziale
    fdata << r << "  " << Y[0] << "  " << Y[1] << endl;

    // risolvo le equazioni del moto usando RK
    for( int i = 0 ; i < nstep ; i++ ){

        RK4Step(r, Y, RHSFunc, dr, neq); // risolvo la ODE
        r += dr;

        // stampo nel file i dati
        fdata << r << "  " << Y[0] << "  " << Y[1] << endl;

        // calcolo l'errore L1 normato
        errL1RK4 += fabs( ExSolution(r) - Y[0] );

    }

    // devo ancora dividere per il numero di punti l'errore
    errL1RK4 /= nstep;

    fdata << endl << endl;


    // Implemento il finite difference method
    int n = nstep + 1; // punti griglia

    // definisco i vettori degli elementi sopra e sotto la diagonale
    double *a, *b, *c, *R, *y;

    // ho chiamato il right-hand side vector del metodo tridiagonal
    // con la lettera maiuscola per differenziarlo dalla variabile raggio
    a = new double [n];
    b = new double [n];
    c = new double [n];
    R = new double [n];
    y = new double [n];

    double alpha , beta;
    alpha = u1;
    beta = u2;

    double h; // incremento
    h = fabs( r2 - r1 )/( (double)( n - 1 ) );

    // condizioni al bordo
    y[0] = alpha;
    y[n-1] = beta;

    // riempio i vettori a, b, c ed r
    for( int i = 0 ; i < n ; i++ ){

        r = r1 + i*h;

        a[i] = 1.0 - 0.5*h/r; // sotto la diagonale (y_{i-1})
        b[i] = - 2.0 - h*h/( r*r ); // diagonale (y_{i})
        c[i] = 1.0 + 0.5*h/r; // sopra la diagonale (y_{i+1})

        R[i] = 0.0;

    }

    // metto le condizioni iniziali nel primo e ultimo elemento di r
    R[1] -= a[1]*y[0];
    R[n-1] -= c[n-1]*y[n-1];

    // richiamo il TridiagSolver per trovare il vettore delle soluzioni
    // schiftiamo di 1 l'argomento, poiche' l'indice 0 lo fissiamo con la
    // condizione al bordo e quindi il primo elemento sara' y2
    TridiagSolver(a+1, b+1, y+1, R+1, c+1, n-2);

    // stampo il vettore nel file
    for ( int i = 0 ; i < n ; i++ ){

        r = r1 + i*h;

        fdata << r << " " << y[i] << endl;

        // calcolo l'errore L1 normato
        errL1FD += fabs( ExSolution(r) - y[i] );

    }

    // devo ancora dividere per il numero di punti l'errore
    errL1FD /= n;

    fdata.close();

    cout << "L1err(Shooting)       = " << errL1RK4 << endl;
    cout << "L1err(Finite Diff)    = " << errL1FD << endl;
    cout << "\nSoluzione salvata in couette.dat\n" << endl;

    // pulisco
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] R;
    delete[] y;

    return 0;

}


// ---------------------- Funzioni ------------------------ //

// Implemento una funzione che calcola la funzione analitica in un punto
double ExSolution(double r){

    // condizioni problema
    double r1 = 1.0;
    double r2 = 10.0;
    double u1 = 1.0;
    double u2 = 0.6;

    // calcolo i parametri della soluzione analitica
    double a = ( r1*r2*( r1*u2 - r2*u1 ) )/( r1*r1 - r2*r2 );
    double b = ( u1*r1 - u2*r2 )/( r1*r1 - r2*r2 );

    double uEx;

    uEx = a/r + b*r;

    return uEx;

}

// definisco il Right-Hand-Side-Function (e' problem dependent). 
// Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore)
// mi faccio dare R
void RHSFunc(double x, double *Y, double *R){

    R[0] = Y[1];
    R[1] = Y[0]/( x*x ) - Y[1]/x;

}

// creo la funzione residuo
// sara' la funzione che diamo alla funzione per la ricerca degli zeri
// e' problem dependent.
// gli do in input la guess sulla derivata
double Residual(double s){

    int neq = 2; // numero equazioni
    int nstep = 200; // numero punti
    double tol = 1.e-8;

    double r; // variabile raggio
    double r1 = 1.0;
    double r2 = 10.0;
    double dr = fabs(r2 - r1) / nstep; // step fisso

    // condizioni iniziali velocità angolare
    double u1 = 1.0;
    double u2 = 0.6;

    // definisco l'array delle soluzioni e imposto le condizioni iniziali
    double Y[neq];
    Y[0] = u1;
    Y[1] = s;

    r = r1; // setto il punto di integrazione iniziale

    // risolvo le equazioni del moto usando RK a 4 step
    for( int i = 0 ; i < nstep ; i++ ){

        RK4Step(r, Y, RHSFunc, dr, neq); // risolvo la ODE
        r += dr;

    }

    // calcolo e ritorno il residuo (valore atteso e' u2 = 0.6)
    return Y[0] - u2;

}

// implemento il metodo Runge-Kutta del quarto ordine.
// gli do in input la variabile di integrazione, il puntatore alle soluzioni,
// il puntatore alla funzione del Right-Hand-Side-Function, l'incremento e
// l'ordine della ODE.
void RK4Step(double t, double *Y, void (*RHSFunc)(double t, double *Y, double *R),
             double h, int neq){
    
    // definisco i vettori per gli step intermedi
    double Y1[NMAX_EQ], k1[NMAX_EQ], k2[NMAX_EQ], k3[NMAX_EQ], k4[NMAX_EQ];
    
    RHSFunc(t,Y,k1); // calcolo k1 con il RSH con t_n e Y_n

    // scrivo il ciclo per determinare Y_n + k1*h/2
    for( int i = 0 ; i < neq ; i++ ){
        
        Y1[i] = Y[i] + 0.5*h*k1[i];

    }

    RHSFunc(t+0.5*h,Y1,k2); // calcolo k2 con il RSH con t_n+h/2 e Y_n+k1*h/2
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + k2*h/2
    for (int i = 0 ; i < neq ; i++){
        
        Y1[i] = Y[i] + h*k2[i]*0.5;

    }

    RHSFunc(t+0.5*h,Y1,k3); // calcolo k3 con il RSH con t_n+h/2 e Y_n+k2*h/2
    
    // scrivo il ciclo per calcolare Y_{n+1} = Y_n + k3*h
    for (int i = 0 ; i < neq ; i++){
        
        Y1[i] = Y[i] + h*k3[i];

    }

    RHSFunc(t+h,Y1,k4); // calcolo k4 con il RSH con t_n+h e Y_n+k3*h
    
    // scrivo il ciclo per calcolare 
    // Y_{n+1} = Y_n + h/6 * ( k1 + 2*k2 + 2*k3 + k4 )
    for (int i = 0 ; i < neq ; i++){
        
        Y[i] += h * ( k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i] ) / 6.0;

    }

}

// metodo della secante
// gli do in input la funzione, gli estremi a e b, la tolleranza su x, 
// uno zero per riferimento e il numero di iterazioni
int secant_method(double (*F)(double), double a, double b, double tol, 
                  double &zero, int &l){

    // definisco le variabili che mi servono per tenere traccia delle
    // varie iterazioni di x
    double xk1 = a, xk = b, xk2 = xk + 1; // dove uso xk come x_k, xk1 come
    // x_{k-1} e xk2 come x_{k+1} ; inizializzo gli zeri sugli estremi 
    // dell'intervallo in cui ricaviamo la retta
    int n = 0; // variabile per contare
    // una variabile di controllo per vedere di quanto miglioriamo la guess
    double xp = 0;

    // definisco le variabili della funzione valutata
    double fa = F(a);
    double fb = F(b);
    double fxk;
    double fxk1;

    // metto i controlli di non avere gia' uno zero
    if( fa == 0.0 ){

        zero = a;

        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n 
             << "    (l'estremo " << a << " e' gia' lo zero)" << endl;

        l = n; // sono le iterazioni
        return 0;

    }
    else if( fb == 0.0 ){

        zero = b;

        // creo l'output delle iterazioni
        cout << "(Secant) # = " << n 
             << "(l'estremo " << b << " e' gia' lo zero)" << endl;

        l = n; // sono le iterazioni
        return 0;

    }
    else{
        
        // metto nel ciclo la condizione sia sulla tolleranza che sul numero 
        // di cicli (messo dopo)
        while( fabs( xk2 - xp ) > tol ){

            n++;

            if( n == 100 ){

                cout << "(Secant) Troppe iterazioni." << endl;

                l = n; // sono le iterazioni
                return 0;

            }

            xp = xk2;
            fxk = F(xk);
            fxk1 = F(xk1);

            // calcolo lo zero x_{k+1}
            xk2 = xk - fxk*( xk - xk1 )/( fxk - fxk1 );

            xk1 = xk;
            xk = xk2;

            // creo l'output voluto (esercizio froot.cpp)
            //cout << "n = " << n << ";   [a,b] = [" << xk1 << ", " << xk 
            //       << "];    x0 = " << xk2 << ";   Deltax = " << fabs(xk2-xk1) 
            //       << ";   f(x0) = " << F(xk2) << endl;

        }

        // creo l'output delle iterazioni
        //cout << "(Secant) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = xk2;
        return 0;

    }

}

// Implemento la funzione di TridiagSolver
// gli do in input i vettori contenenti rispettivamente: gli elementi sotto
// la diagonale, i termini noti, le soluzioni, r, gli elementi sopra la diagonale
// e la dimensione
void TridiagSolver(double *a, double *b, double *x, double *r, double *c, int n){

    // definisco i vettori h e p
    double *h = new double [n];
    double *p = new double [n];

    // calcolo gli elementi h[n] e p[n] per il metodo Tridiag
    // separo i termini patologici
    h[0] = c[0]/b[0];
    p[0] = r[0]/b[0];

    for( int i = 1 ; i < n ; i++ ){

        h[i] = c[i] / ( b[i] - a[i]*h[i-1] );

        p[i] = ( r[i] - a[i]*p[i-1] ) / ( b[i] - a[i]*h[i-1] );

    }

    // applico il metodo di risoluzione
    x[n-1] = p[n-1]; // termine patologico non avendo definito x[n+1]
    
    // applichiamo back-substitution
    for( int i = n-1 ; i >= 0 ; i-- ){

        x[i] = p[i] - h[i]*x[i+1];

    }

    // pulisco
    delete[] h;
    delete[] p;

}