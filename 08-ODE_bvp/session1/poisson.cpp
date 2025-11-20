// we wish to find the electrostatic potential \Phi generated
// by a localized charge distribution ρ(r):
// \nabla^2 \Phi = -4\pi\rho
//

// SESSIONE GUIDATA
// Problema non omogeneo => la derivata è importante. usiamo la derivata s = dphi/dr (0) come parametro da indovinare
// ci mettismo in sferiche così passiamo da una pde ad una ode in r
// evitare le singolarità in coordinate sferiche. Definiamo \phi
// dobbiamo imporre condizioni sul potenziale, a 0 e a grandi distanze
// in un codice non possiamo mettere r=infty, ma possiamo vettere un punto tale per cui lo riteniamo grande abbastanza. Noi prendiamo r = 20
//
// Step delle silide (pag. 10)
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#define NMAX_EQ 64 // numero massimo di eq (sicurezza)

using namespace std;

void RK4Step(double, double *, void (*)(double, double *, double *), double, int);
void RHSFuncPoisson(double, double *, double *);
double ResidualPoisson(double);
int bisection(double (*)(double), double, double, double, double &, int &);

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata, fdata1;
    fdata.open("poisson.dat"); // file per le soluzioni

    //set the significant digits and the scientific notation
    fdata << setprecision(7);
    fdata << setiosflags ( ios::scientific );

    int neq = 2; // numero equazioni
    int npoint = 1000; // numero punti

    double r; // variabile raggio
    double r0 = 0.0; // raggio iniziale
    double rf = 20.0; // raggio finale
    double dr = fabs(rf - r0) / npoint; // step fisso

    // setto la condizioni iniziali
    double phi0 = 0.0;
    double phif = 1.0; // condizione a grandi distanze, ovvero per rf

    // variabile derivata di phi
    double s; // guess, derivata di phi in r
    double s0 = 0.0 , sf = 5.0;

    // definisco l'array delle soluzioni e imposto le condizioni iniziali
    double Y[neq];

    s = s0;

    // ----------------------------- STEP 1 ----------------------------- //
    for( int i = 0 ; i < 6 ; i++ ){

        // setto le condizioni iniziali
        Y[0] = phi0; // phi(r)
        Y[1] = s; // dphi/dr

        r = r0; // setto il punto di integrazione iniziale

        // plotto la condizione iniziale
        fdata << r << "  " << Y[0] << "  " << Y[1] << endl;

        // risolvo le equazioni del moto usando RK a 4 step
        for( int i = 0 ; i < npoint ; i++ ){

            RK4Step(r, Y, RHSFuncPoisson, dr, neq); // risolvo la ODE
            r += dr;

            // stampo nel file i dati
            fdata << r << "  " << Y[0] << "  " << Y[1] << endl;

        }

        fdata << endl << endl;

        s += 0.2; // cambio la guess sulla derivata

    }

    // ----------------------------- STEP 2 ----------------------------- //
    double ds = ( sf - s0 )/(double)npoint;
    double res; // variabile da stampare

    for( int i = 0 ; i < npoint ; i++){

        // richiamo la funzione residuo dando in pasto la s corrente
        res = ResidualPoisson(s);

        s = s0 + i * ds;

        fdata << s << "  " << res << endl;

    }

    fdata.close();

    // ----------------------------- STEP 3 ----------------------------- //
    double szero; // variabile per lo zero
    int l = 0; // variaible per le iterazioni di bisezione
    double tol = 1.e-9;

    bisection(ResidualPoisson, 0.0, 5.0, tol, s0, l);

    cout << s0 << endl;

    // ----------------------------- STEP 4 ----------------------------- //
    // confronto con la soluzione analitica


    return 0;

}

// definisco il Right-Hand-Side-Function (è problem dependent). 
// Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore)
// mi faccio dare R
void RHSFuncPoisson(double r, double *Y, double *R){

    double rho = exp(-r)/(8.0 * M_PI); // Scielded charge

    R[0] = Y[1];
    R[1] = - 4.0 * M_PI * r * rho;

}

// creo la funzione residuo che calcola quanto phi(r) si discosta da phi(20.0) = 1
// sarà la funzione che diamo alla funzione per la ricerca degli zeri
// è problem dependent.
// gli do in input la guess sulla derivata
double ResidualPoisson(double s){

    // ricopio esattamente il main precedente per trovare la soluzione con s
    int neq = 2; // numero equazioni
    int npoint = 1000; // numero punti

    double r; // variabile raggio
    double r0 = 0.0; // raggio iniziale
    double rf = 20.0; // raggio finale
    double dr = fabs(rf - r0) / npoint; // step fisso

    // setto la condizioni iniziali
    double phi0 = 0.0;
    double phif = 1.0; // condizione a grandi distanze, ovvero per rf

    // definisco l'array delle soluzioni e imposto le condizioni iniziali
    double Y[neq];

    // setto le condizioni iniziali
    Y[0] = phi0; // phi(r)
    Y[1] = s; // dphi/dr

    r = r0; // setto il punto di integrazione iniziale

    // risolvo le equazioni del moto usando RK a 4 step
    for( int i = 0 ; i < npoint ; i++ ){

        Y[0] = phi0;
        Y[1] = s;

        RK4Step(r, Y, RHSFuncPoisson, dr, neq); // risolvo la ODE
        r += dr;

    }

    // return il residuo
    return Y[0] - 1;

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

// metodo della bisezione
// gli do in input la funzione, gli estremi a e b, la tolleranza su x, 
// uno zero per riferimento e il numero di iterazioni
int bisection(double (*F)(double), double a, double b, double tol, 
              double &zero, int &l){

    double x; // la mia guess dello zero che aggiorno ad ogni iterazione
    int n = 0; // la variabile che mi permette di contare le iterazioni

    // definisco le variabili della funzione valutata
    double fa = F(a);
    double fb = F(b);
    double fx;

    // metto i controlli di non avere già uno zero
    if( fa == 0.0 ){

        zero = a;

        // creo l'output delle iterazioni
        cout << "(Bisection) # = " << n 
             << "    (l'estremo " << a << " è già lo zero)"<< endl;

        l = n;
        return 0;

    }
    else if( fb == 0.0 ){

        zero = b;

        // creo l'output delle iterazioni
        cout << "(Bisection) # = " << n 
             << "(l'estremo " << b << " è già lo zero)"<< endl;

        l = n;
        return 0;

    }
    else{

        while( fabs(a-b) > tol ){
        
            n++;

            // metto il controllo sul numero di iterazioni
            if( n == 100 ){

                cout << "(Bisection) Troppe iterazioni." << endl;

                l = n;
                return 0;

            }

            // calcolo la prima stima dello zero
            x = ( a+b ) * 0.5;

            // definisco le variabili delle funzioni valutate
            fa = F(a);
            fb = F(b);
            fx = F(x);

            // controllo se è uno zero
            if( fx == 0 ){

                zero = x;

                // creo l'output delle iterazioni
                //cout << "(Bisection) # = " << n << endl;

                l = n;
                return 0;

            }
            // controllo dove si trova x rispetto gli estremi a e b
            else if( fa*fx < 0 ){

                b = x;

            }
            else if ( fa*fx > 0 ){

                a = x;

            }

            // creo l'output voluto (esercizio froot.cpp)
            //cout << "n = " << n << ";   [a,b] = [" << a << ", " << b 
            //       << "];    xm = " << x << ";   Deltax = " << fabs(a-b) 
            //       << ";   f(xm) = " << F(x) << endl;

        }

        // creo l'output delle iterazioni
        //cout << "(Bisection) # = " << n << endl;

        l = n; // sono le iterazioni
        zero = x;
        return 0;

    }

}