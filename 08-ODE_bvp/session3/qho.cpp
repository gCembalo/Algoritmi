// find the eigenvalues of the quantum harmonic oscillator
// -1/2\dv[2]{\psi}{x} + V(x)\psi = E\psi
// con V = 1/2x^2
// Esercizio guidato con gli step nelle slide (pag. 15, 16, 17).
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#define NMAX_EQ 64 // numero massimo di eq (sicurezza)

static double g_E; // variabile globale per autostato fondamentale energia

using namespace std;

void RK4Step(double, double *, void (*)(double, double *, double *), double, int);
void RHSFuncQho(double, double *, double *);
double ResidualQho(double);
int bisection(double (*)(double), double, double, double, double &, int &);
int Bracket(double (*)(double), double, double, double, double *, double *, int &);

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata;
    fdata.open("qho.dat"); // file per le soluzioni

    //set the significant digits and the scientific notation
    fdata << setprecision(7);
    fdata << setiosflags ( ios::scientific );

    int neq = 2; // numero equazioni
    int nstep = 800; // numero punti

    double x; // variabile
    double x0 = -10.0; // variabile iniziale
    double xf = 10.0; // variabile finale
    double dxF = fabs(xf - x0) / nstep; // step avanti (forward)
    double dxB = -fabs(xf - x0) / nstep; // step avanti (Backward)

    x = x0; // setto il punto di integrazione iniziale

    g_E = 0.5; // setto l'autovalore dell'energia

    // definisco l'array delle soluzioni e imposto le condizioni iniziali
    double Y[neq];
    Y[0] = exp( -x*x*0.5 );
    Y[1] = -x*Y[0];

    // Integrazione Forward

    // risolvo le equazioni del moto usando RK a 4 step (Forward)
    for( int i = 0 ; i < nstep ; i++ ){

        RK4Step(x, Y, RHSFuncQho, dxF, neq); // risolvo la ODE
        x += dxF;

        // stampo nel file i dati
        fdata << x << "  " << Y[0] << "  " << Y[1] << endl;

    }

    fdata << endl << endl;

    // Integrazione Backward
    x = xf; // setto il punto di integrazione iniziale

    // imposto le condizioni iniziali
    Y[0] = exp( -x*x*0.5 );
    Y[1] = -x*Y[0];

    // risolvo le equazioni del moto usando RK a 4 step (Backward)
    for( int i = 0 ; i < nstep ; i++ ){

        RK4Step(x, Y, RHSFuncQho, dxB, neq); // risolvo la ODE
        x += dxB;

        // stampo nel file i dati
        fdata << x << "  " << Y[0] << "  " << Y[1] << endl;

    }

    fdata << endl << endl;

    return 0;

}


// definisco il Right-Hand-Side-Function (è problem dependent). 
// Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore)
// mi faccio dare R
void RHSFuncQho(double x, double *Y, double *R){

    // nota che il terzo elemento di Y ( Y[2] ) è k = -2*(E - V(x))
    // V(x) = x*x/2

    R[0] = Y[1];
    R[1] = -2 * ( g_E - 0.5*x*x ) * Y[0];

}

// creo la funzione residuo
// sarà la funzione che diamo alla funzione per la ricerca degli zeri
// è problem dependent.
// gli do in input la guess sulla derivata
double ResidualQho(double E){

    // ricopio esattamente il main precedente per trovare la soluzione con s
    int neq = 2; // numero equazioni
    int nstep = 1200; // numero punti
    int nstepL = 700 , nstepR = nstep - nstepL; // punti nei due intervalli

    double x; // variabile
    double x0 = -10.0; // variabile iniziale
    double xf = 10.0; // variabile finale

    double xm = 0.3; // matching point

    double dxF = fabs(xf - x0) / nstep; // step avanti (forward)
    double dxB = -fabs(xf - x0) / nstep; // step avanti (Backward)

    x = x0; // setto il punto di integrazione iniziale

    g_E = E; // setto l'autovalore dell'energia

    // definisco l'array delle soluzioni e imposto le condizioni iniziali
    double YF[neq] , YB[neq]; // Forward e Backward
    YF[0] = YB[0] = exp( -x*x*0.5 );
    YF[1] = -x*YF[0];
    YB[0] = -x*YB[0];

    // Integrazione Forward

    // risolvo le equazioni del moto usando RK a 4 step (Forward)
    for( int i = 0 ; i < nstep ; i++ ){

        RK4Step(x, Y, RHSFuncQho, dxF, neq); // risolvo la ODE
        x += dxF;

    }

    // Integrazione Backward
    x = xf; // setto il punto di integrazione iniziale

    // imposto le condizioni iniziali
    Y[0] = exp( -x*x*0.5 );
    Y[1] = -x*Y[0];

    // risolvo le equazioni del moto usando RK a 4 step (Backward)
    for( int i = 0 ; i < nstep ; i++ ){

        RK4Step(x, Y, RHSFuncQho, dxB, neq); // risolvo la ODE
        x += dxB;

    }

    // return il residuo (il valore di Y[1] è 1.0 in questo caso)
    return Y[0] - g_E;

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

// bracketing function
// gli do in input la funzione, gli estremi a e b, il numero di intervalli, 
// un array per gli estremi (sinistro e destro) degli intervalli, e un 
// riferimento al numero di zeri
int Bracket(double (*F)(double), double a, double b, double n, double *xL, 
            double *xR, int &nroots){

    // definisco le variabili che uso per contare
    int count = 0, i; // count mi dice quanti zeri ho
    double dx = ( b - a )/(double)n; // spacing
    double aL, aR; // estremi di ogni sotto intervallo

    double fL, fR; // valori delle funzioni valutate agli estremi
    fL = F(a);

    // faccio il loop su tutti gli intervalli
    for( i = 0 ; i < n ; i++ ){

        aL = a + i*dx;
        aR = a + (i+1)*dx;

        // metto la condizione in cui abbiamo un cambiamento di segno e quindi
        // con la quale ci ricordiamo il valore
        fR = F(aR);
        if( fL*fR <= 0.0 ){

            // metto gli estremi nell'array
            xL[count] = aL;
            xR[count] = aR;
            count++; // aggiorno il numero di roots

        }
        
        fL = fR;

    }

    nroots = count;
    return 0;

}