// Solve the equation of motion of a point mass in a central gravitational field:
//
// \dv[2]{x}{t} = -GM/r^2
//
// Use dimensionless units by setting GM=1 and set the initial mass at the point
// (x=4,y=0) with initial velocity in the y-direction.
// Use two dimensions (x-y) only.
//
// What is the maximum velocity for which the orbit is closed ?
// Consider first the case of a circular orbit and integrate for ≈10 orbits
// by counting turning points (1 orbit = 2 turning points). How can we
// safely choose the step size ?
// Now consider an elliptical orbit (try e.g. α = 0.3). 
//Can you devise a strategy to control the time step so that Dθ » const ?
// (Make sure your algorithm produces bounded orbits…)
// How would you scale your results to physical c.g.s units ?
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

void RHSFuncOde3(double, double *,double *);
void RK4Step(double, double *, void (*)(double, double *, double *), double, int);

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata;
    fdata.open("kepler.dat"); // file per le soluzioni

    //set the significant digits and the scientific notation
    fdata << setprecision(7);
    fdata << setiosflags ( ios::scientific );

    double GM = 1.0; // setto le unità adimensionali

    double alpha = 0.3; // definisco la variabile alpha che posso variare per avere orbite di tipo diverso; vedi la definizione della variabile v

    int neq = 4; // ordine della ode
    int norbit = 10; // numero di orbite
    int npoint = 1000; // numero di punti
    int count = 0; // conteggio i punti di inversione

    double t, r, v; // variabile tempo, radiale e velocità
    double dt = 0.4; // incremento

    double E0, E; // definisco le energie
    double dtheta = 0.1; // definisco il DeltaTheta costante
    double vold; // definisco una variabile per salvarmi la velocità, in modo da
    // fare in controllo sui turning point

    // definisco e inizializzo le condizioni iniziali
    double x0 , y0 , vx0 , vy0;
    x0 = 4.0;
    y0 = 0.0;
    r = sqrt( x0*x0 + y0*y0 );
    vx0 = 0.0;
    vy0 = sqrt( alpha / r ); // in generale possiamo scrivere v=sqrt(alpha / r)
    // e se alpha < 2 abbiamo orbite ellittiche
    v = sqrt( vx0*vx0 + vy0*vy0 );
    E0 = v*v*0.5 - 1.0/r; // calcolo energia meccanica iniziale (se negativa 
    // allora le orbite sono chiuse)

    // definisco l'array delle soluzioni e imposto le condizioni iniziali
    double Y[neq];
    Y[0] = x0;
    Y[1] = y0;
    Y[2] = vx0;
    Y[3] = vy0;

    // controllo che tipo di orbita ho e se il sistema può essere legato
    cout << "\n+----------------------------------------+\n" << endl;
    cout << "E = " << E << endl;
    cout << "v = " << v << endl;
    cout << endl;
    if ( v==sqrt(2./r) ) cout << "orbita parabolica" << endl;
    else if ( v>sqrt(2./r) ) cout << "orbita iperbolica" << endl;
    else if ( v==sqrt(1./r) ) cout << "orbita circolare" << endl;
    else cout << "orbita ellittica" << endl;
    //cout << "\n+----------------------------------------+\n" << endl;

    // risolvo le equazioni del moto usando RK a 4 step
    for( int i = 0 ; i < npoint ; i++ ){

        // posso definire gli step temporali in base al punto dell'orbita in 
        // cui mi trovo
        dt = dtheta*r/v; // in cui abbiamo usato:
        // v = ds/dt  =>  ds = v*dt  =>  ds = r*dtheta = v*dt  =>  
        // => dt = r*dtheta / v

        // salvo la vecchia velocità per poter fare il controllo del 
        // turning point
        vold = Y[3];

        RK4Step(t, Y, RHSFuncOde3, dt, neq); // risolvo la ODE
        t += dt;

        // controllo e conteggio il turning point (se presente)
        if( Y[3]*vold < 0 ){

            count ++;
            //cout << "RK step = " << i << " turning point #" << count << endl;

            // fermo il ciclo se subero i 2 turning point per orbita
            if( count == 2*norbit + 1 ){

                break;

            }

        }

        // calcolo l'energia
        r = sqrt( Y[0]*Y[0] + Y[1]*Y[1] );
        v = sqrt( Y[2]*Y[2] + Y[3]*Y[3] );
        E = v*v*0.5 - 1/r;

        // stampo nel file i dati
        fdata << Y[0] << "  " << Y[1] << "  " << Y[2] << "  " << Y[3] << endl;

    }

    fdata.close();

    // mi faccio dire quanti turning point ci sono
    cout << "\n Ci sono #" << count << " turning point in " << norbit 
         << " orbite." << endl;
    cout << "\n+----------------------------------------+\n" << endl;

    return 0;

}

// definisco il Right-Hand-Side-Function (è problem dependent). 
// Gli do in input t e il puntatore ad Y e in uscita (tramite il puntatore) mi
// faccio dare R
void RHSFuncOde3(double t, double *Y, double *R){

    // setto le condizioni iniziali
    double x = Y[0];
    double y = Y[1];
    double vx = Y[2];
    double vy = Y[3];

    double r = sqrt( x*x + y*y );

    R[0] = vx;
    R[1] = vy;
    R[2] = -x/( r*r*r );
    R[3] = -y/( r*r*r );

}

// implemento il metodo Runge-Kutta del quarto ordine.
// gli do in input la variabile di integrazione, il puntatore alle soluzioni,
// il puntatore alla funzione del Right-Hand-Side-Function, l'incremento e
// l'ordine della ODE.
void RK4Step(double t, double *Y, void (*RHSFunc)(double t, double *Y, double *R),
             double h, int neq){
    
    // definisco i vettori per gli step intermedi
    double Y1[neq], k1[neq], k2[neq], k3[neq], k4[neq];
    
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