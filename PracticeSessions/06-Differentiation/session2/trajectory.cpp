// Given the particle trajectory
// alpha*t*t - t*t*t*( 1 - exp( -alpha*alpha/t ) )
// produce a plot of the velocity and acceleration in the range 0<t<α. How many inversion points are present ? (try α=10 to begin with).
//
// To this purpose, divide the range [0, α] into N equally spaced intervals Δt and use this spacing when computing the derivatives (that is, h = Δt).
// Note: x(t) has different L/R limits at t=0. Thus, central 1st and 2nd derivatives are ill-defined. At this point, replace them with the corresponding one-sided approximation.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

// variabile globale
double alpha = 10.0; // definita dal problema

double position(double);
double velocity(double);
double SecondDerivative(double (*)(double), double, double, double, double);
double derFD(double (*)(double), double, double, double);
double derBD(double (*)(double), double, double, double);
double derCD(double (*)(double), double, double, double);
double der4th(double (*)(double), double, double, double, double, double);

int main(){

    cout << setprecision(7);
    cout << setiosflags ( ios::scientific );

    ofstream fdata1, fdata2, fdata3 ,fdata4, fdata5, fdata6;
    fdata1.open("trajectory.dat");
    // per la derivata con FD
    fdata2.open("velocity1.dat");
    // per la derivata con BD
    fdata3.open("velocity2.dat");
    // per la derivata con CD
    fdata4.open("velocity3.dat");
    // per la derivata con 4therr
    fdata5.open("velocity4.dat");
    // accelerazione
    fdata6.open("acceleration.dat");

    //set the significant digits and the scientific notation
    fdata1 << setprecision(7);
    fdata1 << setiosflags ( ios::scientific );
    fdata2 << setprecision(7);
    fdata2 << setiosflags ( ios::scientific );
    fdata3 << setprecision(7);
    fdata3 << setiosflags ( ios::scientific );
    fdata4 << setprecision(7);
    fdata4 << setiosflags ( ios::scientific );
    fdata5 << setprecision(7);
    fdata5 << setiosflags ( ios::scientific );
    fdata6 << setprecision(7);
    fdata6 << setiosflags ( ios::scientific );

    // stampo la posizione in funzione del tempo in [0,alpha]
    for( double t = 0.0 ; t <= alpha ; t += 1.e-2 ){

        fdata1 << t << "     " << position(t) << endl;

    }
    fdata1.close();

    // definisco la spaziatura
    double dt = 1.e-2; // poi lo dimezzerò ad ogni ciclo
    // definisco le mie x_i , x_{i-1} e x_{i+1}
    double ti , tm , tp , tmm , tpp; // t_i è t_{i} ; tp è t_{i+1} ; tm è t_{i-1} ; tmm è t_{i-2} ; tpp è t_{i+2}

    // definisco le variabili per velocità e accelerazione
    double vel , acc;

    // stampo la posizione in funzione del tempo in [0,alpha]
    for( ti = 0.0 ; ti <= alpha ; ti += dt ){

        tm = ti - dt;
        tp = ti + dt;
        tmm = ti - 2.0*dt;
        tpp = ti + 2.0*dt;

        // calcolo la velocità con le diverse derivate
        vel = derFD(position, ti, tp, dt);
        fdata2 << ti << "     " << vel << endl;

        vel = derBD(position, ti, tm, dt);
        fdata3 << ti << "     " << vel << endl;

        // stampo la velocità con la derivata CD tenendo conto dei due limiti diversi a sinistra e destra
        if( ti*tm < 0 || ti*tmm < 0 ){

            // se abbiamo dei t<0 stampa la velocità con la derivata in avanti
            double velA = derFD(position, ti, tp, dt);
            fdata4 << ti << "     " << velA << endl;

        }
        else{

            vel = derCD(position, tm, tp, dt);
            fdata4 << ti << "     " << vel << endl;
        
        }

        // stampo la velocità con la 4thDer tenendo conto dei due limiti diversi a sinistra e destra
        if( ti*tm < 0 || ti*tmm < 0 ){

            // se abbiamo dei t<0 stampa la velocità con la derivata in avanti
            double velB = derFD(position, ti, tp, dt);
            fdata4 << ti << "     " << velB << endl;

        }
        else{

            vel = der4th(position, tmm, tm ,tp, tpp, dt);
            fdata5 << ti << "     " << vel << endl;
        
        }
        
        // stampo l'accelerazione tenendo conto dei due limiti diversi a sinistra e destra
        if( ti*tm < 0 || ti*tmm < 0 ){
            
            // se abbiamo dei t<0 stampa la velocità con la derivata in avanti
            double accA = derFD(velocity, ti, tp, dt);
            fdata6 << ti << "     " << accA << endl;

        }
        else{
            acc = SecondDerivative(position, ti, tm, tp, dt);
            fdata6 << ti << "     " << acc << endl;
        }

    }
    
    fdata2.close();
    fdata3.close();
    fdata4.close();
    fdata5.close();
    fdata6.close();

    return 0;

}

double position(double t){

    // elimino la singolarità
    if( t == 0 ){
        return 0.;
    }
    else{
        return alpha*t*t - t*t*t*( 1 - exp( -alpha*alpha/t ) );
    }

}

// definisco la velocità che mi serve solo per la derivata seconda in 0 avendo il problema di due limiti diversi
double velocity(double t){

    // elimino la singolarità
    if( t == 0 ){
        return 0.;
    }
    else{
        return 2*alpha*t - 3*t*t*( 1 - exp( -alpha*alpha/t ) ) + ( alpha*alpha*t )*exp( -alpha*alpha/t ) ;
    }

}

// x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}
double SecondDerivative(double (*F)(double x), double xi, double xm, double xp, double h){

    // definisco
    double fprimeprime = ( F(xp) - 2.0*F(xi) + F(xm) )/( h*h );

    return fprimeprime;
    
}

// x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}
double derFD(double (*F)(double), double xi, double xp, double h){

    // definisco 
    double fiPrime = ( F(xp) - F(xi) )/h;

    return fiPrime;

}

// x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}
double derBD(double (*F)(double), double xi, double xm, double h){

    // definisco 
    double fiPrime = ( F(xi) - F(xm) )/h;

    return fiPrime;

}

// x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}
double derCD(double (*F)(double), double xm, double xp, double h){

    // definisco 
    double fiPrime = ( F(xp) - F(xm) )/( 2.0*h );

    return fiPrime;

}

// x_i è x_{i} ; xp è x_{i+1} ; xm è x_{i-1} ; xmm è x_{i-2} ; xpp è x_{i+2}
double der4th(double (*F)(double), double xmm, double xm, double xp, double xpp, double h){

    // definisco 
    double fiPrime = ( F(xmm) - 8.0*F(xm) + 8.0*F(xp) - F(xpp) )/( 12.0*h );

    return fiPrime;

}