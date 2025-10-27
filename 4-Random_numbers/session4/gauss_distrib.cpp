// Using the acceptance-rejection method, generate random deviates following a Gaussian-like distribution. with μ=0, σ=0.5.
//
// Choose C to be a constant, e.g. C ≥ W(x) = f(x) and choose x in the range [-5,5];
//
// Generate pairs of random numbers (xk,yk) for k = 0…N-1 (take N = 105) and construct the distribution.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double gaussiana(const double&, const double&);

int main(){

    ofstream fdata; // declare Output stream class to operate on files
    fdata.open("gaussian.dat"); // open output file

    //set the significant digits and the scientific notation
    fdata << setprecision(7);
    fdata << setiosflags ( ios::scientific );

    double sigma = 0.5;

    double x, y;
    double C = 2. / ( sigma * sqrt(2*3.1415926)); // la nostra funzione limite

    srand48(time(NULL));
    for( int i = 0 ; i < 1e5 ; i++ ){

        x = -5 + drand48()*10; // genero i numeri casuali tra -5 e 5
        y = drand48()*C;

        if( y < gaussiana(x, sigma) ){
            fdata << x << "     " << y << endl;
        }

    }
    fdata.close();

    return 0;
}




double gaussiana(const double &x, const double &sigma){
    return 1./(sigma*sqrt(2.*3.1415926)) * exp( -0.5 * (x*x) / (sigma*sigma) );
}