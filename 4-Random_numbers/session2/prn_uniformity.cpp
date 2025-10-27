// test a random-number generator in the range [0,1] (use double precision !) to obtain a numerical measure of its uniformity and randomness before you stake your scientific reputation on it. 
//
// Generate a quick visual test by plotting ri as a function of i (0 ≤ i < 103) : a uniform distribution of uncorrelated point values should appear.
//
// Evaluate the kth -moment of the distribution. If the number are distributed uniformly, then the previous equation can be approximated by (something in the slide).
//
// Check that the deviation varies as 1/sqrt(N) by plotting the error as a function of N (double N at each iteration to produce a log-log plot). Limit to k = 1 and k = 2.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

int main(){

    ofstream fdata1, fdata2; // declare Output stream class to operate on files
    fdata1.open("uniformity.dat"); // open output file
    fdata2.open("prn_uniformity.dat"); // open output file

    //set the significant digits and the scientific notation
    fdata1 << setprecision(7);
    fdata1 << setiosflags ( ios::scientific );
    fdata2 << setprecision(7);  //set the significant digits
    fdata2 << setiosflags ( ios::scientific );  //set the scientific notation
    
    double x; // variabile che rendo casuale

    // creo il file di numeri semi-randomici per test
    srand48(time(NULL));
    for( int k=0 ; k<1e3 ; k++){
        x = drand48();

        // stampo i numeri
        fdata1 << k << "     " << x << endl;
    }
    fdata1.close();

    // questi sono le approssimazioni delle medie
    double funck1 = 1./2.; // k=1
    double funck2 = 1./3.; // k=2

    for( int i=4 ; i<1e5 ; i*=2 ){

        double x1 = 0., x2 = 0.; // definisco i momenti

        for( int j=0 ; j<i ; j++){

            x = drand48();

            // calcolo i momenti
            x1 += x;
            x2 += x*x;
        }

        x1 /= (double)i;
        x2 /= (double)i;

        double errx1 = fabs( x1 - funck1 );
        double errx2 = fabs( x2 - funck2 );

        fdata2 << i << "    " << errx1 << "     " << errx2 << endl;
    }
    fdata2.close();  // Chiudo il file

    return 0;
}