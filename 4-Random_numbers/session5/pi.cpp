#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double gaussiana(const double&, const double&);

int main(){

    cout << setprecision(7);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    ofstream fdata, fdata1, fdata2; // declare Output stream class to operate on files
    fdata.open("pi.dat"); // open output file
    fdata1.open("out_pi.dat"); // open output file
    fdata2.open("err_pi.dat"); // open output file

    //set the significant digits and the scientific notation
    fdata << setprecision(7);
    fdata << setiosflags ( ios::scientific );
    fdata1 << setprecision(7);
    fdata1 << setiosflags ( ios::scientific );
    fdata2 << setprecision(7);
    fdata2 << setiosflags ( ios::scientific );

    double x, y;
    int n = 4, n_in = 0;
    double tol = 1e-4;
    double I; // area del cerchio che calcolo
    double err = 3.0; // errore sull'area

    srand48(time(NULL));

    // faccio il ciclo che continua finché l'errore sull'area non è minore della tolleranza
    while( err > tol and n < 1e7 ){

        n_in = 0;

        for( int i=0 ; i<n ; i++ ){

            // genero coppie di numeri tra -1 e +1
            x = -1 + drand48()*2;
            y = -1 + drand48()*2;

            if( x*x + y*y <= 1){
                n_in ++;
                fdata << x << "     " << y << endl;
            }
            else{
                fdata1 << x << "     " << y << endl;
            }
        }

        // calcolo l'area del cerchio con l'approx I = A*n_in/n
        I = 4.0 * (double)n_in / (double)n;
        err = fabs( I/M_PI - 1.0 );

        fdata2 << n << "    " << err << endl;

        n *= 2;
        cout << "sono a n = " << n << endl;

    }
    fdata.close();
    fdata1.close();
    fdata2.close();

    cout << "\n+------------------------------------------------+\n" << endl;
    cout << "I punti totali sono: " << n << endl;
    cout << "Vengono all'interno del cerchio unitario: " << n_in << endl;
    cout << "\nIl rapporto N_in/N è = " << (double)n_in / (double)n << endl;
    cout << "\nL'area del cerchio è = " << I << "con errore = " << err << endl;
    cout << "\n+------------------------------------------------+\n" << endl;

    return 0;
}