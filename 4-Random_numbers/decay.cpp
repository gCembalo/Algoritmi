#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

int main(){

    ofstream fdata; // declare Output stream class to operate on files
    fdata.open("decay.dat"); // open output file

    //set the significant digits and the scientific notation
    fdata << setprecision(7);
    fdata << setiosflags ( ios::scientific );

    // definisco costanti
    double lambda = 0.01; // decay rate
    int delta_t = 1; // incremento tempo
    int ni = 1000; // atomi iniziali
    int nd, nr = 1000; // atomi decaduti e rimanenti
    double x; // probabilità di decadere
    double res;

    // ciclo temporale
    srand48(time(NULL));
    for( int t = 0 ; t < 500 ; t += delta_t){

        nd = 0; // inizializzo i decadimenti così li conto ad ogni ciclo

        for( int i = 0 ; i < nr ; i++){

            x = drand48(); // questa è la probabilità di decadere

            if( x <= lambda ){
                nd += 1;
            }

        }

        nr -= nd;

        if( nr <= 0 ) break;

        res = (double)nr/(double)ni; // stampo la percentuale di atomi rimasti
        fdata << t << "     " << res << endl;

    }
    fdata.close();

    return 0;
}