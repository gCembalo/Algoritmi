// using random deviates in [0,1] simulate the radioactive decay of an initial distribution of N atoms having a decay rate λ (use λ = 0.01).
// We increase time in discrete steps ∆t (= 1), and for each time interval we count the number of nuclei that have decayed during that ∆t.
//
//The simulation quits when there are no nuclei left to decay or when time exceeds a given threshold (e.g. t < 500). Such being the case, we have an outer loop over the time steps ∆t and an inner loop over the remaining nuclei for each time step.
//
// Produce a plot.
//

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