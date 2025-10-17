#include <iostream>
using namespace std;

int main(){

    // Variabili per i conteggi
    int n = 0, m = 0;

    // Ciclo con i double
    double i = 1.0;
    double epsilon = 0.1;

    while(i+epsilon > i){
        epsilon = epsilon/10;
        n++;
    }

    //cout << setprecision(10);
    cout << "+---------------------------------------+" << endl;
    cout << "L'ultima cifra (float) aggiunta è: " << epsilon << endl;
    cout << "Dunque la precisione è: " << n << endl;
    cout << "+---------------------------------------+" << endl;

    // Ciclo con i float
    float j = 1.0;
    float psi = 0.1;

    while(j+psi > j){
        psi = psi/10;
        m++;
    }

    //cout << setprecision(10);
    cout << "L'ultima cifra (double) aggiunta è: " << psi << endl;
    cout << "Dunque la precisione è: " << m << endl;
    cout << "+---------------------------------------+" << endl;

    return 0;
}