// write a program that extract an integer random number in the range [1,100] and try to guess it. Based on the most recent guess, have the code suggest the interval containing the random number. Count the number of guesses. 
//
// For instance, suppose the number to guess is 35, then your code should produce an output similar to [input from keyboard is colored in red] (in the slide).
//

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

int main(){

    int xin = 0, xend = 100; // definisco i due numeri in cui genererò il numero casuale
    int n = 0; // la mia guess
    int number; // il numero che renderò casuale
    int count = 0; // conta i tentativi

    // Programma per indovinare un numero intero casuale
    srand(time(NULL));
    for(int i=0 ; i<1 ; i++){
        number = rand()%xend+xin;
    }

    while( n != number ){
        cout << "n in [" << xin << "," << xend << "]:\nyour guess: ";
        cin >> n;
        count += 1;
        if( n < xin || n > xend ){
            continue;
        }
        if( number > n ){
            xin = n;
        }
        if( number < n){
            xend = n;
        }
        if(number == n){
            cout << "Hai indovinato in " << count << " tentativi." << endl;
        }
    }

    return 0;
}