#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

int main ()
{
    int i, j;

    float m, k = 1.;                                    //    m = uno per la macchina    k = costante matematica uno
    
    double mm, kk = 1.;                                 //    doppia lettera per la doppia precisione
  
    while ( m!=k )
    {
        m = 1+1/pow(10,i);
        i++;
    } 

    while ( mm!=kk )
    {
        mm = 1+1/pow(10,j);
        j++;
    }

    cout << endl;
    cout << "precisione float: \t" << i-1 << endl;
    cout << "precisione double: \t" << j-1 << endl;
    cout << endl;

    return 0;
}