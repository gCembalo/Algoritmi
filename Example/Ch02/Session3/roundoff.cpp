#include <iostream>                                             // informazioni di input e output
#include <cmath>                                                // uso delle funzioni matematiche
#include <cstdlib>                                              // uso dei generatori casuali
#include <iomanip>                                              // manipolazione dell'output

using namespace std;

int main ()
{
    cout << setiosflags ( ios::scientific );                    // inserisce la notazione scientifica: e+# = 10^#

    int i;
    
    float x_r = 1.e+03, x_c = 1.;
    float fx1_r, fx2_r, f_taylor_r; 
    float fx1_c, fx2_c, f_taylor_c;

    cout << endl;
    cout << "Funzione sqrt(x^2+1)-x per x->00" << endl;
    cout << "-------------------------------------------------------------------------------------------------"  << endl;
           
    
    for ( i=0; i<7; i++ )
    {
        x_r = x_r*10;
        fx1_r = sqrt(x_r*x_r+1.)-x_r;
        fx2_r = (1.0)/(sqrt(x_r*x_r+1)+x_r);                    // moltiplicando e dividendo per sqrt(x*x+1)+x
        f_taylor_r = 0.5/x_r;

        cout << "x = " << x_r << "\t fx1 = " << fx1_r << "\t fx2 = " << fx2_r << "\t f(taylor) = " << f_taylor_r << endl;
    }

    cout << endl;
    cout << "Funzione 1-cos(x) per x->0" << endl;
    cout << "-------------------------------------------------------------------------------------------------"  << endl;


    for ( i=0; i<8; i++ )
    {
        x_c = x_c*0.1;
        fx1_c = 1-cos(x_c);
        fx2_c = 2*sin(x_c/2)*sin(x_c/2);                        // invertendo la formula di bisezione del seno
        f_taylor_c = x_c*x_c*0.5;
        
        cout << "x = " << x_c << "\t fx1 = " << fx1_c << "\t fx2 = " << fx2_c << "\t f(taylor) = " << f_taylor_c << endl;
    }
    
    cout << endl;

    return 0;
}