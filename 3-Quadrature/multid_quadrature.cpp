#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double func(double, double);
double func2(double, double);
double Gauss2D(double (*)(double, double), double, double, double, double, int, int);
void ConvergenceGauss(double (*F)(double, double), double , double , double , double , int , double);

int main(){
    cout << setprecision(12);  //set the significant digits
    cout << setiosflags ( ios::scientific );  //set the scientific notation

    // Dominio di integrazione
    double x0 = -1. , y0 = -1.;
    double x1 = 1. , y1 = 1.;

    // Risultato esatto
    double ris = 412.0/45.0;
    double tol = 1e-5;

    cout << "\n+-----------------------------------------+\n" << endl;
    cout << "         x^4y^2 + 2x^2y^2 - x^2y + 2       \n" << endl;
    cout << "Risultato esatto = " << ris << endl;
    cout << "Gauss            = " << Gauss2D(func, x0, x1, y0, y1, 1, 3) << endl;

    cout << "\n+-----------------------------------------+\n" << endl;
    cout << "         valore pi funzione a gradino       \n" << endl;
    cout << "pi  = 3.1415926535897932384626433\n" << endl;
    cout << "Con una tolleranza = " << setprecision(0) << tol << " : " << endl;
    cout << setprecision(12);  //set the significant digits
    ConvergenceGauss(func2, x0, x1, y0, y1, 4, tol);

    cout << "\n+-----------------------------------------+\n" << endl;

}


double func(double x, double y){
    return x*x*x*x*y*y + 2*x*x*y*y - x*x*y + 2;
}

double func2(double x, double y){
    if( sqrt(x*x + y*y) <= 1 ){
        return 1;
    }
    else{
        return 0;
    }
}

double Gauss2D(double (*F)(double, double), double x0, double x1, double y0, double y1, int n, int ng){

    // Implemento la funzione per funzionare fino a 5 punti
    double w[ng], x[ng]; // array di pesi e zeri del polinomio di Legendre

    // Riempio gli array in base al numero di punti che abbiamo
    if(ng == 1){
        x[0] = 0.;
        w[0] = 2.;
    }
    else if(ng == 2){
        x[1] = sqrt(1./3.);
        x[0] = -sqrt(1./3.);
        w[0] = 1.;
        w[1] = 1.;
    }
    else if(ng == 3){
        x[0] = -sqrt(3./5.);
        x[1] = 0.;
        x[2] = sqrt(3./5.);
        w[0] = 5./9.;
        w[1] = 8./9.;
        w[2] = 5./9.;
    }
    else if(ng == 4){
        x[0] = -sqrt( 3./7. - 2./7.*sqrt(6./5.) );
        x[1] = -sqrt( 3./7. + 2./7.*sqrt(6./5.) );
        x[2] = sqrt( 3./7. - 2./7.*sqrt(6./5.) );
        x[3] = sqrt( 3./7. + 2./7.*sqrt(6./5.) );
        w[0] = ( 18. + sqrt(30.) )/( 36. );
        w[1] = ( 18. - sqrt(30.) )/( 36. );
        w[2] = ( 18. + sqrt(30.) )/( 36. );
        w[3] = ( 18. - sqrt(30.) )/( 36. );
    }
    else if(ng == 5){
        x[0] = -1./3.*sqrt( 5. - 2.*sqrt(10./7.) );
        x[1] = -1./3.*sqrt( 5. + 2.*sqrt(10./7.) );
        x[2] = 0.;
        x[3] = 1./3.*sqrt( 5. - 2.*sqrt(10./7.) );
        x[4] = 1./3.*sqrt( 5. + 2.*sqrt(10./7.) );

        w[0] = ( 322. + 13.*sqrt(70.) )/( 900. );
        w[1] = ( 322. - 13.*sqrt(70.) )/( 900. ); 
        w[2] = 128./225.;
        w[3] = ( 322. + 13.*sqrt(70.) )/( 900. );
        w[4] = ( 322. - 13.*sqrt(70.) )/( 900. );
    }
    else{
        cout << "codice non implementato per Ng>5." << endl;
        return 0.;
    }

    // Calcolo gli integrali
    double dx = fabs(x0-x1)/(double)n;
    double dy = fabs(y0-y1)/(double)n;

    double a , b , s1 , s2;
    double c , d , s3 , s4;

    double sum = 0.0;
    for( int i = 0 ; i < n ; i++ ){  // Ciclo per le x

        // Definisco gli estremi dell'intervallo
        a = x0 + i*dx;
        b = a + dx;
        s1 = (b - a)/2;
        s2 = (b + a)/2;

        for( int j = 0 ; j < n ; j++ ){   // Ciclo per le y
            
            // Definisco le variabili dell'intervallo
            c = y0 + j*dy;
            d = c + dy;
            s3 = (d - c)/2;
            s4 = (d + c)/2;

            for( int ik = 0 ; ik < ng ; ik ++ ){   // Ciclo per variare x

                for( int jk = 0 ; jk < ng ; jk ++){   // Ciclo per variare y
                    sum += w[ik] * s1 * s3 * w[jk] * F( s1*x[ik] + s2 , s3*x[jk] + s4 );
                }

            }

        }

    }

    return sum ;
}

void ConvergenceGauss(double (*F)(double, double), double x0, double x1, double y0, double y1, int ng, double tol){
    // Controllo quanti intervalli servono
    int i = 1;
    while ( fabs( Gauss2D(func2, x0, x1, y0, y1, i, ng) - 3.1415926535897932384626433 ) > tol ){
        i ++ ;
    }
    cout << "Gauss = " << Gauss2D(func2, x0, x1, y0, y1, i, ng) << "\ncon " << i << " intervalli." << endl;
}