#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

ofstream fdata;

void Boundary ( double**, double*, double*, int, double );

int Jacobi ( double**, double**, double*, double*, double, double, double, int );
int Gauss_Seidel ( double**, double**, double*, double*, double, double, double, int );
int SOR ( double**, double**, double*, double*, double, double, double, int );

double Solution ( double, double, double );

int main ()
{   
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 2 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 10 );

    int i, j, k;
    int n = 32;
    int kj, kg, ks;

    double h, S, tol = 1.e-7;
    double xi = 0., xf = 1.;
    double yi = 0., yf = 1.;
    double x[n], y[n];
    double **phi0, **phi1;

    phi0 = new double*[n];
    phi0[0] = new double[n*n];

    h = fabs(xf-xi)/(double)(n-1);
    
    for ( i=1; i<n; i++ )
    {
        phi0[i] = phi0[i-1]+n;
    }

    phi1 = new double*[n];
    phi1[0] = new double[n*n];

    for ( i=1; i<n; i++ )
    {
        phi1[i] = phi1[i-1]+n;
    }

    cout << endl;

    fdata.open ( "elliptic.dat" );

    for ( k=0; k<2; k++ )
    {
        for ( i=0; i<n; i++ )
        {
            x[i] = xi+i*h;
            y[i] = yi+i*h;
            
            for ( j=0; j<n; j++ )
            {
                phi0[i][j] = 0.;
                phi1[i][j] = 0.; 
            }
        }
        
        S = (double)k*2;
        kj = 0;
        kj = Jacobi ( phi0, phi1, x, y, S, h, tol, n );
        
        for ( i=0; i<n; i++ )
        {
            for ( j=0; j<n; j++ )
            {
                fdata << x[i] << " " << y[j] << " " << phi1[i][j] << endl;
            }
            
            fdata << endl << endl;
        }
        
        cout << "S = " << S << "\t Jacobi: \t k = " << kj << endl;
    }

    for ( k=0; k<2; k++ )
    {
        for ( i=0; i<n; i++ )
        {
            x[i] = xi+i*h;
            y[i] = yi+i*h;
            
            for ( j=0; j<n; j++ )
            {
                phi0[i][j] = 0.;
                phi1[i][j] = 0.; 
            }
        }
        
        S = (double)k*2;
        kg = 0;
        kg = Gauss_Seidel ( phi0, phi1, x, y, S, h, tol, n );
        
        for ( i=0; i<n; i++ )
        {
            for ( j=0; j<n; j++ )
            {
                fdata << x[i] << " " << y[j] << " " << phi1[i][j] << endl;
            }
            
            fdata << endl << endl;
        }
        
        cout << "S = " << S << "\t GaussSeidel: \t k = " << kg << endl;
    }

    for ( k=0; k<2; k++ )
    {
        for ( i=0; i<n; i++ )
        {
            x[i] = xi+i*h;
            y[i] = yi+i*h;
            
            for ( j=0; j<n; j++ )
            {
                phi0[i][j] = 0.;
                phi1[i][j] = 0.; 
            }
        }
        
        S = (double)k*2;
        ks = 0;
        ks = SOR ( phi0, phi1, x, y, S, h, tol, n );
        
        for ( i=0; i<n; i++ )
        {
            for ( j=0; j<n; j++ )
            {
                fdata << x[i] << " " << y[j] << " " << phi1[i][j] << endl;
            }
            
            fdata << endl << endl;
        }
        
        cout << "S = " << S << "\t SOR: \t\t k = " << ks << endl;
    }

    cout << endl;

    fdata.close ();

    delete[] phi0;
    delete[] phi1;

    return 0;
}

double Solution ( double x, double y, double S )
{
    return exp(-M_PI*x)*sin(-M_PI*y)+S*0.25*(x*x+y*y);
}

void Boundary ( double** phi, double* x, double* y, int n, double S )
{
    int i;
    
    for( i=0; i<n; i++ )
    {
        phi[i][0] = Solution(x[i],y[0],S);
        phi[i][n-1] = Solution(x[i],y[n-1],S);
        phi[0][i] = Solution(x[0],y[i],S);
        phi[n-1][i] = Solution(x[n-1],y[i],S);
    }
}

int Jacobi ( double** phi0, double** phi1, double* x, double* y, double S, double h, double tol, int n )
{
    int i, j, k = 0;

    double sum = 0., error = 1., dx = 0., dy = 0.;

    while ( error>tol )
    {
        sum = 0.;

        Boundary ( phi0, x, y, n, S );
        Boundary ( phi1, x, y, n, S );

        for ( i=0; i<n; i++ )
        { 
            for ( j=0; j<n; j++ )
            {
                phi0[i][j] = phi1[i][j]; 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1;  j<n-1;  j++ )
            {
                phi1[i][j] = 0.25*(phi0[i+1][j]+phi0[i-1][j]+phi0[i][j+1]+phi0[i][j-1]-h*h*S); 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1; j<n-1; j++ )
            {
                dx = phi1[i+1][j]-2*phi1[i][j]+phi1[i-1][j];
                dy = phi1[i][j+1]-2*phi1[i][j]+phi1[i][j-1];
                sum = sum+fabs(dx+dy-h*h*S);
            }
        }
   
        error = sum;
        k = k+1;
    }

    return k;
}

int Gauss_Seidel ( double** phi0, double** phi1, double* x, double* y, double S, double h, double tol, int n )
{
    int i, j, k = 0;

    double sum = 0., error = 1., dx = 0., dy = 0.;

    while ( error>tol )
    {
        sum = 0.;

        Boundary ( phi0, x, y, n, S );
        Boundary ( phi1, x, y, n, S );

        for ( i=0; i<n; i++ )
        { 
            for ( j=0; j<n; j++ )
            {
                phi0[i][j] = phi1[i][j]; 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1; j<n-1; j++ )
            {
                phi1[i][j] = 0.25*(phi0[i+1][j]+phi1[i-1][j]+phi0[i][j+1]+phi1[i][j-1]-h*h*S); 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1; j<n-1; j++ )
            {
                dx = phi1[i+1][j]-2.*phi1[i][j]+phi1[i-1][j];
                dy = phi1[i][j+1]-2.*phi1[i][j]+phi1[i][j-1];
                sum = sum + fabs(dx+dy-h*h*S);
            }
        }
   
        error = sum;
        k = k+1;  
    }

    return k;
}

int SOR ( double** phi0, double** phi1, double* x, double* y, double S, double h, double tol, int n )
{
    int i, j, k = 0;
    
    double sum = 0., error = 1., dx = 0., dy = 0., omega;

    omega = 2./(1.+M_PI/(double)n);
    
    while ( error>tol )
    {
        sum =0.;
        
        Boundary ( phi0, x, y, n, S );
        Boundary ( phi1, x, y, n, S );

        for ( i=0; i<n; i++ )
        { 
            for ( j=0; j<n; j++ )
            {
                phi0[i][j] = phi1[i][j]; 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1; j<n-1; j++ )
            {
                phi1[i][j] = (1-omega)*phi0[i][j]+0.25*omega*(phi1[i-1][j]+phi0[i+1][j]+phi0[i][j+1]+phi1[i][j-1]-h*h*S); 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1; j<n-1; j++ )
            {
                dx = phi1[i+1][j]-2.*phi1[i][j]+phi1[i-1][j];
                dy = phi1[i][j+1]-2.*phi1[i][j]+phi1[i][j-1];
                sum = sum + fabs(dx+dy-h*h*S);
            }
        }
   
        error = sum;
        k = k+1;  
    }

    return k;
}