#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

ofstream fdata;

void Boundary ( double**, double, double, int, int );
int Jacobi ( double**, double**, double*, double*, double, double, double, double, int, int );
int Gauss_Seidel ( double**, double**, double*, double*, double, double, double, double, int, int );
int SOR ( double**, double**, double*, double*, double, double, double, double, int, int );

int main ()
{   
    cout << setiosflags ( ios::scientific );
    cout << setprecision ( 2 );

    fdata << setiosflags ( ios::scientific );
    fdata << setprecision ( 10 );

    int i, j, k;
    int n = 128, m = 64;
    int kj, kg, ks;

    double dx, dy;
    double S = 0., tol = 1.e-7;
    double xi = 0., xf = 2.;
    double yi = 0., yf = 1.;
    double x[n], y[m];
    double **phi0, **phi1;

    phi0 = new double*[n];
    phi0[0] = new double[n*n];

    dx = fabs(xf-xi)/(double)(n-1);
    dy = fabs(yf-yi)/(double)(m-1);

    for ( i=1; i<n; i++ )
    {
        phi0[i] = phi0[i-1]+n;
    }

    phi1 = new double*[n];
    phi1[0] = new double[n*n];

    for( i=1; i<n; i++ )
    {
        phi1[i] = phi1[i-1]+n;
    }

    cout << endl;

    fdata.open ( "conductivity.dat" );

    for ( i=0; i<n; i++ )
    {
        x[i] = xi+i*dx;
    }

    for ( j=0; j<m; j++ )
    {
        y[j] = yi+j*dy;
    }   

    for ( i=0; i<n; i++ )
    { 
        for ( j=0; j<m; j++ )
        {
            phi0[i][j] = 0.;
            phi1[i][j] = 0.; 
        }
    }

    kj = Jacobi ( phi0, phi1, x, y, S, dx, dy, tol, n, m );
        
    for ( i=0; i<n; i++ )
    {
        for ( j=0; j<m; j++ )
        {
            fdata << x[i] << " " << y[j] << " " << phi1[i][j] << endl;
        }
            
        fdata << endl << endl;
    }
        
    cout << "Jacobi: \t k = " << kj << endl;

    for ( i=0; i<n; i++ )
    {
        x[i] = xi+i*dx;
    }

    for ( j=0; j<m; j++ )
    {
        y[j] = yi+j*dy;
    }   

    for ( i=0; i<n; i++ )
    { 
        for ( j=0; j<m; j++ )
        {
            phi0[i][j] = 0.;
            phi1[i][j] = 0.; 
        }
    }

    kg = Gauss_Seidel ( phi0, phi1, x, y, S, dx, dy, tol, n, m );
        
    for ( i=0; i<n; i++ )
    {
        for ( j=0; j<m; j++ )
        {
            fdata << x[i] << " " << y[j] << " " << phi1[i][j] << endl;
        }
            
        fdata << endl << endl;
    }
        
    cout << "GaussSeidel: \t k = " << kg << endl;

    for ( i=0; i<n; i++ )
    {
        x[i] = xi+i*dx;
    }

    for ( j=0; j<m; j++ )
    {
        y[j] = yi+j*dy;
    }   

    for ( i=0; i<n; i++ )
    { 
        for ( j=0; j<m; j++ )
        {
            phi0[i][j] = 0.;
            phi1[i][j] = 0.; 
        }
    }

    ks = SOR ( phi0, phi1, x, y, S, dx, dy, tol, n, m );

    for ( i=0; i<n; i++ )
    {
        for ( j=0; j<m; j++ )
        {
            fdata << x[i] << " " << y[j] << " " << phi1[i][j] << endl;
        }
            
        fdata << endl << endl;
    }
        
    cout << "SOR: \t\t k = " << ks << endl;

    cout << endl;

    fdata.close ();

    delete[] phi0;
    delete[] phi1;

    return 0;
}

void Boundary ( double** phi, double dx, double dy, int n, int m )
{
    int i, j;
    
    for ( i=0; i<n; i++ )
    {
        phi[i][0] = 0.;
        phi[i][m-1] = 2-i*dx;
    }

    for ( j=0; j<m; j++ )
    {
        phi[0][j] = phi[1][j];
        phi[n-1][j] = phi[n-2][j]+3*dx;
    }
}

int Jacobi ( double** phi0, double** phi1, double* x, double* y, double S, double dx, double dy, double tol, int n, int m )
{
    int i, j, k = 0;

    double sum = 0., error = 1.;

    while ( error>tol )
    {
        sum = 0.;

        Boundary ( phi0, dx, dy, n, m );
        Boundary ( phi1, dx, dy, n, m );

        for ( i=0; i<n; i++ )
        { 
            for ( j=0; j<m; j++ )
            {
                phi0[i][j] = phi1[i][j]; 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1;  j<m-1;  j++ )
            {
                phi1[i][j] = 0.25*(phi0[i+1][j]+phi0[i-1][j]+phi0[i][j+1]+phi0[i][j-1]-dx*dy*S); 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1; j<m-1; j++ )
            {
                //dx = phi1[i+1][j]-2*phi1[i][j]+phi1[i-1][j];
                //dy = phi1[i][j+1]-2*phi1[i][j]+phi1[i][j-1];
                sum = sum+fabs(phi1[i][j]-phi0[i][j])*dx*dy;
            }
        }
   
        error = sum;
        k = k+1;
    }

    return k;
}

int Gauss_Seidel ( double** phi0, double** phi1, double* x, double* y, double S, double dx, double dy, double tol, int n, int m )
{
    int i, j, k = 0;

    double sum = 0., error = 1.;

    while ( error>tol )
    {
        sum = 0.;

        Boundary ( phi0, dx, dy, n, m );
        Boundary ( phi1, dx, dy, n, m );

        for ( i=0; i<n; i++ )
        { 
            for ( j=0; j<m; j++ )
            {
                phi0[i][j] = phi1[i][j]; 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1; j<m-1; j++ )
            {
                phi1[i][j] = 0.25*(phi0[i+1][j]+phi1[i-1][j]+phi0[i][j+1]+phi1[i][j-1]-dx*dy*S); 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1; j<m-1; j++ )
            {
                //dx = phi1[i+1][j]-2*phi1[i][j]+phi1[i-1][j];
                //dy = phi1[i][j+1]-2*phi1[i][j]+phi1[i][j-1];
                sum = sum+fabs(phi1[i][j]-phi0[i][j])*dx*dy;
            }
        }
   
        error = sum;
        k = k+1;  
    }

    return k;
}

int SOR ( double** phi0, double** phi1, double* x, double* y, double S, double dx, double dy, double tol, int n, int m )
{
    int i, j, k = 0;
    
    double sum = 0., error = 1., omega;

    omega = 2/(1+M_PI/(double)n);
    
    while ( error>tol )
    {
        sum =0.;
        
        Boundary ( phi0, dx, dy, n, m );
        Boundary ( phi1, dx, dy, n, m );

        for ( i=0; i<n; i++ )
        { 
            for ( j=0; j<m; j++ )
            {
                phi0[i][j] = phi1[i][j]; 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1; j<m-1; j++ )
            {
                phi1[i][j] = (1-omega)*phi1[i][j]+0.25*omega*(phi0[i+1][j]+phi1[i-1][j]+phi0[i][j+1]+phi1[i][j-1]-dx*dy*S); 
            }
        }

        for ( i=1; i<n-1; i++ )
        { 
            for ( j=1; j<m-1; j++ )
            {
                //dx = phi1[i+1][j]-2*phi1[i][j]+phi1[i-1][j];
                //dy = phi1[i][j+1]-2*phi1[i][j]+phi1[i][j-1];
                sum = sum+fabs(phi1[i][j]-phi0[i][j])*dx*dy;
            }
        }
   
        error = sum;
        k = k+1;  

        cout.flush ();
    }

    return k;
}