#include "my_header.h"

int Bisection ( double (*F)(double), double a, double b, double x_tol, double y_tol, double& zero_b, int& kb, int kmax )
{
    int err;

    double xm = 0., dx = 0.;
    double fa = F(a), fb = F(b), fxm = 0.;
    
    do
    {
        xm = (b+a)*0.5;
        dx = fabs(b-a);
        fxm = F(xm);

        if ( fa==0 )
        {
            dx = 0.;
            zero_b = a;
        }
        else if ( fb==0 )
        {
            dx = 0.;
            zero_b = b;
        }
        else if ( fa*fxm<0 )
        {
            b = xm;
            fb = fxm;
            zero_b = xm;
        }
        else if ( fa*fxm>0 )
        {
            a = xm;
            fa = fxm;
            zero_b = xm;
        }
        else
        {
            dx = 0.;
            zero_b = xm;
        }

        kb++;

        if ( kb>kmax) break;
    } 
    while ( dx>x_tol );

    if ( dx<=x_tol && fabs(fxm)<=y_tol )
    {
        err = 0;                                                                             
    }
    else
    {
        err = 1;                                                                                
    }

    return err;    
}

int FalsePos ( double (*F)(double), double a, double b, double x_tol, double y_tol, double& zero_f, int& kf, int kmax )
{
    int err;

    double del = 0., xm = 0.;
    double fa = F(a), fb = F(b), fxm = 0.;

    do
    {
        xm = F(a)*(a-b)/(fb-fa)+a;                                                        
        del = fabs(xm-zero_f);
        fxm = F(xm);

        if ( fa==0 )
        {
            del = 0.;
            zero_f = a;
        }
        else if ( fb==0 )
        {
            del = 0.;
            zero_f = b;
        }
        else if ( fa*fxm<0 )
        {
            b = xm;
            fb = fxm;
            zero_f = xm;
        }
        else if ( fa*fxm>0 )
        {
            a = xm;
            fa = fxm;
            zero_f = xm;
        }
        else
        {
            del = 0.;
            zero_f = xm;
        }

        kf++;

        if ( kf>kmax) break;
    } 
    while ( del>x_tol );

    if ( del<=x_tol && fabs(fxm)<=y_tol )
    {
        err = 0;                                                                             
    }
    else
    {
        err = 1;                                                                                
    }

    return err;    
}

int Secant ( double (*F)(double), double a, double b, double x_tol, double y_tol, double& zero_s, int& ks, int kmax )
{
    int err;

    double fa = F(a), fb = F(b);
    double dx = 0.;
    
    do
    {
        dx = fb*(b-a)/(fb-fa);
        a = b;
        fa = fb;
        b = b-dx;
        zero_s = b;
        fb = F(b);
        ks++;

        if ( ks>kmax) break;
    }
    while ( fabs(dx)>x_tol );

    if ( dx<=x_tol && fabs(fb)<=y_tol )
    {
        err = 0;                                                                             
    }
    else
    {
        err = 1;                                                                                
    }

    return err;  
}

int Newton ( double (*F)(double), double (*D)(double), double c, double x_tol, double y_tol, double& zero_n, int& kn, int kmax )
{
    int err;

    double f = F(c), d = D(c);
    double dx = 0.;
    
    do
    {
        dx = f/d;
        c = c-dx;
        zero_n = c;
        f = F(c);
        d = D(c);
        kn++;

        if ( kn>kmax) break;
    }
    while ( fabs(dx)>x_tol );

    if ( dx<=x_tol && fabs(f)<=y_tol )
    {
        err = 0;                                                                             
    }
    else
    {
        err = 1;                                                                                
    }

    return err;
}