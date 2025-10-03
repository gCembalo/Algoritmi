#include <iostream> 

using namespace std;

double AddOne ( double );

int main ()             
{
    double a;

    cout << endl;
    cout << "Inserire un numero reale" << endl;

    cin >> a;
    
    cout << "------------------------" << endl;
    cout << "n+1 = " << AddOne(a) << endl;
    cout << endl;

    return 0;
}

double AddOne ( double a )
{
    return a+1;
}