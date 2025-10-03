#include <iostream> 

using namespace std;

void Quoziente ( int, int, int&, int& );

int main ()             
{
    int a, b;
    int q = 0, r = 0;

    cout << endl;
    cout << "Inserire due numeri interi" << endl;

    cin >> a >> b;

    cout << "--------------------------" <<  endl;

    if ( b==0 )
    {
        cout << "Non si può dividere per zero!" << endl;
    }
    else
    {
        Quoziente ( a, b, q, r );

        cout << a << "/" << b << " = " << q << " resto " << r << endl;
    }
    
    cout << endl;

    return 0;
}

void Quoziente ( int a, int b, int& q, int& r )
{
    q = a/b;
    r = a%b;
}