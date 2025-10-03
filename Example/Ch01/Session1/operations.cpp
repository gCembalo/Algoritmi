#include <iostream>

using namespace std;

int main ()             
{
    int a, b;

    cout << endl;
    cout << "Inserire due numeri interi" << endl;

    cin >> a >> b;

    cout << "-------------------" << endl;
    cout << "somma: \t\t\t" << a+b << endl;
    cout << "differenza: \t\t" << a-b << endl;
    cout << "prodotto: \t\t" << a*b << endl;

    if ( b==0 )
    {
        cout << "Non si può dividere per zero!" << endl;
    }
    else
    {
        cout << "quoziente intero: \t" << a/b << endl;
        cout << "quoziente reale: \t" << (float)a/(float)b << endl;       
    }

    float c, d;

    cout << endl;
    cout << "Inserire due numeri reali" << endl;

    cin >> c >> d;
    
    cout << "-------------------" << endl;
    cout << "somma: \t\t\t" << c+d << endl;
    cout << "differenza: \t\t" << c-d << endl;
    cout << "prodotto: \t\t" << c*d << endl;

    if ( d==0 )
    {
        cout << "Non si può dividere per zero!" << endl;
    }
    else
    {
        cout << "quoziente intero: \t" << (int)(c/d) << endl;
        cout << "quoziente reale: \t" << c/d << endl;       
    }

    cout << endl;

    return 0;
}