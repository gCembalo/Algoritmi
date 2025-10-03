#include <iostream> 

using namespace std;

float Sum ( float, float );

int main ()             
{
    float a, b;

    cout << endl;
    cout << "Inserire due numeri reali" << endl;

    cin >> a >> b;
    
    cout << "-------------------------" << endl;
    cout << "somma: " << Sum(a,b) << endl;
    cout << endl;

    return 0;
}

float Sum ( float a, float b )
{
    return a+b;
}