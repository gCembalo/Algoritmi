#include <iostream> 

using namespace std;

int main ()             
{
    int i;

    cout << endl;
    
    while ( i<10 )
    {
        if ( i%2==0 )
        {
            cout << i+1 << endl;  
        }

        i++;
    }

    cout << endl;

    return 0;
}