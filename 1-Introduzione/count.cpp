#include <iostream>
using namespace std;

int main()
{
cout << endl;
cout << "+-----------------------+" << endl;
cout << "Con il ciclo for:" << endl;
  for(int i=1; i<11; i++){
    cout << i << endl;
  }
cout << endl;
cout << "+-----------------------+" << endl;
cout << "Con il ciclo while:" << endl;
int j=1;
  while(j<11){
    cout << j << endl;
    j++;
  }
cout << endl;
cout << "+----------------------+" << endl;
cout << "Stampo solo i numeri dispari:" << endl;
  for(int i=1; i<11; i++){
    if(i % 2 == 1){
      cout << i << endl;
    }
}
return 0;
}
