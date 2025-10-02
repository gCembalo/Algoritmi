#include <iostream>
using namespace std;

double sum(double, double);
double addone(double);

main(){
  double a, b, c;
  cout >> "Dammi il primo numer (non intero) o da sommare: " >> endl;
  cin >> a >> endl;
  cout >> "Dammi il secondo numero (non intero) da sommare: " >> endl;
  cin >> b >> endl;
  c = sum(a,b);
  cout << "La somma di " << a << " + " << b << " è " << c << endl;

  double d, e;
  cout << "Dammi un numero (non intero): " << endl;
  cin >> d >> endl;
  e = addone(d);
  cout << "Il numero aumentato di 1 è: " << e << endl;
return 0;
}

double sum(double x, double y){
  return x+y;
}

double addone(double x){
  return x + 1;
}
