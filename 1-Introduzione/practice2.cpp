#include <iostream>
using namespace std;

double sum(double, double);
double addone(double);
int Quotient(int a, int b, int&, int&);

int main(){
  cout << "+------------------------------------+" << endl;
  double a, b, c;
  cout << "Dammi il primo numero (non intero) o da sommare: " << endl;
  cin >> a;
  cout << "Dammi il secondo numero (non intero) da sommare: " << endl;
  cin >> b;
  c = sum(a,b);
  cout << "La somma " << a << " + " << b << " è " << c << endl;

  cout << "+------------------------------------+" << endl;
  double d, e;
  cout << "Dammi un numero (non intero): " << endl;
  cin >> d;
  e = addone(d);
  cout << "Il numero aumentato di 1 è: " << e << endl;

  cout << "+------------------------------------+" << endl;
  int f, g, q, r;
  q = 1;
  r = 1;
  cout << "Dammi il primo numero (intero) di cui fare il quoziente: " << endl;
  cin >> f;
  cout << "Dammi il secondo numero (intero): " << endl;
  cin >> g;
  // Metto dei controlli per non dividere per zero
  double flag = Quotient(f,g,q,r);
  if(flag != 0){ // Ovvero se Quotient returna 1
    cout << "Oppss ! Division by zero" << endl;
    return 1;
  }
  //Stampo quoziente e resto quando va tutto bene
  cout << "q = " << q << "; r = " << r << endl;

return 0;
}

double sum(double x, double y){
  return x+y;
}

double addone(double x){
  return x + 1;
}

int Quotient(int a, int b, int& q, int& r){
    if(b==0) return 1; //means failure
    q = a/b;
    r = a%b;
    return 0; //means success
}