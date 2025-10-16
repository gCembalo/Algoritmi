#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int main(){
  cout << setiosflags ( ios::scientific );
  int N=10; // numero di step

  double xend=10;
  double dx = xend/N;
  double xbeg = 0.0;
  double x;
  int n = 0;
  cout << endl;
  cout << "+-------------------------+" << endl;
  cout << "Ciclo con il float:" << endl;
  for(double z=xbeg; z<xend; z+=dx){
    cout << n << ": " << x << endl;
    n+=1;
    x = z;
  }
  cout << "Visto che x = " << x << "è il metodo sbagliato." << endl;
  cout << "+-------------------------+" << endl;

  double dy = 10/N;
  double y = 0.0;
  int m = 0;
  cout << endl;
  cout << "+------------------------+" << endl;
  cout << "Ciclo con l'intero:" << endl;
  for(int i=0; i<N; i++){
    y+=dx;
    cout << i << ": " << y << endl;
  }
  cout << "Visto che y = " << y << "è il metodo corretto." << endl;
  cout << "+------------------------+" << endl;
  cout << endl;

return 0;
}
