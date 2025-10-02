#include <iostream>
using namespace std;

int main(){

double xend=10;
double dx = 10/xend;
double x = 0.0;
int n = 0;
cout << endl;
cout << "+-------------------------+" << endl;
cout << "Ciclo con il float:" << endl;
for(double x; x<=xend; x+=dx){
  cout << n << ": " << x << endl;
  n+=1;
}
cout << "+-------------------------+" << endl;


int N=10;
double dy = 10/N;
double y = 0.0;
int m = 0;
cout << endl;
cout << "+------------------------+" << endl;
cout << "Ciclo con l'intero:" << endl;
for(int i=0; i<=N; i++){
  y+=dy;
  cout << i << ": " << y << endl;
}
cout << "+------------------------+" << endl;
cout << endl;

return 0;
}
