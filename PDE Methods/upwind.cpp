#include <stdlib.h>
#include <iostream>
#include <math.h>
using namespace std;

void Upwind(int N, double C, double *uold, double *unew);
void LW(int N, double C, double *uold, double *unew);
double func(double x, char mode);
double exact(double x, double); 
int ind(int i,int N);

int main(int argc, char *argv[])
{
  int N;
  char S;
  char mode;
  cout << "Input N:  " ;
  cin >> N ;
  while(1){
    cout << "\n[U]pwind or [L]ax-Wendroff : ";
    cin >> S;
    if(S=='U'||S=='u'||S=='L'||S=='l') break;
  }
  cout << "[S]moooth or [D]iscontinuous :";
  cin >> mode;
  
  int xmin = 0;
  int xmax = 1;
  double h = (xmax - xmin)/(double)N;

  double a = 1;
  double v = 0.9; 
  int M = (int)(N/v);  

  double dt = v*h/a;
  double T = M*dt;
  
  double C = v/2; // times another v in LW function
  double * uold = new double[N+1];
  double * unew = new double[N+1];
  double * ERR = new double[M];

  // Initialize our first 'uold'
  for (int i = 0; i < N; i++){
    uold[i] = func(i*h,mode); // to do each step(*dx)
  }

  for (int i = 0; i < M ; i++){
   
    if(S=='U'||S=='u')
      Upwind(N,C,uold,unew);
    else if(S=='L'||S=='l')
      LW(N,C,uold,unew);
    
    //    cout << "\n----- Round [" << i << "] -----" << endl;
    for (int j = 0; j < N; j++) uold[j] = unew[j];
    unew[N] = unew[0];     // Set the boundary conditions
   
  }

  if(S=='U'||S=='u')
    cout << "\nUpwind\n";
  else if(S=='L'||S=='l')
    cout << "\nLax-Wendroff\n";
  for (int k = 0; k <= N; k++)  cout << unew[k] << endl;
  cout << "\n Exact Sol" << endl;
  for (int k = 0; k <= N; k++)  cout<< exact(k*h,T) << endl;
  if(mode=='S'||mode=='s') {
    double e=0.0;
    for(int i=0; i< N; i++) {
      if( fabs( unew[i] - exact(i*h,T) ) > e )
        e = fabs( unew[i] - exact(i*h,T)); }
    cout << "\n----- ERROR -----\n";
    cout << e << endl;
  }


  return 0;
}

void Upwind(int N, double C, double *uold, double *unew){
  for (int i =0; i< N; i++)
    unew[i] =
        uold[ind(i,N)]-C*(uold[ind(i+1,N)]-uold[ind(i-1,N)])
        +C*(uold[ind(i+1,N)]-2*uold[ind(i,N)]+uold[ind(i-1,N)]);
}

void LW(int N, double C, double *uold, double *unew){
  for (int i =0; i< N; i++)
    unew[i] =
        uold[ind(i,N)]-C*(uold[ind(i+1,N)]-uold[ind(i-1,N)])
        +C*0.9*(uold[ind(i+1,N)]-2*uold[ind(i,N)]+uold[ind(i-1,N)]);
}

double func(double x, char mode){
  double U;
  switch(mode){
    case 's':
      U = sin(2.0*M_PI*x);
      break;
    case 'd':
      if(x>=0.0 && x <= 0.5)
        U=1.0;
      else if(x>=0.5 && x<1.0)
        U=0.0;
      break;
  }
  return U;
}

double exact(double x, double t){
  double Uex;
  Uex = sin(2.0*M_PI*(x-t));
  return Uex;
}

int ind(int i, int N){
  if (i < 0)
    i = N - abs(i);
  else if (i> N-1)
    i = i - N;
  return i;
}
