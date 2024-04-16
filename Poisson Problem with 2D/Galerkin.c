#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const double e = 2.73;  //Need a finer definition
const int N = 100; // 100 or 1000
double f(double x){
  return pow(M_PI,2) * sin (M_PI*x);
}
double hat1(double x, double x1, double x2) {
  return (x - x1)/(x2 - x1);
}
double hat2(double x, double x1, double x2) {
  return (x2 - x)/(x2 - x1);
}
double Cal_hat1(double x1, double x2) {
  double xm;
  xm = (x2-x1)*0.5;
  return xm * (f(x2)*hat1(x2,x1,x2) + f(x1)*hat1(x1,x1,x2));
}
double Cal_hat2(double x1, double x2) {
  double xm;
  xm = (x2-x1)*0.5;
  return xm * (f(x2)*hat2(x2,x1,x2) + f(x1)*hat2(x1,x1,x2));
}
void Plot(FILE * pFILE, double * l, double * v, double * p);

int main(void) {
  int i;
  FILE * pFILE;
  pFILE = fopen ("datei.dat","w");
  /*construct a grid*/
  double * x = (double *) malloc((N+1) * sizeof(double));
  double alpha = 10.0; 
  for (i = 0; i<=N; i++) {
    x[i] = (pow(e,alpha*i/N)-1)/(pow(e,alpha)-1); 
    /* printf("x[%d]=%.07lf\n",i,x[i]); */
  }
  /* calculate h for each element */
  double * h = (double *) malloc(N * sizeof(double));
  for (i = 0; i<N; i++) {
    h[i] = x[i+1] - x[i];
    /* printf("h[%d]=%.07lf\n",i,h[i]); */
  }
  /* initialize a stiff matrix */
  double * m = (double *) malloc((N+1) * sizeof(double));
  double * s = (double *) malloc((N+1) * sizeof(double));
  double * F = (double *) malloc((N+1) * sizeof(double));
  double * U = (double *) malloc((N+1) * sizeof(double));
  for (i=0;i<=N;i++) {
    m[i]=0.0; s[i]=0.0; F[i]=0.0; U[i]=0.0;
  }
  m[0] = 1.0;    F[0]=0.0;
  m[N] = 1.0;    F[N]=0.0;
  m[1] = 1/h[0]; F[1] = Cal_hat1(x[0],x[1]);
  for(i=1;i<N-1;i++) {
    m[i] = m[i] + 1/h[i];
    s[i] = - 1/h[i];
    m[i+1] = 1/h[i];

    F[i]    = F[i] + Cal_hat2(x[i],x[i+1]);
    F[i+1]  = Cal_hat1(x[i],x[i+1]);
  }
  m[N-1] = m[N-1] + 1/h[N-1];
  F[N-1]     = F[N-1] + Cal_hat2(x[N-1],x[N]);

  // Tridiagonal solver
  for (i = 2; i < N; i++) {
    m[i] = m[i]-(s[i-1]*s[i-1]/m[i-1]);
    F[i] = F[i]-(s[i-1]*F[i-1]/m[i-1]);
  }
  // Backward substitution
  U[N-1] = F[N-1]/m[N-1];
  for (i = N-2; i > 0; i--) {
    U[i] = (F[i]-s[i]*U[i+1])/m[i];
  }
  // Error Calculation
  double hmax = 0.0;
  double ERR = 0.0;
  for (i=0;i<N;i++) {
    ERR += (x[i+1]-x[i])*0.5*(pow((U[i+1]-sin(M_PI*x[i+1])),2)+pow((U[i]-sin(M_PI*x[i])),2));
    if (hmax<=h[i])
      hmax = h[i];
  } 
  ERR = sqrt(ERR);
  printf("ERROR = %lf\t h_max = %lf\n\n",ERR,hmax);
  for (i = 0; i <= N; i++)
    printf("[%d] u:%lf sol:%lf \n",i,U[i],sin(M_PI*x[i]));

  for(i=0; i<=N; i++)
    fprintf(pFILE,"%lf\t%lf\n",x[i],U[i]);
  fprintf(pFILE,"\n"); // Between two set of data
  for(i=0; i<=N; i++)
    fprintf(pFILE,"%lf\t%lf\n",x[i],sin(M_PI*x[i]));
  
}














