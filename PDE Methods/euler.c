#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double l_init(double x);
double E_init(double x);
int a(int i, int j);
int I(int i);

void INIT_SYS(double * v, double * l, double * E, double h);
void EVAL_U(double * v, double * l, double * E, double * u);
void EVAL_P(double * v, double * l, double * E, double * p);
void EVAL_F(double * v, double * l, double * E, double * f, double * p);
void Flux_Calc(double * fluxP, double * fluxM, double * alpha, double * u, double * f);
void Update_U(double * u, double * unew);
void Update_Components(double * v, double * l, double * E, double * u);
void showArray(double * u);
void showComponentfromMatrix(FILE * pFILE, double * u, double * l);

const int N =100;
int row, col, count;
const double g = 1.4;

int main(void)
{
  FILE * pFILE;
  pFILE = fopen ("datei.dat","w");
  
  row=3;
  col=N;
  count = row * col; 
  double h = 1/(float)N; // x is [0,1]
  double T=0.2;
  
  double * v = (double *) malloc(col * sizeof(double));
  double * l = (double *) malloc(col * sizeof(double));
  double * E = (double *) malloc(col * sizeof(double));
  double * p = (double *) malloc(col * sizeof(double));
  double * c = (double *) malloc(col * sizeof(double));

  double * alpha = (double *) malloc(col * sizeof(double));
  double * fluxP  = (double *) malloc(row * col * sizeof(double));
  double * fluxM  = (double *) malloc(row * col * sizeof(double));

  double * u = (double *) malloc(row * col * sizeof(double));
  double * f = (double *) malloc(row * col * sizeof(double));
  double * unew = (double *) malloc(row * col * sizeof(double));

  //Main Procedures
  INIT_SYS(v, l, E, h);
  EVAL_U(v, l, E, u);
  EVAL_P(v, l, E, p);
  EVAL_F(v, l, E, f, p);
  showComponentfromMatrix(pFILE,u ,l); //wirte the first plot
  double t=0.0;
  int times =0; 
  // Time step loop starts here
  for(double j = 0.0; j < T; j+=t) {

    for(int i = 0; i < N; i++) 
      c[i] = sqrt(g * p[i] / l[i]); // calc c as in (8)

    for (int i=0; i < N; i++) {
      double C = 1.0;
      if(fabs(v[i+1])+c[i+1] < fabs(v[i])+c[i])
	alpha[i] = C * (fabs(v[i])+c[i]);
      else
	alpha[i] = C * (fabs(v[i+1])+c[i+1]);
    }

    Flux_Calc(fluxP, fluxM, alpha, u, f);
  
    double lamda = 0.0; // find lamda_max for time-step
    for(int i = 0; i < N; i++) 
      if(lamda < fabs(v[i])+c[i])
	lamda = fabs(v[i])+c[i];

    /* find a time step here in order to calc. new sol. */
    double CFL = 1.0;
    t = CFL * h / lamda;
  
    // find next solution
    for (int i=0; i < N; i++) {
      unew[a(0,i)] = u[a(0,i)] - (t/h)*(fluxP[a(0,i)] - fluxM[a(0,i)]);
      unew[a(1,i)] = u[a(1,i)] - (t/h)*(fluxP[a(1,i)] - fluxM[a(1,i)]);
      unew[a(2,i)] = u[a(2,i)] - (t/h)*(fluxP[a(2,i)] - fluxM[a(2,i)]);
    }

    Update_U(u, unew);
    Update_Components(v, l, E,u);
    EVAL_P(v, l, E, p);
    EVAL_F(v, l, E, f, p); 
    
    showComponentfromMatrix(pFILE, u, l);
    times += 1;
  }  
  printf("\n%d times",times);
  fclose (pFILE);
}

void Flux_Calc(double * fluxP, double * fluxM, double * alpha, double * u, double * f) {
  for (int i=0; i < col; i++) {

    fluxP[a(0,i)] = (0.5)*(f[a(0,I(i+1))]+f[a(0,I(i))]) - (0.5)*alpha[I(i+1)]*(u[a(0,I(i+1))]-u[a(0,I(i))]);
    fluxP[a(1,i)] = (0.5)*(f[a(1,I(i+1))]+f[a(1,I(i))]) - (0.5)*alpha[I(i+1)]*(u[a(1,I(i+1))]-u[a(1,I(i))]);
    fluxP[a(2,i)] = (0.5)*(f[a(2,I(i+1))]+f[a(2,I(i))]) - (0.5)*alpha[I(i+1)]*(u[a(2,I(i+1))]-u[a(2,I(i))]);
    
    fluxM[a(0,i)] = (0.5)*(f[a(0,I(i))]+f[a(0,I(i-1))]) - (0.5)*alpha[I(i)]*(u[a(0,I(i))]-u[a(0,I(i-1))]);
    fluxM[a(1,i)] = (0.5)*(f[a(1,I(i))]+f[a(1,I(i-1))]) - (0.5)*alpha[I(i)]*(u[a(1,I(i))]-u[a(1,I(i-1))]);
    fluxM[a(2,i)] = (0.5)*(f[a(2,I(i))]+f[a(2,I(i-1))]) - (0.5)*alpha[I(i)]*(u[a(2,I(i))]-u[a(2,I(i-1))]);

  } 
}

void EVAL_P(double * v, double * l, double * E, double * p){
  for (int i = 0; i < col; i++)
    p[i] = (g - 1) * (E[i]-(0.5)*l[i]*pow(v[i],2));
}
void EVAL_U(double * v, double * l, double * E, double * u){
  for (int i = 0; i < col; i++) {
    u[a(0,i)] = l[i];
    u[a(1,i)] = l[i] * v[i];
    u[a(2,i)] = E[i];
  }
}
void EVAL_F(double * v, double * l, double * E, double * f, double * p){
  for (int i=0; i<col; i++){
    f[a(0,i)] = l[i] * v[i];
    f[a(1,i)] = l[i]*pow(v[i],2) + p[i];
    f[a(2,i)] = v[i]*(E[i]+p[i]);
  }
}

void Update_U(double * u, double * unew){
  for (int i=0; i < col; i++) {
    u[a(0,i)] = unew[a(0,i)];
    u[a(1,i)] = unew[a(1,i)];
    u[a(2,i)] = unew[a(2,i)];
  }
}
void Update_Components(double * v, double * l, double * E, double * u){
  for (int i=0; i < col; i++) {
    l[i] = u[a(0,i)];
    v[i] = u[a(1,i)] / l[i];
    E[i] = u[a(2,i)];
  }
}

double l_init(double x){
  double init = 0.0;
  if(x<0.5)init=1.0;
  else init = 0.125;
  return init;
}
double E_init(double x){
  double init = 0.0;
  if(x<0.5)init=2.5;
  else init = 0.25;
  return init;
}
void INIT_SYS(double * v, double * l, double * E, double h){
  for (int i=0; i<col; i++){
    l[i]= l_init(i*h);
    v[i] = 0.0;
    E[i] = E_init(i*h);
  }
}

int I(int i){
  if (i < 0)
    i = N - abs(i);
  else if (i> N-1)
    i = i - N;
  return i;
}
int a(int i, int j){
  int a;
  a = i * col + j;
  return a; 
}

void showArray(double * u){
  for(int c =0; c<col; c++){
    printf("%.03lf\t",u[c]);
    if((c+1)%col==0)printf("\n");
  }
  printf("--------\n\n");
}
void showComponentfromMatrix(FILE * pFILE, double * u, double * l){
  double h = 1/(float)N; // x is [0,1]
  for(int c=0; c<col; c++){
    fprintf(pFILE,"%.03lf\t %.03lf\t %.03lf\t %.03lf\t\n",h*c,u[c],u[col+c]/l[c],u[2*col+c]);
    //if((c+1)%row==0)printf("\n");
  }
  fprintf(pFILE,"\n");
}
