#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double l_init(double x);
double E_init(double x);
int a(int i, int j);
int m(int i, int j); // for matrix mul 3 x 3

void INIT_SYS(double * v, double * l, double * E, double h);
void EVAL_U(double * v, double * l, double * E, double * u);
void EVAL_P(double * v, double * l, double * E, double * p);
void EVAL_F(double * v, double * l, double * E, double * f, double * p);
void Update_U(double * u, double * unew);
void Update_Components(double * v, double * l, double * E, double * u);
void showArray(double * u);
void showMatrix(double * u, int s);
void Plot(FILE * pFILE, double * l, double * v, double * p);

const int N = 1000;
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
  double * fluxP  = (double *) malloc(row * col * sizeof(double));
  /* double * fluxM  = (double *) malloc(row * col * sizeof(double)); */

  double * u = (double *) malloc(row * col * sizeof(double));
  double * f = (double *) malloc(row * col * sizeof(double));
  double * unew = (double *) malloc(row * col * sizeof(double));

  // Aroe computation
  double * enthalpy = (double *) malloc(col * sizeof(double));
  double * uavg = (double *) malloc(col * sizeof(double));
  double * havg = (double *) malloc(col * sizeof(double));
  double * cavg = (double *) malloc(col * sizeof(double));
  double * Tm = (double *) malloc(col * sizeof(double));
  // Matrix Multiplication
  double * R = (double *) malloc(3 * 3 * sizeof(double));
  double * L = (double *) malloc(3 * 3 * sizeof(double));
  double * InvR = (double *) malloc(3 * 3 * sizeof(double));
  double * tmpMatrix = (double *) malloc(3 * 3 * sizeof(double));
  double * tmpMatrix2 = (double *) malloc(3 * 3 * sizeof(double));
  
  double * udiff = (double *) malloc(row * sizeof(double));
  double * A = (double *) malloc(row * sizeof(double));
  //Main Procedures
  INIT_SYS(v, l, E, h);
  EVAL_U(v, l, E, u);
  EVAL_P(v, l, E, p);
  EVAL_F(v, l, E, f, p);
  /* Plot(pFILE,l ,v, p); //wirte the first plot */

  double t = 0.0;
  int times = 0; 
  // Time step loop starts here
  for(double j = 0.0; j < T; j+=t) {
    
    for (int i = 0; i < col; i++){
      enthalpy[i] = 0.0; uavg[i] = 0.0;
      havg[i] = 0.0; cavg[i] = 0.0; Tm[i]   = 0.0;
    }
    
    for (int i = 0; i < col; i++)
      enthalpy[i] = (E[i] + p[i]) / l[i];
  
    for (int i = 0; i < col-1; i++) {
      uavg[i] = (sqrt(l[i])*v[i] + sqrt(l[i+1])*v[i+1]) / (sqrt(l[i]) + sqrt(l[i+1]));
      havg[i] = (sqrt(l[i])*enthalpy[i] + sqrt(l[i+1])*enthalpy[i+1]) / (sqrt(l[i]) + sqrt(l[i+1]));
      cavg[i] = sqrt((g-1)*(havg[i]-(0.5)*pow(uavg[i],2)));
      Tm[i]   = (g-1)/pow(cavg[i],2);
    }

    // SUPER Flux calculation starts here
    for (int r = 0; r < col-1; r++) {
      
      for(int i=0; i<3; i++)
	for(int j=0; j<3; j++) {
	  R[m(i,j)] = 0.0;
	  L[m(i,j)] = 0.0;
	  InvR[m(i,j)] = 0.0;
	}
      
      // initialize matrices
      R[m(0,0)] = 1.0; R[m(0,1)] = 1.0; R[m(0,2)] = 1.0;
      R[m(1,0)] = uavg[r]-cavg[r];
      R[m(1,1)] = uavg[r]; 
      R[m(1,2)] = uavg[r]+cavg[r];
      R[m(2,0)] = havg[r]-uavg[r]*cavg[r]; 
      R[m(2,1)] = (0.5)*pow(uavg[r],2);
      R[m(2,2)] = havg[r]+uavg[r]*cavg[r];
      
      L[m(0,1)] = 0.0; L[m(0,2)] = 0.0; L[m(1,0)] = 0.0;
      L[m(1,2)] = 0.0; L[m(2,0)] = 0.0; L[m(2,1)] = 0.0;

      if (uavg[r]-cavg[r]< cavg[r])
	L[m(0,0)] = fabs((0.5) * (cavg[r]+pow(uavg[r]-cavg[r],2)/cavg[r]));
      else
	L[m(0,0)] = fabs(uavg[r]-cavg[r]);

      if (uavg[r]< cavg[r])
	L[m(1,1)] = fabs((0.5) * (cavg[r]+pow(uavg[r],2)/cavg[r]));
      else
	L[m(1,1)] = fabs(uavg[r]);
      
      if (uavg[r]+cavg[r] < cavg[r])
	L[m(2,2)] = fabs((0.5) * (cavg[r]+pow(uavg[r]+cavg[r],2)/cavg[r]));
      else
	L[m(2,2)] = fabs(uavg[r]+cavg[r]);
      
      InvR[m(0,0)] = (0.5)*(1+Tm[r]*(pow(uavg[r],2)-havg[r])+uavg[r]/cavg[r]); 
      InvR[m(0,1)] = (-0.5)*(Tm[r]*uavg[r]+1/cavg[r]); 
      InvR[m(0,2)] = (0.5)*Tm[r];
      InvR[m(1,0)] = -Tm[r]*(pow(uavg[r],2)-havg[r]); 
      InvR[m(1,1)] = Tm[r]*uavg[r]; 
      InvR[m(1,2)] = -Tm[r];
      InvR[m(2,0)] = (0.5)*(1+Tm[r]*(pow(uavg[r],2)-havg[r])-uavg[r]/cavg[r]); 
      InvR[m(2,1)] = (-0.5)*(Tm[r]*uavg[r]-1/cavg[r]); 
      InvR[m(2,2)] = (0.5)*Tm[r];
    
      // initialize matrices
      for(int i=0; i<3; i++)
	for(int j=0; j<3; j++) {
	  tmpMatrix[m(i,j)] = 0.0;
	  tmpMatrix2[m(i,j)] = 0.0;
	}
      // first multiplication
      for(int i=0; i<3; i++)  // maybe it should be i++
	for(int j=0; j<3; j++)
	  for(int k=0; k<3; k++)
	    tmpMatrix[m(i,j)]+=L[m(i,k)]*InvR[m(k,j)];
      // second multiplication 
      for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
	  for(int k=0; k<3; k++)
	    tmpMatrix2[m(i,j)]+=R[m(i,k)]*tmpMatrix[m(k,j)];

      for(int i=0; i<3; i++)
	udiff[i] = A[i] = 0.0;
     
      udiff[0] = u[a(0,r+1)]-u[a(0,r)];
      udiff[1] = u[a(1,r+1)]-u[a(1,r)];
      udiff[2] = u[a(2,r+1)]-u[a(2,r)]; 

      // one last matrix multiplication
      for(int i=0; i<3; ++i)
	for(int j=0; j<1; ++j)
	  for(int k=0; k<3; ++k)
	    A[i]+=tmpMatrix2[m(i,k)]*udiff[k];
      
      fluxP[a(0,r)] = (0.5)*(f[a(0,r+1)]+f[a(0,r)]) - (0.5)*A[0];
      fluxP[a(1,r)] = (0.5)*(f[a(1,r+1)]+f[a(1,r)]) - (0.5)*A[1];
      fluxP[a(2,r)] = (0.5)*(f[a(2,r+1)]+f[a(2,r)]) - (0.5)*A[2];
    }
 
    /* find a time step here in order to calc. new sol. */
    double lamda = 0.0; // find lamda_max for time-step
    for(int i = 0; i < N-1; i++) 
      if(lamda < fabs(uavg[i])+cavg[i])
	lamda = fabs(uavg[i])+cavg[i];
      
    double CFL = 0.9;
    t = CFL * h / lamda;
  
    // Boundary condition
    unew[a(0,0)] = u[a(0,0)];
    unew[a(1,0)] = u[a(1,0)];
    unew[a(2,0)] = u[a(2,0)];
    unew[a(0,N-1)] = u[a(0,N-1)];
    unew[a(1,N-1)] = u[a(1,N-1)];
    unew[a(2,N-1)] = u[a(2,N-1)];
    // find next solution
    for (int i=1; i < col-1; i++) {
      unew[a(0,i)] = u[a(0,i)] - (t/h)*(fluxP[a(0,i)] - fluxP[a(0,i-1)]);
      unew[a(1,i)] = u[a(1,i)] - (t/h)*(fluxP[a(1,i)] - fluxP[a(1,i-1)]);
      unew[a(2,i)] = u[a(2,i)] - (t/h)*(fluxP[a(2,i)] - fluxP[a(2,i-1)]);
    }
    
    Update_U(u, unew);
    Update_Components(v, l, E,u);
    EVAL_P(v, l, E, p);
    EVAL_F(v, l, E, f, p); 
    
    /* Plot(pFILE, l, v, p); */
    times += 1;
  }  
  printf("\n%d times",times);
  
  Plot(pFILE, l, v, p);
  fclose (pFILE);
  // FREE
  free(v); free(l); free(E); free(p); free(c); free(fluxP);
  free(u); free(f); free(unew);
  free(enthalpy); free(uavg); free(havg); free(cavg); free(Tm);
  free(R); free(L); free(InvR); free(tmpMatrix); free(tmpMatrix2);
  free(udiff); free(A);
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

int a(int i, int j){
  int a;
  a = i * col + j;
  return a; 
}
int m(int i, int j){
  int m;
  m = i * 3 + j;
  return m; 
}

void showMatrix(double * u,int s){
  for(int c =0; c<s; c++){
    printf("%.03lf\t",u[c]);
    if((c+1)%9==0)printf("\n");
  }
  printf("--------\n\n");
}
void showArray(double * u){
  for(int c =0; c<col; c++){
    printf("%.03lf\t",u[c]);
    if((c+1)%col==0)printf("\n");
  }
  printf("--------\n\n");
}
void Plot(FILE * pFILE, double * l, double * v, double * p){
  double h = 1/(float)N; // x is [0,1]
  for(int c=0; c<col; c++){
    fprintf(pFILE,"%.03lf\t %.03lf\t %.03lf\t %.03lf\t\n",h*c,l[c],v[c],p[c]);
    //if((c+1)%row==0)printf("\n");
  }
  fprintf(pFILE,"\n");
}
