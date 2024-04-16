//
//  main.c
//  PDE1
// 
//  Created by WM on 02/11/15.
//  Copyright © 2015 WM. All rights reserved.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void dgesv_(const int *N, const int *nrhs, double *A, const int *lda, int
	    *ipiv, double *b, const int *ldb, int *info);

void dgels_(const char *trans, const int *M, const int *N, const int *nrhs,
	    double *A, const int *lda, double *b, const int *ldb, double *work,
	    const int * lwork, int *info);

unsigned factorial(unsigned n);


int main(void) {
  int m,p,q;
  printf("\n\n[Start]\nInsert the order of derivative:");
  scanf("%d",&m);
  printf("Insert stencils p and q:");
  scanf("%d%d",&p,&q);
  p = abs(p);
  q = abs(q);
  while(p+q+1<m+2) {
    printf("Please insert p and q again:");
    scanf("%d%d",&p,&q);
    p = abs(p);
    q = abs(q);
  }
    
  printf("Order:%d p:%d q:%d\n\n",m,p,q); //End of user's input
    
  //adjust the matrix size
  int COL;
  COL=p+q+1;
    
    
  //Maybe define Malloc here to reserve the memory
    
  //Taylor series expansion
  float M[COL][COL],Mt[COL][COL];
  double A[COL*COL];
  int r=0; // row of the matrix
  int h=0,o=0;
  for (h=(-p); h<=q; h++)
    {
      printf("U(x+(%dh)) = ",h);
      for (o=0; o<COL; o++)
	{ //also use 'o' as a col.
	  M[r][o]= pow(h, o)/factorial(o);
	  printf("(%.3f)U_%d  ",M[r][o],o);
	  if (o==COL-1) printf("\n");
	}
      r++;
    }
    
  //Transpose matrix A
  for (int i=0; i<COL; i++)
    {
      for (int j=0; j<COL; j++)
	{
	  Mt[i][j]=M[j][i];
	}
    }
    
  //store the matrix in Fortran format
  int c=0;
  for (int j=0; j<COL; j++)
    {
      for (int i=0; i<COL; i++)
	{
	  A[c]=Mt[i][j];
	  c++;
	}
    }
    
  //check the values
  printf("\n\nA = [");
  for (int i=0;i<(COL*COL); i++)
    {
      printf("%.3f  ",A[i]);
      if ((i+1)%COL==0) printf(",");
    }
  printf("]\n");
    
    
  printf("\nRead matrix A\n");
  //Read out our matrix
  for (int i=0; i<COL; i++)
    for (int j=0; j<COL;j++)
      { printf("%.3f   ",Mt[i][j]);
	if (j==COL-1) printf("\n");}
    
  //Build matrix B
  double B[COL];
  for (int i=0; i<COL; i++) //set all b to 0.
    {B[i]=0;}
  B[m]=1; //set the derivative
    
  printf("\nRead matrix B\n");
  for (int i=0; i<COL; i++) //check values in B
    printf("%.3f\n",B[i]);
    
    
  /*=================================================*/
  //Linear algebra solver
    
  int N = p+q+1;
  int nrhs = 1;
  int lda = p+q+1;
  int ipiv[p+q+1];
  int ldb = p+q+1;
  int info;
    
  dgesv_(&N, &nrhs, A, &lda, ipiv, B, &ldb, &info);
    
  if(info == 0) /* succeed */
    {       printf("\nSolution for #%d derivative\n",m);
      for (int i=0; i<COL; i++)
	{printf("%lf \n", B[i]);}
    }
  else
    fprintf(stderr, "dgesv_ fails %d\n", info);

  printf("\nAcc. of Apprx. = %d\n\n",COL-m+1); //Print out the ACC
    
}//end

unsigned factorial(unsigned n)
{
  if (n == 1)
    return 1;
  else if (n == 0)
    return 1;
  else
    return n * factorial(n - 1);
}
