#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int iter_max = 1000;

double sol_ref(double x, double y) 
{
   return sin(M_PI*x)*exp(-M_PI*y); 
}

void boundary(double **A, double *x, double *y, int n)
{
   /*Boundary conditions along the x axis for all processes */
      for ( int i=0; i<=n; i++)
      {
        A[i][0] = sol_ref(x[i], y[0]);
        A[i][n] = sol_ref(x[i], y[n]);
      }

    /*Boundary conditions along the y axis for all processes */
      for ( int j=0; j<=n; j++)
      {
        A[0][j] = sol_ref(x[0], y[j]);
        A[n][j] = sol_ref(x[n], y[j]);
      }
}


int main()
{   
  
  
  /* domain parameters */
  double a = 1.;       // Length of domain in x-direction
  double b = 1.;       // Length of domain in y-direction
  int n = 4096;         // Number of cells in x-direction
  int m = 4096;         // Number of cells in y-direction


  /* spatial discretization*/ 
  double dx = a/n;
  double dy = b/m;

 
  
  /* domain discretisatio */
  double *x;
  double *y;
  x=malloc((n+1)*sizeof(double));
  y=malloc((m+1)*sizeof(double)) ; 

  for ( int i=0; i<=n; i++)
  { 
     x[i] = i*dx;
     y[i] = i*dy; 
  } 
  
 /* solution variable */  
  double **A;
  double **Anew; 
 
 /*allocation*/
  A = malloc((n+1)*sizeof(double*));
  for (int i=0; i <= n; i++)
  { 
      A[i]=malloc((m+1)*sizeof(double));
  }
  
  Anew = malloc((n+1)*sizeof(double*));
  for (int i=0; i <= n; i++)
  {
      Anew[i]=malloc((m+1)*sizeof(double));
  }

 /* convergence parameters */  
  double tol = 10e-7;

  int iter=0;  
  double error;
  
  /*Boundary conditions */
  boundary(A, x, y, n);
  
  
  printf("------------%s----------\n", "Debut");            
  while (iter  < iter_max )
  {   
     
      error = 0.0;

      for( int j  = 1; j <= n-1; j++)
      {
         for( int i  = 1; i <= m-1; i++) 
         {  
            Anew[j][i] = 0.25 * ( A[j][i+1] + A[j][i-1] + A[j-1][i] + A[j+1][i] );
            error = fmax(error, fabs(Anew[j][i] - A[j][i]));  
            //error = error + pow( (A[j][i] - sol_ref(x[j],y[i])), 2);
           //printf("%f\n", Anew[j][i]);
         }
      } 
       
       /* Keep new temperature field */
       for( int j  = 1; j <= n-1; j++)
       {
         for( int i  = 1; i <= m-1; i++)
         {
            A[j][i] = Anew[j][i];
         }
       }
     
       if (iter % 100 == 0 ) printf("%d, %0.6f\n", iter, error);   
       
       if (error < tol) break;    

       iter++;
  }
  printf("----------%s----------\n", "Fin");

 /*delocate*/
  for (int i=0; i <= n; i++)
  {   
     free(Anew[i]);
     free(A[i]);
  }
  free(Anew);
  free(A);
  free(x);
  free(y);

  return 0;
  
}
