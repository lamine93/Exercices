#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// problem parameters
const double a = 1.;
const double b = 1.;

const int nx = 1024; //number of node points along y
const int ny = 1024; //number of node points along x

//convergence parameters
double tol = 1e-4;
int iter_max = 1000;

double sol_ref(double x, double y)
{
   return sin(M_PI*x)*exp(-M_PI*y);
}

void discretisation(double *x, double *y)
{
  double dx = a/nx; // spatial step along x
  double dy = b/ny; // spatial step along y

  for (int i=0; i<=nx; i++)
  {
      x[i] = i*dx;
      y[i] = i*dy;
  }

}

void boundary(double *T, double *x, double *y)
{
      /*Boundary conditions along the x axis for all processes */
      for ( int i=0; i<=nx; i++)
      {
        T[i*nx] = sol_ref(x[i], y[0]);
        T[i*ny+ny] = sol_ref(x[i], y[ny]);
      }

      /*Boundary conditions along the y axis for all processes */
      for ( int j=0; j<=ny; j++)
      {
        T[j] = sol_ref(x[0], y[j]);
        T[nx*ny+j] = sol_ref(x[nx], y[j]);
      }

}

void laplace2d(double *T, double *Tnew, double *error)
{

   for( int j  = 1; j <= nx-1; j++)
   {
     for( int i  = 1; i <= ny-1; i++)
     {
        Tnew[j*nx + i] = 0.25 * ( T[j*nx + (i+1)] + T[j*nx + (i-1)] + T[(j-1)*nx + i] + T[(j+1)*nx + i] );
        *error = fmax(*error, fabs(Tnew[j*nx + i] - T[j*nx + i]));
     }
   }

   for( int j = 1; j <= nx-1; j++)
   {
     for( int i  = 1; i <= ny-1; i++)
     {
        T[j*nx + i] = Tnew[j*nx + i];
     }
   }

}

int main(int argc,char **argv)
{
  double tol = 1e-4;
  int iter = 0;
  double error;
  
  double *T    = (double*) malloc(sizeof(double) * (nx+1) * (ny+1));
  double *Tnew = (double*) malloc(sizeof(double) * (nx+1) * (ny+1));
   
  double *x =  (double*) malloc(nx * sizeof(double));
  double *y =  (double*) malloc(ny * sizeof(double));
  
  if(!x || !y || !T || !Tnew ) return 0;
 
  discretisation(x, y);

  boundary(T, x, y);
    
  while (iter  < iter_max )
  { 
    error = 0.0;
    
    laplace2d(T, Tnew, &error);
    
    if (iter % 100 == 0 ) printf("%d, %0.6f\n", iter, error);

    if (error < tol) break;

    iter++;
  }

  free(T);
  free(Tnew);
  free(x);
  free(y);
  
  return 0;
}
