#include <stdio.h>

// problem parameters
const int a = 1
const int b = 1

const int nx = 512 //number of node points along y
const int ny = 512 //number of node points along x

//convergence parameters
double tol = 10e-7;
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
        T[i*n] = sol_ref(x[i], y[0]);
        T[i*n+n] = sol_ref(x[i], y[ny]);
      }

      /*Boundary conditions along the y axis for all processes */
      for ( int j=0; j<=ny; j++)
      {
        T[j] = sol_ref(x[0], y[j]);
        T[n*n+j] = sol_ref(x[nx], y[j]);
      }

}

void laplace2d(double *T, double *Tnew)
{

  while (iter  < iter_max )
  {

        error = 0.0;

        for( int j  = 1; j <= nx-1; j++)
        {
           for( int i  = 1; i <= ny-1; i++)
           {
              Tnew[j*nx + i] = 0.25 * ( T[j*nx + (i+1)] + T[j*nx + (i-1)] + T[(j-1)*nx + i] + T[(j+1)*nx + i] );
              error = Tmax(error, fabs(Tnew[j*nx + i] - T[j*nx + i]));
           }
        }

        for( int j = 1; j <= nx-1; j++)
        {
           for( int i  = 1; i <= ny-1; i++)
           {
              T[j*nx + i] = Tnew[j*nx + i];
           }
         }

         if (iter % 100 == 0 ) printf("%d, %0.6f\n", iter, error);

         if (error < tol) break;

         iter++;
   }  

}

int main(int argc,char **argv)
{
  double *T    = (double*) malloc(sizeof(double) * nx * ny);
  double *Tnew = (double*) malloc(sizeof(double) * nx * ny);
   
  double *x = malloc((nx+1)*sizeof(double));
  double *y = malloc((ny+1)*sizeof(double));

  discretisation(x, y);

  boundary(T, x, y);
  
  laplace2d(T, Tnew);

  free(T);
  free(Tnew);
  free(x);
  free(y);
  
  return 0;
}
