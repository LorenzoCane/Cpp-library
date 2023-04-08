/* Header file with solution methods for elliptic pde using Dirichlet b. c. 
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

void Jacobi_Method(double **phi_old, double *x, double *y, void (*BoundaryConditions)(double **,double *, double *, int, int), double tol, int* iterations);