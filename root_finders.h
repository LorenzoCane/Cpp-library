// header file for chapter 05
// include all the root finders definitions, implementations can be find in .cpp 


#ifndef  ROOT_FINDERS_H
#define  ROOT_FINDERS_H

#define  _USE_MATH_DEFINES    

#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>

using namespace std;

int Bisection ( double (*F) (double) , double , double , double , double , double& , int , int&);
int FalsePosition ( double (*F) (double) , double , double , double ,double , double&, int ,int&);
int Secant ( double (*F)(double), double , double , double , double& , int, int&);
int Newton(double (*F)(double) , double (*Fprime)(double) ,double a, double b,double tol, double &, int nmax, int&);

void Bracketing(double (*F) (double), double a, double b , int N ,double* xL , double*  xR , int& nzero);

#endif