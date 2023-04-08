// header file for chapter 03
// include all the quadrature methods definitions, implementations can be find in .cpp 


#ifndef  QUADRATURE_H
#define  QUADRATURE_H

#define  _USE_MATH_DEFINES    

#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>

using namespace std;

double RectangularRule(double (*F) (double), double , double, int);
double MidPointRule(double (*F) (double), double , double, int);
double TrapezoidalRule(double (*F) (double), double , double, int);
double SimpsonRule(double (*F) (double), double , double, int);
double GaussQuadratureRule (double (*F) (double) , double , double , int , int );

#endif