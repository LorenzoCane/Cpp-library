// header files with differentiation methods. It includes FD,BD,CD


#ifndef DIFFERENTITION_METHOD_H

#define DIFFERENTITION_METHOD_H

#define _USE_MATH_DEFINE

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;


double ForwardDifference (double (*F)(double x), double x, double h );

double BackwardDifference (double (*F)(double x), double x, double h);

double CentralDifference (double (*F)(double x), double x, double h);

double FourthDifference (double (*F)(double x), double x, double h );

double SecDifference (double (*F)(double x), double a, double h );

double OneSideSecDifference (double (*F)(double x), double a, double h );

#endif