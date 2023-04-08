/* this header includes all the ode resolver definitions 
  implementations can be found in .cpp file
*/


#ifndef ODE_METHODS_H
#define ODE_METHODS_H

#define _USE_MATH_DEFINE

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;
 
void EulerStep      (double t , double *Y , void (*RHSFunc)(double , double * , double *) , double dt , int neq);  
void RK2Step        (double t , double *Y , void (*RHSFunc)(double , double * , double *) , double dt , int neq);
void RK4Step        (double t , double *Y , void (*RHSFunc)(double , double * , double *) , double dt , int neq); 
void PositionVerlet (double t , double *x , double *v , void(*acc)(double , double * , double *) , double dt , int npart);

#endif