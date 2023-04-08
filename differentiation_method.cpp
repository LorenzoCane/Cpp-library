// differentiation method implementations FD, BD, CD


#include "differentiation_method.h"



double ForwardDifference(double (*F)(double a), double x, double h){

 // Forward difference using h as increment, store value of f' in f prime

    double fPrime = ( F(x+h) - F(x)) / h ;

    return fPrime;
}

//-----------------------------------------------------------------------------------------

double BackwardDifference (double (*F)(double x), double x, double h ){

 // Backward difference using h as increment, store value of f' in f prime

    double fPrime = ( F(x) - F(x-h)) / h ;

    return fPrime;
}

//-----------------------------------------------------------------------------------------

double CentralDifference (double (*F)(double x), double x, double h ){

 // Central difference using h as increment, store value of f' in f prime

    double fPrime = 0.5 * ( F(x+h) - F(x-h)) / h ;

    return fPrime;
}

//---------------------------------------------------------------------------------------

double FourthDifference(double (*F)(double x), double x, double h){

 // Fourth order difference using h as increment, store value of f' in f prime

    double fPrime = ( F(x- 2*h) - 8 * F(x-h) + 8 * F(x+h) - F(x +2*h)) / (12*h)  ;

    return fPrime;
}

//---------------------------------------------------------------------------------------

double SecDifference (double (*F)(double x), double a, double h){  

 // Second derivate for F function use it only in centered function

    double fSecond = (F(a+h) - 2* F(a) + F(a-h))/ (h*h);

    return fSecond;
}


//---------------------------------------------------------------------------------------

double OneSideSecDifference (double (*F)(double x), double a, double h ){

 // Second derivate for function F, it can be use one patological function with only a right side

    double fSecond = (F(a+h+h) + F(a) - 2 * F(a+h)) / (h*h) ;

    return fSecond;
}

