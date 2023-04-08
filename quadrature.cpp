// include implementations of quadrature methods functions 

#include "quadrature.h"

double RectangularRule(double (*F)(double x), double a , double b , int n){

    double h= ( b - a)/(double)n ;
    double sum = 0.0;

    for (int i=0 ; i<n ; i++){

        sum = sum + F(a+i*h);     // add the n values of the interval

    }

    sum *= h ;    // *h at the end to save time

    return sum;

}



//----------------------------------------------------------------------------------------------
double MidPointRule( double (*F)(double x), double a , double b , int n){

    double h= ( b - a)/(double)n ;
    double sum = 0.0;

    for(int i= 1 ; i<n ; i++ ) {   //mid points start in a+h/2 and end in a+ (n-1) h/2

        sum = sum + F(a + 0.5*(2*i-1)*h);  //only 2n-1 * i
    }

    sum *= h ;

    return sum;
}

//---------------------------------------------------------------------------------------------

double TrapezoidalRule( double (*F)(double x) , double a , double b , int n) {

    double h= ( b - a)/(double)n ;
    double sum = 0.5*(F(a)+F(b));        //add initial and final values before

    for(int i = 1 ; i<n ; i++){

        sum = sum + F(a+i*h);
    }

    sum *= h;

    return sum;
}

//---------------------------------------------------------------------------------------------

double SimpsonRule( double (*F)(double x) , double a , double b , int n) {

    double h= ( b - a)/(double)n ;
    double sum = (F(a)+F(b))/3;        //add initial and final contributes before
    double w = 4.0;                    // this parameter is used to avoid if else structure

    for(int i = 1 ; i<n ; i++){

        sum = sum + w*F(a+i*h)/3;    // use w = 4.0 and do w= 6.0- w for each cycle
        w = 6.0 - w ;                // if w = 4 the next one will have w = 2 and viceversa --> if you have to move between two value use this strategy
    }

    sum *= h;

    return sum;
}

//---------------------------------------------------------------------------------------------

double GaussQuadratureRule(double (*F)(double) , double a , double b , int N ,  int Ng){

    double x[10];  // gaussian points
    double w[10];  // weights
    double sum = 0.0; // global sum
    double h = fabs(b - a)/N ;


    if (Ng==2) {
        x[0] = 1.0/sqrt(3.0);
        x[1] = - x[0];
        w[0] = 1.0 ;
        w[1] = 1.0;
    }

    else if (Ng == 3){
        x[0] = 0.0;
        x[1] = sqrt( 3.0/5.0);
        x[2] = -x[1];
        w[0] = 8.0/9.0 ;
        w[2] = 5.0/9.0 ;  
        w[1] = w[2] ;     
    }

    else if (Ng == 4){
        x[0] = sqrt(3.0/7.0 - 2.0/7.0 * sqrt(0.2*6.0)) ;
        x[1] = - x[0] ;
        x[2] = sqrt (3.0/7.0 + 2.0/7.0*sqrt(0.2*6.0)) ;
        x[3] = - x[2] ; 
        w[0] = w[1] = (18.0 + sqrt(30.0))/36.0 ;
        w[2] =  w[3] = (18.0 - sqrt(30.0))/36.0 ; 
    }

    for( int i =0 ; i < N ; i++){
        double x0 = a + i*h ;
        double x1 = a + (i+1) *h ;
        double sumk = 0.0 ;   //interval sum
        double k = 2.0 / Ng ;
        
        for(int j=0 ; j < Ng ; j++){

            sumk += w[j] * F(0.5*(x1-x0)*x[j] + 0.5*(x0 + x1));
            
        }

        sumk *= 0.5* (x1-x0);
        sum += sumk;

        
    }
    return sum;
}

