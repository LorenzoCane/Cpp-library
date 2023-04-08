// header file of matrix and linear system of equations useful function

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

void PrintMatrix(double **a, int nrow, int ncol);

void MatrixVectMult(double **a , double **b , double nrow, double **c);

void Rows_swap(double **a,double *b, int j , int k, int ncol);

void Gauss_Elim(double **A , double *B , int nrow);

void Gauss_Solver(double **A, double *B, double *x, int nrow);

void Tridiag_solver(double *a , double *b , double *c , double *r , double *x , int nrow, bool show);


