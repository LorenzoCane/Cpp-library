// definitions of useful functions for matrices and linear system of eqations

#include "matrix_tools.h"

void PrintMatrix(double **a, int nrow , int ncol){
    // print matrix of nrow rows and ncol columns, it works on double pointers
     for(int i=0 ; i < nrow ; i++){
        for(int j = 0 ; j < ncol ; j++){
            cout << setw(10) << right << a[i][j] << "  ";
        }
        cout << endl;
    }
}

//------------------------------------------------------------------------------------------------------------------
void MatrixVectMult(double **a , double **b , double nrow, double **c){
    //multipliies a *b and store the result in c 
    for (int k = 0 ; k < nrow ; k++){
        c[0][k] = 0.0;
        for ( int j = 0 ; j < nrow ; j++){                 // da rifare
            c[0][k] += a[k][j] * b[j][0];
        }
       
    }
}

//--------------------------------------------------------------------------------------------------------------------

void Rows_swap(double **a, double *b , int j , int k, int ncol){
    //swap the rows j and k of a linear system of eq, given by the **a double matrix with ncol columns and the vector b
    double tmp[ncol]; 
    for (int i = 0 ; i < ncol ; i++){
        tmp[i] = a[j][i];
        a[j][i] = a[k][i];
        a[k][i] = tmp[i];
    }
     
    double temp = b[k];
    b[k] = b [j];
    b[j] = temp;
}

//--------------------------------------------------------------------------------------------------------------------------------
void Gauss_Elim(double **A , double *B , int nrow){
    // use gauss elimination process (with partial pivoting) to symplify a linear system of equation composed by square matrix a and vector b
    
    double g;
    int ncol     = nrow;
    for (int k = 0 ; k < nrow - 1 ; k++){  // we have to modify n-1 rows, the first one is untouched so we have n-1 G matrices
        int tmax = k;
        double Amax = fabs(A[k][k]);
        for (int t = k+1 ; t < nrow ; t++){
            double temp = fabs(A[t][k]) ; 
            if (temp > Amax){
                tmax = t;
                Amax = temp; 
            } 
        }
        if (k != tmax) Rows_swap(A,B, k , tmax, ncol);
       
        for (int i= k+1 ; i < nrow ; i++){    // loop over the rows
            g = A[i][k]/A[k][k] ;
            for( int j = k+1 ; j< nrow ; j++ ){
                A[i][j] -= g*A[k][j];
            } 
            A[i][k] = 0.0 ;  //first element put manually to zero
            B[i] -= g*B[k] ;          
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------------
void Gauss_Solver(double **A, double *B, double *x, int nrow){
    /*Solves linear equations system that have been previously simplified with Gauss_Elim
     the solution is stored in the vector x , A is the square matrix , B is the known factors vector 
    */
    for (int i = nrow -1 ; i >=0 ; i--){
        double tmp;
        tmp = B[i];
        for ( int j = nrow-1 ; j > i ; j--) tmp -= x[j] *  A[i][j] ;
        x [i] = tmp/A[i][i];
    }
}
 //---------------------------------------------------------------------------------------------------------------------------------
void Tridiag_solver(double *a , double *b , double *c , double *r , double *x, int nrow, bool show){
   //tridiagonal solver a is the under diagonal , b principal diagonal, c upper diagonal , r knowed term , x solution vector
   // show = true print the result of the solver

    double h[nrow] , p[nrow] ;
    
    h[0] = c[0]/b[0];
    for(int i = 1 ; i < nrow ; i++){
        h[i] = c[i] / (b[i] - a[i] * h[i-1]);    // h has an extra element that is not used
    }

    p[0] = r[0]/b[0];
    for(int i = 1 ; i < nrow ; i++){
        p[i] = (r[i] - a[i] * p[i-1]) / (b[i] - a[i] * h[i-1]);
    }

    //solving

    x[nrow-1] = p[nrow -1];
    for(int i = nrow -1 ; i >= 0  ; i--){  // the process goes backward 
        x[i] = p[i] - h[i] * x[i+1] ;
    }

    //print
    if (show){
        for(int j = 0 ; j < nrow ; j++){
            cout << "Solution of the tridiagonal problem:" << setw(10) << right << x[j] << "  " << endl;
        }
        cout << endl;
    }    
}