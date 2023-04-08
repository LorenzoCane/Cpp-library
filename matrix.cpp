/* 24/11/2022
 write a simple program that creates a NxN matrix A using either standard native or dynamical method (use, e.g. N = 4).
*/

// to be compiled with: g++ matrix.cpp matrix_tools.cpp -o matrix
#include "matrix_tools.h"

void Func (double **arr);
void MatrixVectMult(double **a , double **b, double nrows, double **c) ;

int main(){
    cout << fixed << setprecision(4);

    int nrow , ncol ;
    nrow = ncol = 4 ; 

    double **A ;   // double pointer
    A = new double* [nrow] ; // a contains nrows pointers
    A[0] = new double [nrow*ncol] ; // the first pointer contains all the matrix elements
    for (int i = 1 ; i < nrow ; i++){
        A[i] = A[i-1] + ncol ;  // pointer algebra : every pointer starts ncol element after : points to the following row 
    }

    Func(A);

    //print the matrix
    PrintMatrix(A , nrow , ncol);

    double **B;
    int brow = 4;
    int bcol = 1 ;
    B = new double *[brow] ;
    B[0] = new double [bcol*brow];
    for (int i = 1 ; i < brow ; i++){
        B[i] = B[i-1] + bcol ;  // pointer algebra : every pointer starts ncol element after : points to the following row 
    }

    B[0][0] = 1.0 , B[1][0] = 0.0 ,  B[2][0] = 3.0 ; B[3][0] = 2.0 ;

    double **c;
    int crow = nrow;
    int ccol = bcol ;
    c = new double *[nrow] ;
    c[0] = new double [bcol*nrow];
    for (int i = 1 ; i < nrow ; i++){
        c[i] = c[i-1] + bcol ;  // pointer algebra : every pointer starts ncol element after : points to the following row 
    }
    
    MatrixVectMult(A, B , nrow, c);

    cout << endl << endl ;

    for(int i=0 ; i < nrow ; i++){
        for(int j = 0 ; j < bcol ; j++){
            cout << setw(10) << right << c[i][j] << "  ";
        }
        cout << endl;
    }

    delete A[0];
    delete A;
    delete B[0];
    delete B;

    return 0;
}

//******************************************************************************

void Func(double **arr){
    arr[0][0] = 1.0,  arr[0][1] = 3.0 , arr[0][2] = 2.0 , arr [0][3] = -4.0;
    arr[1][0] = 7.0 ,arr[1][1] = 2.0 , arr[1][2] = 4.0 , arr[1][3] =1.0 ;
    arr[2][0] = 0.0 , arr[2][1]= -1.0 , arr[2][2] = 2.0, arr[2][3] = 2.0 ;
    arr[3][0] = 6.0 , arr[3][1] = 3.0 , arr[3][2] = 0.0 , arr[3][3] =  1.0 ;
}

//-------------------------------------------------------------------------------------
