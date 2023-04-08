/* definition of functions with solution methods for elliptic pde using Dirichlet b. c. */

#include "elliptic_pde.h"


void Jacobi_Method(double **phi_old, double *x, double *y, void (*BoundaryConditions)(double **,double *, double *, int, int), double tol, int* iterations){
    double **psi_new;
    psi_new = new double* [NX] ;
    psi_new[0] = new double [NX*NY];
    for(int i = 0 ; i<ny ; i++){
        psi_new[i] = psi_new[i-1]+ny;
    }

    //boundary conditions
    Boundary_Conditions(psi_old,x, y, NX, NY);


    for (int i = 0 ; i < nx ; i++){
        for (int j = 0; j < ny; j++){
            psi_new[i][j] = psi_old[i][j] ;
        }
    }

    //iterating jacobi process
    double eps = 1.0; 
    int iterations = 0 ;

    while(eps >= tol){
        for(int i = IBEG ; i <=IEND ; i++){
            for (int j = JBEG ; j <= JEND ; j++){
                psi_new [i][j] = 0.25 *(psi_old[i+1][j] + psi_old[i-1][j] + psi_old[i][j+1] + psi_old[i][j-1] - h*h*S);
            }
        }

        for (int i = 0 ; i < nx ; i++){
            for (int j = 0; j < ny; j++){
                psi_old[i][j] = psi_new[i][j] ;
            }
        }
        eps = 0.0;
        double psi_x;
        double psi_y ;
        for (int i = IBEG ; i <= IEND ; i++){
            for (int j = JBEG; j <= JEND; j++){
                psi_x = psi_old[i+1][j]-2.0*psi_old[i][j]+psi_old[i-1][j];
                psi_y = psi_old[i][j+1]-2.0*psi_old[i][j]+psi_old[i][j-1];
                eps += fabs(psi_x + psi_y - h*h*S);
            }
        
        }
        iterations++;
        cout << "Err = " << eps << endl;
        if(iterations > 5000) break;    
    }

    cout << "Number of iterations: N = " << iterations << endl;


    delete[] psi_new , psi_old;
    delete psi_new[0], psi_old[0]; 
}
