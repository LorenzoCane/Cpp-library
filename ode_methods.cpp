/* Implementations of ode resolver methods

*/

void EulerStep(double t , double *Y, void (*RHSFunc)(double,double *,double *) , double dt , int neq){
    
    double rhs[neq];      // create a neq dimensions array
    RHSFunc(t, Y , rhs);  //compute rhsfunc and store the value in rhs

    for (int i = 0 ; i< neq ; i++ ){
        Y[i]  += dt*rhs[i];
    }

}

//-------------------------------------------------------------------------------------------------------------

void RK2Step ( double t , double *Y , void (*RHSFunc)(double , double * ,double *) , double dt , int neq){

    double Y1[neq] , k1[neq] , k2[neq];

    RHSFunc(t, Y , k1);

    for(int i = 0 ; i < neq ; i++)  Y1[i] = Y[i]+0.5*dt*k1[i];

    RHSFunc(t+0.5*dt , Y1 , k2);

    for(int i = 0 ; i < neq ; i++) Y[i] += dt*k2[i];


}

//---------------------------------------------------------------------------------------------------------------

void RK4Step  (double t , double *Y , void (*RHSFunc)(double , double * , double *) , double dt , int neq){

    double Y1[neq] , Y2[neq] , Y3[neq] , k1[neq] , k2[neq], k3[neq] , k4[neq] ;

    RHSFunc(t , Y , k1);

    for (int i = 0 ; i < neq ; i++) Y1[i] = Y[i] + 0.5*dt*k1[i];
    RHSFunc(t+0.5*dt , Y1 , k2);

    for (int i = 0 ; i < neq ; i++) Y2[i] = Y[i] + 0.5*dt*k2[i];
    RHSFunc(t+0.5*dt , Y2 , k3);

    for (int i = 0 ; i < neq ; i++) Y3[i] = Y[i] + dt*k3[i];
    RHSFunc(t+dt , Y3 , k4);

    for(int i = 0 ; i < neq ; i++) Y[i] += dt/6.0*(k1[i] +2.0*k2[i]+ 2.0*k3[i] +k4[i]);

}

//-----------------------------------------------------------------------------------------------------------------

void PositionVerlet (double t , double *x , double *v , void(*acc)(double , double * , double *) , double dt , int npart){

    double x_half[npart] ;
    for (int i = 0 ; i < npart ; i++ ) x_half[i] = x[i] + 0.5*dt*v[i] ;  // compute x at n+1/2 time step

    double accel[npart];   // accelleration array
    acc(t , x_half , accel);

    for (int j = 0 ; j < npart ; j++ ) v[j] += dt *accel[j] ;   // v at time n+1
    for (int k = 0 ; k < npart ; k++ ) x[k] = x_half[k] + 0.5 * dt * v[k] ; // space at time n+1

}

