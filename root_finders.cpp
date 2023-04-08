// include implementations of root finders functions and bracketing function

#include "root_finders.h"


int Bisection(double (*F)(double), double a , double b , double xtol, double ftol , double &z, int nmax, int &k){ //possible to return different error types as different integers

    double dx = fabs(b-a); //interval lenght
    z = a + 0.5*dx; // mid point
    k = 0 ; //iteration counter
    double fa = F(a); // function evaluated in a 
    double fb = F(b) ; //function evaluated in b         
    double fz = F(z);

    if(fa==0) {
        z = a;
        cout << "Extreme a is a root" << endl;
        return 0;
    }

    if (fb==0){
        z = b;
        cout << "Extreme b is a root" << endl;
        return 0;
    }
    while (/*fabs(fz) >=ftol ||*/  dx >= xtol || fz == 0.0){     //choose how to use tollerance: whit xtol and ftol > 0 both are use as condition/ using a large value for xtol only ftol condition is used (and viceversa)     
    
        k++;

        if (k>nmax){
            cout << "Error: too many iterations" << endl;
            return 1;
        }
        dx = fabs(b-a);
        z = a + 0.5*dx;
        fz = F(z) ;              //function evaluated in z
        //cout << "Bisection: k = " << k << "; [a,b] = [" << a << ", " << b << "]; xm = " << z << "; dx = " << dx << "; fm = " << fz << endl;
        
        
        if (fa*fz< 0) {
            b = z;
            fb = fz;
        }
        else {
            a = z;
            fa = fz ;
        }
        
      
       

    }

    return 0;
}


//----------------------------------------------------------------------------------------------------------

int FalsePosition(double (*F)(double) , double a , double b , double xtol ,double ftol,  double &z, int nmax, int& k){     // using a y= m*x+q

    double dx = fabs (b-a);
    k = 0;
    double fa = F(a);
    double fb = F(b);
    z = (a*fb - b*fa) / (fb-fa);
    double fz= F(z);
    double delta = b-a; // measure how much the extreme moves 


    if(fa==0) {
        z = a;
        cout << "Extreme a is a root" << endl;
        return 0;
    }

    if (fb==0){
        z = b;
        cout << "Extreme b is a root" << endl;
        return 0;
    }

    while (fabs(delta) > xtol){  //choose how to use tollerance: whit xtol and ftol > 0 both are use as condition/ using a large value for xtol only ftol condition is used (and viceversa)

        k++;

        if (k>nmax){
            cout << "Error: too many iterations" << endl;
            return 1;
        }
        
        z = (a*fb - b*fa)/(fb-fa);
        fz = F(z);
        //cout << "False Position: k = " << k << "; [a,b] = [" << a << ", " << b << "]; xm = " << z << "; delta = " << delta << "; fm = " << fz << endl;
        
        
        if (fa*fz< 0.0) {
            delta = fabs(b-z);
            b = z;
            fb = fz;
        }
        else {
            delta = fabs(a-z);
            a = z;
            fa = fz ;
        }


    }
    return 0;
}

 //era possibile utilizzare un ciclo for su un massimo di interazioni arbitrarie (es 100) ogni volta si fa un if con le codizioni appropriate ed eventialmente si ritorna 


//-------------------------------------------------------------------------------------------------------------------------------------------

int Secant(double (*F)(double), double a , double c, double tol , double &z, int nmax, int& k){

    //double b = a + fabs(c-a);
    double b = c;
    double fa = F(a);
    double fb = F(b);
    double dx = b-a ;
    k = 0;


    if(fa==0) {
        z = a;
        cout << "Extreme a is a root" << endl;
        return 0;
    }

    if (fb==0){
        
        z = b;
        cout << "Extreme b is a root" << endl;
        return 0;
    }

    while(fabs(dx) >= tol){      //b is the therical zero

        k++;
        if (k>nmax){
            cout << "Error: too many iterations" << endl;
            return 1;
        }
        dx = fb*(b-a)/(fb-fa); //theoric increment

        //cout << "Secant: k = " << k << "; a = " << a << "; b = " << b << "; dx = " << dx << "; f(b) = " << fb << endl;

        a = b;
        fa = fb;
        b = b - dx;       //redefine b and F(b)
        fb = F(b);

    }

    z = b;

    return 0;
}

//--------------------------------------------------------------------------------------------------

int Newton(double (*F)(double) , double (*F_prime)(double), double a, double c,double tol, double &z, int nmax, int& k){

    double b = a + fabs(c-a);
    //double b = c;
    double fa = F(a);
    double fb = F(b);
    double fb_prime = F_prime(b);
    double dx = b-a ;
    k = 0;


    if(fa==0) {
        z = a;
        cout << "Extreme a is a root" << endl;
        return 0;
    }

    if (fb==0){
        
        z = b;
        cout << "Extreme b is a root" << endl;
        return 0;
    }

    while(fabs(dx) >= tol){

        k++;
        if (k>nmax){
            cout << "Error: too many iterations" << endl;
            return 1;
        }

        dx = fb/fb_prime ;
        //cout << "Newton: k = " << k << "; a = " << a << "; b = " << b << "; dx = " << dx << "; f'(b) =" <<fb_prime << "; f(b) = " << fb << endl;

        a = b;
        b -=dx ;
        fb = F(b);
        fb_prime = F_prime(b);

    }

    z = b; 

    return 0;
}

//--------------------------------------------------------------------------------------------------

void Bracketing( double (*F) (double) , double a, double b , int N , double* xL , double*  xR , int& nzero){

    double dx =fabs(b-a)/(double)N;
    double xi = a  , xf = a +dx ;
    double Fi , Fe ;
    nzero = 0;
    Fi = Fe =  F(xi);
    
    for(int i = 0; i<N; i++){
        xi = a +i*dx;
        xf = a + (i+1)*dx;
        Fi = Fe ;
        Fe = F(xf);
     
        if( Fi*Fe<0){    // posso usare una chiamata a F all'inizio e solo una durante il for , provo
            xL[nzero]=xi;
            xR[nzero]=xf;
            nzero++;
        }

    }

}
