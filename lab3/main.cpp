#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>

double x0 = 0.01;
double v0 = 0.0;
double dt0 = 1.0;
double S = 0.75;
double p = 2.0;
double tmax = 40;
double a = 5;

//by nie zwracac tablicy
struct dane
{
    double x;
    double v;
};


//funkcje do zad 1
double f(double vn)
{
    return vn;
}

double g(double xn, double vn)
{
    return (a*(1-xn*xn)*vn - xn);
}

dane m_RK2(double xn, double vn, double dt)
{
    dane tmp;
    
    double k1x;
    double k1v;
    
    double k2x;
    double k2v;

    k1x = f(vn);
    k1v = g(xn, vn);

    k2x = vn + dt * k1v;
    k2v = g(xn + dt*k1x, vn + dt*k1v);

    tmp.x = xn + (dt/2)*(k1x + k2x);
    tmp.v = vn + (dt/2)*(k1v + k2v);
    return tmp;
}

void kontrola_kroku_RK2(FILE * file, double TOL)
{
    double t = 0.0;
    double  dt = dt0;
    double xn = x0;
    double vn = v0;
    double tmax = 40.0;

    double Ex = 0.0;
    double Ev = 0.0;

    dane tmp1;
    dane tmp2;

   
    fprintf(file, "%f %f %f %f \n", t, dt, xn, vn);
    do{
        tmp2 = m_RK2(xn,vn,dt);
        tmp2 = m_RK2(tmp2.x,tmp2.v,dt);


        tmp1 = m_RK2(xn,vn,2*dt);


        Ex = (tmp2.x - tmp1.x)/(pow(2,p) - 1);
        Ev = (tmp2.v - tmp1.v)/(pow(2,p) - 1);

        if(std::max(fabs(Ex),fabs(Ev)) < TOL ){
            t =t+2*dt;
            xn = tmp2.x;
            vn = tmp2.v;
            fprintf(file, "%f %f %f %f \n", t, dt, xn, vn);

        }
        dt = (pow(S*TOL/(std::max(fabs(Ex),fabs(Ev))),(1/(p+1)))*dt);

    }while(t < (tmax - dt));

}

//operacje na macierzach - zad 2
double a11(){
    return 1.;
}

double a12(double dt){
    return (-1)*(dt/2);
}

double a21(double dt, double xn1k, double vn1k){
    return (-1)*(dt/2)*((-2)*a*xn1k*vn1k-1);
}

double a22(double dt, double xn1k){
    return (1 - (dt/2)*a*(1-xn1k*xn1k));
}

//funkcje do zad 2
double F(double xn1, double xn, double vn1, double vn, double dt)
{
    return xn1 - xn - (dt/2)*(f(vn)+f(vn1));
}

double G(double xn1, double xn, double vn1, double vn, double dt)
{
    return vn1 - vn - (dt/2)*(g(xn,vn)+g(xn1,vn1));
}

dane m_trapezow(double xn, double vn, double dt)
{

    dane tmp;
    double delta = pow(10,-10); //krok

    double xn1k = xn; //Przypisanie poczatkowych
    double vn1k = vn;
    double xn1 = xn;
    double vn1 = vn;
    double dx, dv;

    do
    {
        
        dx = ((-1)*F(xn1,xn, vn1, vn, dt)*a22(dt,xn1k)-(-1)*G(xn1,xn,vn1,vn,dt)*a12(dt))/
                                            (a11()*a22(dt,xn1k) - a12(dt)*a21(dt,xn1k,vn1k));

        dv = (a11()*(-1)*G(xn1,xn,vn1,vn,dt) - a21(dt,xn1k,vn1k)*(-1)*F(xn1,xn,vn1,vn,dt))/
                                            (a11()*a22(dt,xn1k) - a12(dt)*a21(dt,xn1k,vn1k));

        xn1 = xn1k + dx;
        vn1 = vn1k + dv;

    }while(fabs(dx) < delta && fabs(dv) < delta); //warunek stopu

    
    tmp.x = xn + (dt/2)*(f(vn) + f(vn1)); //przypisanie ostatecznych wartości 
    tmp.v = vn + (dt/2)*(g(xn,vn) + g(xn1,vn1));

    return tmp;
}


void kontrola_kroku_trapez(FILE * file, double TOL)
{
    double t = 0.0;
    double  dt = dt0;
    double xn = x0;
    double vn = v0;

    double Ex = 0.0;
    double Ev = 0.0;

    dane tmp1;
    dane tmp2;

   
    fprintf(file, "%f %f %f %f \n", t, dt, xn, vn);
    do{
        tmp2 = m_trapezow(xn,vn,dt);
        tmp2 = m_trapezow(tmp2.x,tmp2.v,dt);


        tmp1 = m_trapezow(xn,vn,2*dt);


        Ex = (tmp2.x - tmp1.x)/(pow(2,p) - 1);
        Ev = (tmp2.v - tmp1.v)/(pow(2,p) - 1);

        if(std::max(fabs(Ex),fabs(Ev)) < TOL ){
            t =t+2*dt;
            xn = tmp2.x;
            vn = tmp2.v;
            fprintf(file, "%f %f %f %f \n", t, dt, xn, vn);

        }
        dt = (pow(S*TOL/(std::max(fabs(Ex),fabs(Ev))),(1/(p+1)))*dt);

    }while(t < (tmax - dt));

}

int main()
{

    FILE * file;
//  Metoda RK2
    double TOL = pow(10,-2);
    file = fopen("RK2_pow-2.dat", "w");
    kontrola_kroku_RK2(file, TOL);
     fprintf(file, "\n"); 
   
     fclose(file);


    file = fopen("RK2_pow-5.dat", "w");
    TOL = pow(10,-5);
    kontrola_kroku_RK2(file, TOL);

    fclose(file);


// Metoda trapezow

    file = fopen("trapezy_-2.dat", "w");
    TOL = pow(10,-2);
    kontrola_kroku_trapez(file, TOL);

    fprintf(file, "\n");
    
     fclose(file);


    file = fopen("trapezy_-5.dat", "w");
    TOL = pow(10,-5);
    kontrola_kroku_trapez(file, TOL);

    fclose(file);
}
