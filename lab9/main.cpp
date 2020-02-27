#include <iostream>
#include <math.h>

const int nx = 150;
const int nt = 1000;
const double delta = 0.1;
const double deltaT = 0.05;
const double xA = 7.5;
const double sigma = 0.5;
const double xF = 2.5;

double alfa = 0;
double beta = 0;

//dyskretyzacja

double x = 0;

//tablice
//sprawdzić rozmiar
double u0[nx +1] = {0};
double u[nx+1]={0};
double v[nx+1]={0};
double vp[nx+1]={0};
double a[nx+1]={0};

double aF(double x, double t)
{
    double toRet = 0.0;
	if(x == xF){
		double tmax = deltaT*nt;
		toRet = cos(50.0*t/tmax);
	}
	return toRet;
}

void warunkiBrzegowe()
{
    u[0] = 0;
    u[nx] = 0;

    v[0] = 0;
    v[nx] = 0;
}

void warunekPoczatkowy()
{
    for(int i = 1; i < nx ; i++) //pomijamy miejsce złączenia na poczatku i koncu
    {
        x = delta * i;
        u[i] = exp((-1) * (x - xA) * (x - xA)/ 2 / sigma/ sigma);
        v[i] = 0;
    }
}

void algorytm(FILE *wykres, FILE *mapa)
{
    //zachowanie poprzedniego wyniku
    for(int i = 0; i < nx; i++)
    {
        u0[i] = u[i];
    }

    //inicjalizacja wzor 11
    if(alfa < 0.8)
    {
        for(int i = 1; i < nx; i++)
        {
            a[i] = (u[i+1] - 2 * u[i] + u[i-1])/delta /delta - beta * (u[i] - u0[i])/ deltaT + alfa * aF(delta * i, 0);
        }
    }
    for(int k = 1; k <= nt; k++)
    {
        double t = 0.0;
        double E = 0.0;

        for(int i = 0; i < nx; i++){
			vp[i] = v[i] + (deltaT/2.0) * a[i];
		}

		for(int i = 0; i < nx; i++){
			u0[i] = u[i];
		}

		for(int i = 0; i < nx; i++){
			u[i] = u[i] + deltaT*vp[i];
		}

		for(int i = 1; i < nx; i++){
			a[i] = ((u[i+1] - 2.0*u[i] + u[i-1]) / pow(delta, 2.0)) - (beta*(u[i] - u0[i])/deltaT) + alfa*aF(delta*i, deltaT*k);
		}

		for(int i = 0; i < nx; i++){
			v[i] = vp[i] + (deltaT/2.0)*a[i];
		}

		for(int i = 1; i < nx; i++)
		{
			t += pow(v[i], 2.0) + pow((u[i+1] - u[i - 1]) / (delta * 2.0), 2.0);
		}
		E = 0.25*delta*(pow((u[1]-u[0])/delta, 2.0) + pow((u[nx]-u[nx-1])/delta, 2.0)) + 0.5*t*delta;

        fprintf(wykres, "%f %f\n", deltaT * k, E);
		for(int i = 0; i < nx; i++)
		{
			fprintf(mapa, "%f %f %f\n", deltaT*k, delta*i, u[i]);
		}
		fprintf(mapa, "\n");

    }
    
}



int main()
{
    FILE * wykres[4];
    wykres[0] = fopen("wykres_0.txt", "w");
    wykres[1] = fopen("wykres_0.1.txt", "w");
    wykres[2] = fopen("wykres_1.txt", "w");
    wykres[3] = fopen("wykres_1alfa.txt", "w");

    FILE * mapa[4];

    mapa[0] = fopen("mapa_0.txt", "w");
    mapa[1] = fopen("mapa_0.1.txt", "w");
    mapa[2] = fopen("mapa_1.txt", "w");
    mapa[3] = fopen("mapa_1alfa.txt", "w");

    warunkiBrzegowe();
    warunekPoczatkowy();

    algorytm(wykres[0], mapa[0]);

    beta = 0.1;
    //warunekPoczatkowy();
    algorytm(wykres[1], mapa[1]);

    beta = 1;
    //warunekPoczatkowy();
    algorytm(wykres[2], mapa[2]);

    alfa = 1;
    algorytm(wykres[3], mapa[3]);



    for(int i = 0; i < 3; i++)
    {
        fclose(mapa[i]);
        fclose(wykres[i]);
    }

    return 0;
}
