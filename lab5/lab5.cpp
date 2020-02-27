#include <cstdio>
#include <cstdlib>
#include <cmath>

//predefiniowane stale

const double delta = 0.2;
const double d_x = 0.2;
const double d_y = 0.2;

const int n_x = 128;
const int n_y = 128;


const double x_max = delta*n_x;
const double y_max = delta*n_y;
const double TOL = pow(10,-8);


//funkcja stopu

double stop(double V[][n_y+1], int k) 
{
	double s = 0;
	for(int i = 0; i <= n_x-k; i += k)
    {
		for(int j = 0; j <= n_y-k; j += k)
        {
			s += 0.5*pow(k*delta, 2)*(pow((V[i+k][j] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i][j+k])/(2*k*delta), 2) + pow((V[i][j+k] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i+k][j])/(2*k*delta), 2));
		}
	}
	return s;
}


//relaksacja wielosiatkowa
void poisson(FILE * calka, FILE * map)
{
    //zmienne stopu
    double s = 0.;
    double s1 = 0;
    int iter = 0;

    //zaczynamy od najrzadszej siatki
    int k = 16;

    double V[n_x+1][n_y+1];

    //zerowanie tablicy
    for(int i = 0 ; i <= n_x; i++)
    {
        for(int j = 0; j <= n_y; j++)
        {
            V[i][j] = 0;
        }
    }

    //warunki brzegowe Dirichleta
	for(int j = 0; j <= n_y; j++)
    {
		V[0][j] = sin(M_PI*d_y*j/y_max);
		V[n_x][j] = sin(M_PI*d_y*j/y_max);
	}

	for(int i = 0; i <= n_x; i++)
    {
		V[i][n_y] = (-1)*sin(2*M_PI*d_x*i/x_max);
		V[i][0] = sin(2*M_PI*d_x*i/x_max);
	}
   

    while(k >= 1)
    {
        s1 = stop(V, k);
        //dyskretyzacja
        do{
			for(int i = k; i <= n_x-k; i += k)
            {
		        for(int j = k; j <= n_y-k; j += k)
                {
			        V[i][j] = 0.25*(V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);
                }
            }

            //parametry stopu
            s = s1;
			s1 = stop(V, k);

            //przypisanie calki funkcjonalnej
            fprintf(calka, "%d %d %f \n",k, iter, s1 );
            
            iter++;
		}while(fabs((s1-s)/s) > TOL);

        fprintf(calka, "\n \n" );

        //przypisanie map
        for(int i = 0; i <= n_x; i += k)
        {
            for(int j = 0; j <= n_y; j += k)
            {
                fprintf(map, "%f %f %f \n", d_x*i, d_y*j, V[i][j]);
            }
            fprintf(map, "\n");
        }
        fprintf(map, "\n");

        //zageszczanie siatki

	    for(int i = 0; i <= n_x-k; i += k)
        {
		    for(int j = 0; j <= n_y-k; j += k)
            {
                //if potrzebny by nie zniszczyc warunkow dirichleta
                V[i + k/2][j + k/2] = 0.25*(V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
                if(i!=n_x-k)
                {
                    V[i + k][j + k/2] = 0.5*(V[i+k][j] + V[i+k][j+k]);
                }
                if(j!=n_y-k)
                    V[i + k/2][j + k] = 0.5*(V[i][j+k] + V[i+k][j+k]);
                if(j!=0)
                    V[i + k/2][j] = 0.5*(V[i][j] + V[i+k][j]);
                if(i!=0)
                    V[i][j + k/2] = 0.5*(V[i][j] + V[i][j+k]);
		    }
	    }
        //zmniejszanie siatki
	    k = k/2;
	}
}

int main(){
    FILE * calka = fopen("calka.dat", "w");
    FILE * map = fopen("map.dat", "w");
    
    
    poisson( calka, map);
    
    
    fclose(calka);
    fclose(map);

    return 0;
}
