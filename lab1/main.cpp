#include <cstdlib>
#include <cstdio>
#include <cmath>

const double lam = -1.0;
const double max_t = 5.0;
const double min_t =0.0;

double sourceOfTense(double omega, double time)
{
	return 10*sin(omega*time);
}

void euler_1(FILE *results, FILE *errors, double step)
{
	double y = 1;
	double dy= 1;
	double y0 = 1;

	double t = min_t;
	while(t<= max_t)
	{
		dy = exp(lam * t);
		
		fprintf(results, "%f %f %f\n", t, y, dy);
		fprintf(errors, "%f %f\n", t, y - dy);

		y = y + step * y * lam;
		y0 = y;
		t = t + step;
	}

}


void RK_2(FILE *results, FILE *errors, double step)
{
	double y = 1;
	double dy= 1;
	double y0 = 1;
	
	double k1;
	double k2;

	double t = min_t;
	while(t<= max_t)
	{
		dy = exp(lam * t);
		
		fprintf(results, "%f %f %f\n", t, y, dy);
		fprintf(errors, "%f %f\n", t, y - dy);
	
		k1 = lam * y;
		k2 = lam * (y + (step * k1));

		y = y + (step / 2.0) * (k1 + k2);
		y0 = y;
		t = t + step;
	}

}

void RK_4(FILE *results, FILE *errors, double step)
{
	double y = 1;
	double dy= 1;
	double y0 = 1;
	
	double k1;
	double k2;
	double k3;
	double k4;

	double t = min_t;
	while(t <= max_t)
	{
		dy = exp(lam * t);
		
		fprintf(results, "%f %f %f\n", t, y, dy);
		fprintf(errors, "%f %f\n", t, y - dy);
	
		k1 = lam * y;
		k2 = lam * (y + (step * k1/2.0));
		k3 = lam * (y + (step * k2/2.0));
		k4 = lam * (y + (step * k3));

		y = y + (step / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4) ;
		y0 = y;
		t = t + step;
	}

}

void RRZ2(FILE * Q_file, FILE * I_file, double k_omega)
{
	double R =100;
	double L = 0.1;
	double C = 0.001;
	double V;
	double I;
	double Q;
	
	double omega_0 = 1/sqrt(L * C);
	double omega = omega_0 * k_omega;
	
	double T_0 = 2 * M_PI / omega_0;
	double t_0 = 0;
	double t_k = 4 * T_0;
	
	double step = 0.001;

	double k1_Q, k2_Q, k3_Q, k4_Q;
	double k1_I, k2_I, k3_I, k4_I;

	double t = t_0;	
	while( t <= t_k)
	{
		V = sourceOfTense(omega, t);
		
	        k1_Q = I;
    		k1_I = V/L - Q/(L*C) - R * I/L;

    		k2_Q = I + k1_I * step/2.0;
    		k2_I = sourceOfTense(omega, t+(step/2.0))/L + (Q + k1_Q * step/2.0)/(L*C) - (R/L)*(I + k1_I * step/2.0);

    		k3_Q = I + k2_I*step/2.0;
    		k3_I = sourceOfTense(omega, t+(step/2.0))/L + (Q + k2_Q * step/2.0)/(L*C) - (R/L) * (I + k2_I * step/2.0);

    		k4_Q = I + step*k3_I;
    		k4_I = sourceOfTense(omega, t+step)/L + (Q + k3_Q * step)/(L*C) - (R/L) * (I + k3_I * step);


		fprintf(Q_file, "%f %f\n", t, Q);
    		fprintf(I_file, "%f %f\n", t, I);

    		Q = Q + (step/6.0)*(k1_Q + 2.0*k2_Q + 2.0*k3_Q + k4_Q);
    		I = I + (step/6.0)*(k1_I + 2.0*k2_I + 2.0*k3_I + k4_I);
		t = t + step;
	}


}

int main()
{
	FILE * results, *errors;
	//Rozwiazanie zadania 1 
	//krok = 1
	
	results = fopen("results1_1.0.dat", "w");
	errors = fopen("errors1_1.0.dat", "w");

	euler_1(results, errors, 1);

	fclose(errors);
	fclose(results);

	//krok = 0.1
	
	results = fopen("results1_0.1.dat", "w");
	errors = fopen("errors1_0.1.dat", "w");

	euler_1(results, errors, 0.1);

	fclose(errors);
	fclose(results);


	//krok = 0.01
	
	results = fopen("results1_0.01.dat", "w");
	errors = fopen("errors1_0.01.dat", "w");

	euler_1(results, errors, 0.01);

	fclose(errors);
	fclose(results);


	//Rozwiazanie zad 2
	//krok = 1
	results = fopen("results2_1.0.dat", "w");
	errors = fopen("errors2_1.0.dat", "w");

	RK_2(results, errors, 1);

	fclose(errors);
	fclose(results);
	

	//krok = 0.1
	results = fopen("results2_0.1.dat", "w");
	errors = fopen("errors2_0.1.dat", "w");

	RK_2(results, errors, 0.1);

	fclose(errors);
	fclose(results);
	
	
	//krok = 0.01
	results = fopen("results2_0.01.dat", "w");
	errors = fopen("errors2_0.01.dat", "w");

	RK_2(results, errors, 0.01);

	fclose(errors);
	fclose(results);

	//Rozwiazanie zad 3
	//krok = 1
	results = fopen("results3_1.0.dat", "w");
	errors = fopen("errors3_1.0.dat", "w");

	RK_4(results, errors, 1);

	fclose(errors);
	fclose(results);
	

	//krok = 0.1
	results = fopen("results3_0.1.dat", "w");
	errors = fopen("errors3_0.1.dat", "w");

	RK_4(results, errors, 0.1);

	fclose(errors);
	fclose(results);
	
	
	//krok = 0.01
	results = fopen("results3_0.01.dat", "w");
	errors = fopen("errors3_0.01.dat", "w");

	RK_4(results, errors, 0.01);

	fclose(errors);
	fclose(results);

	// Zad 4
	
	FILE * Q_file, *I_file;
	
	//omega = 0.5
	Q_file = fopen("Q_file_0.5.dat", "w");
	I_file = fopen("I_file_0.5.dat", "w");

	RRZ2(Q_file, I_file, 0.5);

	fclose(Q_file);
	fclose(I_file);

	//omega = 0.8
	Q_file = fopen("Q_file_0.8.dat", "w");
	I_file = fopen("I_file_0.8.dat", "w");

	RRZ2(Q_file, I_file, 0.8);

	fclose(Q_file);
	fclose(I_file);

	//omega = 1
	Q_file = fopen("Q_file_1.dat", "w");
	I_file = fopen("I_file_1.dat", "w");

	RRZ2(Q_file, I_file, 1);

	fclose(Q_file);
	fclose(I_file);
	
	//omega = 1.2
	Q_file = fopen("Q_file_1.2.dat", "w");
	I_file = fopen("I_file_1.2.dat", "w");

	RRZ2(Q_file, I_file, 1.2);

	fclose(Q_file);
	fclose(I_file);

}
