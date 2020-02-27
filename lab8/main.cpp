#include <iostream>
#include </usr/include/gsl/gsl_linalg.h>
#include </usr/include/math.h>

const int nx = 40;
const int ny = 40;
const int N = (nx + 1) * (ny + 1);
const double delta = 1;
const double deltaT = 1;
const double tA = 40;
const double tB = 0;
const double tC = 30;
const double tD = 0;
const int IT_MAX = 2000;

int signum;

const double kB = 0.1;
const double kD = 0.6;

//macierze
/*
double A[N][N];
double B[N][N];
double c[N];
double T[N];
double d[N];
*/

gsl_matrix *A = gsl_matrix_calloc(N, N);
gsl_matrix *B = gsl_matrix_calloc(N, N);
gsl_vector *c = gsl_vector_calloc(N);
gsl_vector *T = gsl_vector_calloc(N);
gsl_vector *d = gsl_vector_calloc(N);

gsl_permutation *p = gsl_permutation_alloc(N);
void wypelnianie()
{
    int l;
    for(int i = 1; i < nx; i++)
    {
        for(int j = 1; j < ny; j++)
        {
            l = i + j * (nx + 1); //chyba tutaj
            gsl_matrix_set(A, l, l-1-nx, deltaT/2/delta/delta);
            gsl_matrix_set(A, l, l-1, deltaT/2/delta/delta);
            gsl_matrix_set(A, l, l+nx+1, deltaT/2/delta/delta);
            gsl_matrix_set(A, l, l+1, deltaT/2/delta/delta);

            gsl_matrix_set(A, l, l, (-1) * 2 * deltaT / delta / delta - 1);

            gsl_matrix_set(B, l, l-nx-1, (-1) * deltaT / 2 / delta / delta);
            gsl_matrix_set(B, l, l-1, (-1) * deltaT / 2 / delta / delta);
            gsl_matrix_set(B, l, l+1, (-1) * deltaT / 2 / delta / delta);
            gsl_matrix_set(B, l, l+1+nx, (-1) * deltaT / 2 / delta / delta);

            gsl_matrix_set(B, l, l, 2 * deltaT / delta / delta - 1);

        }
    }

    //WB Dirichleta

    int i = 0;
    for(int j = 0; j <= ny; j++)   //<=ny
    {
        l = i + j * (nx + 1);
        gsl_matrix_set(A, l, l, 1);
        gsl_matrix_set(B, l, l, 1);
        gsl_vector_set(c, l, 0);
    }
    i = nx;
    for(int j = 0; j <= ny; j++)   //<=ny
    {
        l = i + j * (nx + 1);
        gsl_matrix_set(A, l, l, 1);
        gsl_matrix_set(B, l, l, 1);
        gsl_vector_set(c, l, 0);
    }

    

    int j = ny;
    for(i = 1; i < nx; i++)
    {
        l = i + j * (nx + 1);
        gsl_matrix_set(A, l, l-nx-1,(-1) / kB /delta);
        gsl_matrix_set(A, l, l, 1 + 1 / kB / delta);
        gsl_vector_set(c, l, tB);
        
        for(int z = 0; z < N; z++)
        {
            gsl_matrix_set(B, l, z, 0);
        }
    }

    j = 0;
    for(i = 1; i < nx; i++)
    {   
        l = i + j * (nx + 1);
        gsl_matrix_set(A, l, l,1 + 1 / kD / delta);
        gsl_matrix_set(A, l, l+nx+1, (-1) * 1 / kD / delta);
        gsl_vector_set(c, l,tD);

        for(int z = 0; z < N; z++)
        {
            gsl_matrix_set(B, l, z, 0);
        }
    }

    //warunki poczatkowe temp.

    i = 0;
    for(j = 0; j <=ny; j++)
    {
        l = i + j * (nx + 1);
        gsl_vector_set(T, l, tA);
    } 
    
    i = nx;
    for(j = 0; j <=ny; j++)
    {
        l = i + j * (nx + 1);
        gsl_vector_set(T, l, tC);
    }

    for(i = 1; i < nx; i++)
    {
        for(j = 0; j <= ny; j++)
        {
            l = i + j * (nx + 1);
            gsl_vector_set(T, l, 0);
        }
    } 
    
}

void rozkladLU()
{
    gsl_linalg_LU_decomp(A, p, &signum);
}

int alfa = 1;
int beta = 0;
gsl_vector *BT = gsl_vector_calloc(N);
void CN(double iter)
{
    for(int i = 0; i < iter; i++)
    {
        std::cout<<i<<"\n";
        
            gsl_blas_dgemv(CblasNoTrans, alfa, B, T, beta, d);
            gsl_blas_daxpy(alfa, c, d);
            gsl_linalg_LU_solve(A, p, d, T);
        
    }
}


int main()
{
    FILE * mapa[5];

    mapa[0] = fopen("mapa_0.txt", "w");
    mapa[1] = fopen("mapa_1.txt", "w");
    mapa[2] = fopen("mapa_2.txt", "w");
    mapa[3] = fopen("mapa_3.txt", "w");
    mapa[4] = fopen("mapa_4.txt", "w");

    wypelnianie();
    rozkladLU();

    double z[] ={100,200,500,1000,1500};
    for(int g = 1; g < 6; g++)
    {
        CN(z[g-1]);
        for(int i = 0; i < nx; i++)
        {
            for(int j = 0; j < ny; j++)
            {
               int l = i + j * (nx + 1);
               
               fprintf(mapa[g-1], "%d %d %f \n", i, j, gsl_vector_get(T, l));
            }
            fprintf(mapa[g-1], "\n");
        }
        fprintf(mapa[g-1], "\n");
    }



    for(int i = 0; i < 5; i++)
    {
        fclose(mapa[i]);
    }

    return 0;
}
