#include <stdio.h>
#include <stdlib.h>
#include <stdalign.h>
#include <math.h>
#include <string.h>
#include "counter.h"

#define CLS 64
#define MAX_ITER 15000
#define TOL 1e-5

//tamaño do bloque
#define B_SIZE 32

//etiquetas
#define MENOS_INST 0
#define DIV_LAZOS 1
#define DES_LAZOS 2
#define OPER_BLOQUES 3

int opt[] = {0,0,0,0}; //opciones optimizacion

//devolve un número aleatorio entre 0 e 10
double aleatorio()
{
    return ((double)rand() / RAND_MAX) * 10;
}

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        perror("Insuficientes argumentos de entrada\n");
        exit(EXIT_FAILURE);
    }

    const int N = atoi(argv[1]);
    const double suma = 6000 * N;

    srand(N); //fijamos semilla de aleatoriedad

    double** a; //matriz de coeficientes
    double* a_optm;  //matriz de coeficientes optimizada
    double* b = (double*) aligned_alloc(CLS, N * sizeof(double)); //terminoos independientes
    double* x = (double*) aligned_alloc(CLS, N * sizeof(double)); //solución
    double* x_new = (double*) aligned_alloc(CLS, N * sizeof(double)); //nova solución
    double norm2;//norma del vector al cuadrado

    opt[MENOS_INST]   = (argc > 2) ? atoi(argv[2]) : 0;
    opt[DIV_LAZOS]    = (argc > 3) ? atoi(argv[3]) : 0;
    opt[DES_LAZOS]    = (argc > 4) ? atoi(argv[4]) : 0;
    opt[OPER_BLOQUES] = (argc > 5) ? atoi(argv[5]) : 0;

    //reservamos memoria para as filas da matriz
    if(opt[MENOS_INST])
    {
        a_optm = (double*) aligned_alloc(CLS, N * N * sizeof(double));
    }
    else
    {
        a = (double**) aligned_alloc(CLS, N * sizeof(double*));

        for(int i = 0; i < N; i++)
        {
            a[i] = (double*) aligned_alloc(CLS, N * sizeof(double));
        }
    }

    if(opt[DES_LAZOS])
    {
        for(int i = 0; i < N; i+=2)
        {
            b[i] = aleatorio(); //valores aleatorios en b
            x[i] = 0;
            x_new[i] = 0;

            for(int j = 0; j < N; j+=5)
            {
                if(opt[MENOS_INST])
                {
                    a_optm[i * N + j] = aleatorio();
                    a_optm[i * N + j+1] = aleatorio();
                    a_optm[i * N + j+2] = aleatorio();
                    a_optm[i * N + j+3] = aleatorio();
                    a_optm[i * N + j+4] = aleatorio();
                    a_optm[(i+1) * N + j] = aleatorio();
                    a_optm[(i+1) * N + j+1] = aleatorio();
                    a_optm[(i+1) * N + j+2] = aleatorio();
                    a_optm[(i+1) * N + j+3] = aleatorio();
                    a_optm[(i+1) * N + j+4] = aleatorio();
                }
                else
                {
                    a[i][j] = aleatorio();
                    a[i][j+1] = aleatorio();
                    a[i][j+2] = aleatorio();
                    a[i][j+3] = aleatorio();
                    a[i][j+4] = aleatorio();
                    a[i+1][j] = aleatorio();
                    a[i+1][j+1] = aleatorio();
                    a[i+1][j+2] = aleatorio();
                    a[i+1][j+3] = aleatorio();
                    a[i+1][j+4] = aleatorio();
                }
            }
        }
    }
    else
    {
        //inicializamos matrices e vectores
        for(int i = 0; i < N; i++)
        {
            b[i] = aleatorio(); //valores aleatorios en b
            x[i] = 0;
            x_new[i] = 0;

            for(int j = 0; j < N; j++)
            {
                if(opt[MENOS_INST])
                    a_optm[i * N + j] = aleatorio();
                else
                    a[i][j] = aleatorio();
            }
        }
    }

    //aseguramos que a matriz sea diagonal dominante
    for(int i = 0; i < N; i++)
    {
        if(opt[MENOS_INST])
            a_optm[i * N + i] += suma;
        else
            a[i][i] += suma;
    }

    //método de Jacobi
    start_counter();

    for(int iter = 0; iter < MAX_ITER; iter++)
    {
        norm2 = 0;

        //OPERACION POR BLOQUES
        if(opt[OPER_BLOQUES])
        {
            if(opt[MENOS_INST])
            {
                for(int i = 0; i < N; i++)
                {
                    double sigma = 0.0;

                    for(int bj = 0; bj < N; bj += B_SIZE)  //va de bloque en bloque
                    {
                        int fin = (bj + B_SIZE < N) ? bj + B_SIZE : N;

                        for(int j = bj; j < fin; j++) //en cada bloque 
                        {
                            if(i != j)
                                sigma += a_optm[i * N + j] * x[j];
                        }
                    }

                    x_new[i] = (b[i] - sigma) / a_optm[i * N + i];
                    norm2 += (x_new[i] - x[i]) * (x_new[i] - x[i]);
                }
            }
            else
            {
                for(int i = 0; i < N; i++)
                {
                    double sigma = 0.0;

                    for(int bj = 0; bj < N; bj += B_SIZE)
                    {
                        int fin = (bj + B_SIZE < N) ? bj + B_SIZE : N;

                        for(int j = bj; j < fin; j++)
                        {
                            if(i != j)
                                sigma += a[i][j] * x[j];
                        }
                    }

                    x_new[i] = (b[i] - sigma) / a[i][i];
                    norm2 += pow((x_new[i] - x[i]), 2);
                }
            }
        }

        //DIVISION DE LAZOS
        else if(opt[DIV_LAZOS])
        {
            if(opt[MENOS_INST])
            {
                for(int i = 0; i < N; i++)
                {
                    double sigma = 0.0;

                    for(int j = 0; j < i; j++)
                        sigma += a_optm[i * N + j] * x[j];

                    for(int j = i+1; j < N; j++)
                        sigma += a_optm[i * N + j] * x[j];

                    x_new[i] = (b[i] - sigma) / a_optm[i * N + i];
                    norm2 += (x_new[i] - x[i]) * (x_new[i] - x[i]);
                }
            }
            else
            {
                for(int i = 0; i < N; i++)
                {
                    double sigma = 0.0;

                    for(int j = 0; j < i; j++)
                        sigma += a[i][j] * x[j];

                    for(int j = i+1; j < N; j++)
                        sigma += a[i][j] * x[j];

                    x_new[i] = (b[i] - sigma) / a[i][i];
                    norm2 += pow((x_new[i] - x[i]), 2);
                }
            }
        }
        else
        {
            if(opt[MENOS_INST])
            {
                for(int i = 0; i < N; i++)
                {
                    double sigma = 0.0;

                    for(int j = 0; j < N; j++)
                        if(i != j)
                            sigma += a_optm[i * N + j] * x[j];

                    x_new[i] = (b[i] - sigma) / a_optm[i * N + i];
                    norm2 += (x_new[i] - x[i]) * (x_new[i] - x[i]);
                }
            }
            else
            {
                for(int i = 0; i < N; i++)
                {
                    double sigma = 0.0;

                    for(int j = 0; j < N; j++)
                        if(i != j)
                            sigma += a[i][j] * x[j];

                    x_new[i] = (b[i] - sigma) / a[i][i];
                    norm2 += pow((x_new[i] - x[i]), 2);
                }
            }
        }

        if(opt[MENOS_INST])
        {
            double* tmp = x;
            x = x_new;
            x_new = tmp;
        }
        else
        {
            memcpy(x, x_new, sizeof(double) * N); //x = x_new
        }

        //Criterio de parada
        if(!opt[MENOS_INST] && sqrt(norm2) < TOL)
            break;
        else if(norm2 < TOL * TOL)
            break;
    }

    double tiempo = get_counter(); //ciclos de reloj

    printf("Norma: %e\nTiempo: %lf\n", norm2, tiempo); //mostramos la norma

    return 0;
}