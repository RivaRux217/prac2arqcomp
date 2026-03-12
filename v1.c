#include <stdio.h>
#include <stdlib.h>
#include <stdalign.h>
#include <math.h>
#include <string.h>
#include "counter.h"

#define N 1250
#define CLS 64

//Devolve un número aleatorio entre 0 e 10
double aleatorio()
{
    return ((double)rand() / RAND_MAX) * 10;
}

int main()
{
    srand(N); //Fijamos semilla de aleatoriedad

    double** a = (double**) aligned_alloc(CLS, N * sizeof(double*)); //Matriz de coeficientes
    double* b = (double*) aligned_alloc(CLS, N * sizeof(double)); //Vector de termos independientes
    double* x = (double*) aligned_alloc(CLS, N * sizeof(double)); //Vector solución
    double tol = 1e-5; //Tolerancia para a converxencia
    int max_iter = 15000; //Numero de iteracións
    double* x_new = (double*) aligned_alloc(CLS, N * sizeof(double)); //Nova solución
    double norm2;//norma del vector al cuadrado

    //Reservamos memoria para as filas da matriz
    for(int i = 0; i < N; i++)
    {
        a[i] = (double*) aligned_alloc(CLS, N * sizeof(double));
    }

    //Inicializamos matrices e vectores
    for(int i = 0; i < N; i++)
    {
        b[i] = aleatorio(); //Poñemos valores aleatorios en B
        x[i] = 0;
        x_new[i] = 0;

        for(int j = 0; j < N; j++)
        {
            a[i][j] = aleatorio();
        }
    }

    //Aseguramos que a matriz sea diagonal dominante
    for(int i = 0; i < N; i++)
    {
        double suma = 0;

        for(int j = 0; j < N; j++)
        {
            suma += a[i][j];
        }

        a[i][i] += suma;
    }

    //Método de Jacobi
    start_counter();

    for(int iter = 0; iter < max_iter; iter++)
    {
        norm2 = 0;
        
        for(int i = 0; i < N; i++)
        {
            double sigma = 0.0;
            
            for(int j = 0; j < N; j++)
            {
                if(i != j){
                    sigma += a[i][j] * x[j];
                }
            
            }
            
            x_new[i] = (b[i] - sigma) / a[i][i];
            norm2 += pow((x_new[i] - x[i]), 2);
        }

        memcpy(x, x_new, sizeof(x_new) * N); //x = x_new

        if(sqrt(norm2) < tol)
        {
            break;
        }
    }

    double tiempo = get_counter(); //Obtemos ciclos de reloj

    printf("Norma: %e\nTiempo: %lf\n", norm2, tiempo); //mostramos el valor de la norma

    return 0;
}

