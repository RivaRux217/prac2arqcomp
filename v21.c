#include <stdio.h>
#include <stdlib.h>
#include <stdalign.h>
#include <math.h>
#include <string.h>
#include "counter.h"

#define CLS 64
#define MAX_ITER 15000
#define TOL 1e-5

//Tamaño do bloque
#define B_SIZE 32


//Devolve un número aleatorio entre 0 e 10
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

    srand(N); //Fijamos semilla de aleatoriedad

    double** a; //Matriz de coeficientes
    double* a_optm;  //Matriz de coeficientes optimizada
    double* b = (double*) aligned_alloc(CLS, N * sizeof(double)); //Vector de termos independientes
    double* x = (double*) aligned_alloc(CLS, N * sizeof(double)); //Vector solución
    double* x_new = (double*) aligned_alloc(CLS, N * sizeof(double)); //Nova solución
    double norm2;//norma del vector al cuadrado


    //Reservamos memoria para as filas da matriz
    a_optm = (double*) aligned_alloc(CLS, N * N * sizeof(double));


    //Inicialización
    for(int i = 0; i < N; i++)
    {
        b[i] = aleatorio();
        x[i] = 0.0;  
        x_new[i] = 0.0;

        for(int j = 0; j < N; j++)
        {
            a_optm[i * N + j] = aleatorio();
        }
    }

    // Aseguramos que la matriz sea diagonal dominante (Esto ya lo tenías bien)
    for(int i = 0; i < N; i++)
    {
        a_optm[i * N + i] += suma;
    }

    //Aseguramos que a matriz sea diagonal dominante
    for(int i = 0; i < N; i++)
    {
        a_optm[i * N + i] += suma;
    }

    //Método de Jacobi
    start_counter();

    for(int iter = 0; iter < MAX_ITER; iter++)
    {
        norm2 = 0;

        
    //DESENRROLLAMIENTO DE LAZOS

        for(int i = 0; i < N; i++)
        {
            double sigma = 0.0;
            int j;

            for(j = 0; j < N-4; j+=4){
                if(i != j)  sigma += a_optm[i * N + j] * x[j];
                if(i != (j+1))  sigma += a_optm[i * N + j + 1] * x[j + 1];
                if(i != (j+2))  sigma += a_optm[i * N + j + 2] * x[j + 2];
                if(i != (j+3))  sigma += a_optm[i * N + j + 3] * x[j + 3];               
            }

            for(; j < N; j++){
                if(i != j)   sigma += a_optm[i * N + j] * x[j];
            }

            x_new[i] = (b[i] - sigma) / a_optm[i * N + i];
            norm2 += (x_new[i] - x[i]) * (x_new[i] - x[i]);
        }
            
        double* tmp = x;
        x = x_new;
        x_new = tmp;
        

        //Criterio de parada
        if(sqrt(norm2) < TOL)
            break;

    }

    double tiempo = get_counter(); //Obtemos ciclos de reloj

    printf("Norma: %e\nTiempo: %lf\n", norm2, tiempo); //mostramos el valor de la norma

    return 0;
}
