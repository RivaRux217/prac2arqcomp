#include <stdio.h>
#include <stdlib.h>
#include <stdalign.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>
#include "counter.h"

#define CLS 64 //Aquí alineamos a 32, necesario para AVX
#define DBL 4 //Número de doubles que podemos tratar cunha sola instrucción AVX
#define MAX_ITER 15000
#define TOL 1e-5

//Devolve un número aleatorio entre 0 e 10
double aleatorio()
{
    return ((double)rand() / RAND_MAX) * 10;
}

int main(int argc, char** argv)
{

    /**
     * FLAGS
     * N -> Número de filas e columnas da matriz
     */
    const int N = (argc == 2) ? atoi(argv[1]) : 1250;


    printf("Iniciando con N = %d\n\n", N);

    srand(N); //Fijamos semilla de aleatoriedad

    int nFila = ((N + 3) / 4) * 4; //Forzamos que las filas de la matriz estén alineadas con 32
    double* a = (double*) aligned_alloc(CLS, nFila * N * sizeof(double)); //Matriz de coeficientes
    double* b = (double*) aligned_alloc(CLS, N * sizeof(double)); //Vector de termos independientes
    double* x = (double*) aligned_alloc(CLS, N * sizeof(double)); //Vector solución
    double* x_new = (double*) aligned_alloc(CLS, N * sizeof(double)); //Nova solución
    double norm2;//norma del vector al cuadrado
    int limit = (N / DBL) * DBL; //Variable necesaria para evitar acceder a posiciones inválidas con N no múltiplo de 4

    //Inicializamos matrices e vectores
    for(int i = 0; i < N; i++)
    {
        b[i] = aleatorio(); //Poñemos valores aleatorios en B
        x[i] = 0;
        x_new[i] = 0;
        double suma = 0;
        
        for(int j = 0; j < N; j++)
        {
            a[i * nFila + j] = aleatorio();
        }

        suma += 6000 * N;

        a[i * nFila + i] += suma;  //Aseguramos que a matriz sea diagonal dominante
    }

    //Método de Jacobi
    start_counter();
    
    for(int iter = 0; iter < MAX_ITER; iter++)
    {
        norm2 = 0;

        for(int i = 0; i < N; i++)
        {
            double sigma = 0.0;
            __m256d vTemp = _mm256_setzero_pd(); //Seteamos registro de guardado temporal de sumas a 0

            /**
             * DESENROLLO DE LAZOS
             * Con instrucciones escalares y desenrollo de lazos realizamos 4 iteraciones en una
             * Con instrucciones vectoriales sin desenrollo de lazos operamos 4 valores en una instrucción
             * Con AVX + desenrollo de lazos realizamos 4 sumas vectoriales de 4 doubles cada una
             * Entonces, j avanza de 16 en 16 para realizar las 4 operaciones en 1
             */
            for(int j = 0; j < limit; j += 16)
            {
                __m256d vA, vX;
                int coef = 0; //Valor sumado a j para el desenrollo de lazos

                /* Iteración 1 */
                vA = _mm256_load_pd(a + (i * nFila + j)); //Guardamos 4 posiciones de a en registro vectorial
                vX = _mm256_load_pd(x + j); //Guardamos 4 posiciones de x en registro vectorial

                //Realizamos la multiplicación vectorial y la sumamos a los resultados guardados en vTemp
                vTemp = _mm256_add_pd(vTemp, _mm256_mul_pd(vA, vX)); 

                /* Iteración 2 */
                coef += 4;
                vA = _mm256_load_pd(a + (i * nFila + j + coef));
                vX = _mm256_load_pd(x + j + coef);

                vTemp = _mm256_add_pd(vTemp, _mm256_mul_pd(vA, vX)); 

                /* Iteración 3 */
                coef += 4;
                vA = _mm256_load_pd(a + (i * nFila + j + coef));
                vX = _mm256_load_pd(x + j + coef);

                vTemp = _mm256_add_pd(vTemp, _mm256_mul_pd(vA, vX)); 

                /* Iteración 4 */
                coef += 4;
                vA = _mm256_load_pd(a + (i * nFila + j + coef));
                vX = _mm256_load_pd(x + j + coef);

                vTemp = _mm256_add_pd(vTemp, _mm256_mul_pd(vA, vX)); 
            }


            //Bucle de limpieza (se suman los valores sobrantes del vector si N no es múltiplo de 4)
            for(int j = limit; j < N; j++)
            {
                sigma += a[i * nFila + j] * x[j];
            }

            sigma -= a[i * nFila + i] * x[i]; //Restamos el elemento diagonal

            //Este bloque se realiza fuera del bucle para evitar afectar a la eficiencia de AVX
            alignas(CLS) double tmp[4]; //Array estático para guardar la suma en sigma
            _mm256_store_pd(tmp, vTemp); //Guardamos valores de vTemp en tmp
            sigma += (tmp[0] + tmp[1] + tmp[2] + tmp[3]); //Sumamos las 4 posiciones de vTemp

            x_new[i] = (b[i] - sigma) / a[i * nFila + i];

            norm2 += (x_new[i] - x[i]) * (x_new[i] - x[i]);
        }

        //x = x_new
        double* tmp = x;
        x = x_new;
        x_new = tmp;
        
        if(norm2 < TOL * TOL)
        {
            break;
        }
    }


    double tiempo = get_counter(); //Obtemos ciclos de reloj

    printf("Norma: %e\nTiempo: %lf\n", norm2, tiempo); //mostramos el valor de la norma

    return 0;
}

