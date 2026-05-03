#include <stdio.h>
#include <stdlib.h>
#include <stdalign.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "counter.h"

#define CLS 64
#define MAX_ITER 15000
#define TOL 1e-5

//Devolve un número aleatorio entre 0 e 10
double aleatorio()
{
    return ((double)rand() / RAND_MAX) * 10;
}

int main(int argc, char** argv)
{
    char* sch[] = {
        "Static",
        "Dynamic",
        "Guided",
        "Auto"
    };
    
    char* crt[] = {
        "Reduction",
        "Critical",
        "Atomic"
    };

    /**
     * FLAGS
     * N -> Número de filas e columnas da matriz
     * TH -> Número de fíos
     * SC -> Tipo de schedule. 0, default: static. 1: dynamic. 2: guided. 3: auto.
     * CRIT -> Decide uso de critical, reduction ou atomic en rexións críticas. 
     *         <1 ou non especificado: usa reduction
     *         1: usa critical
     *         >1: usa atomic
     */
    const int N    = (argc >= 2) ? atoi(argv[1]) : 1250;
    int TH, SC, CRIT;

    if(argc <= 2)
    {
        TH = 32;
        SC = 2;
        CRIT = 0;
    }
    else
    {
        TH   = (argc >= 3) ? atoi(argv[2]) : 1;
        SC   = (argc >= 4) ? atoi(argv[3]) : 0;
        CRIT = (argc >= 5) ? atoi(argv[4]) : 0; //Por defecto usamos reduction
    }

    int critIx = (CRIT == 0) ? 0 : ((CRIT == 1) ? 1 : 2);

    printf("Iniciando con:\n\tN: %d\n\tHilos: %d\n\tModo de schedule: %s\n\tRexións críticas con %s\n\n", 
        N, TH, sch[SC], crt[critIx]);

    srand(N); //Fijamos semilla de aleatoriedad

    double* a = (double*) aligned_alloc(CLS, N * N * sizeof(double)); //Matriz de coeficientes
    double* b = (double*) aligned_alloc(CLS, N * sizeof(double)); //Vector de termos independientes
    double* x = (double*) aligned_alloc(CLS, N * sizeof(double)); //Vector solución
    double* x_new = (double*) aligned_alloc(CLS, N * sizeof(double)); //Nova solución
    double norm2;//norma del vector al cuadrado

    //Inicializamos matrices e vectores
    for(int i = 0; i < N; i++)
    {
        b[i] = aleatorio(); //Poñemos valores aleatorios en B
        x[i] = 0;
        x_new[i] = 0;
        double suma = 0;
        
        for(int j = 0; j < N; j++)
        {
            a[i * N + j] = aleatorio();
        }

        suma += 6000 * N;

        a[i * N + i] += suma;  //Aseguramos que a matriz sea diagonal dominante

    }

    //Elegimos modo de schedule (seteamos chunk_size a 0 para que el compilador decida)
    switch (SC)
    {
        //Modo static
        case 0:
        default:
            omp_set_schedule(omp_sched_static, 0);
        break;
    
        case 1:
            omp_set_schedule(omp_sched_dynamic, 0);
        break;

        case 2:
            omp_set_schedule(omp_sched_guided, 0);
        break;

        case 3:
            omp_set_schedule(omp_sched_auto, 0);
        break;
    }


    //Método de Jacobi
    start_counter();
    
    
    for(int iter = 0; iter < MAX_ITER; iter++)
    {
        norm2 = 0;
            
        /**
         * Paralelizamos bucle interno
         * schedule(runtime) -> Deseamos introducir por comandos el tipo de schedule a usar
         * reduction(+ : norm2) -> Privatizamos la variable norm2, haciendo que cada hilo calcule una parte y se sumen todas
         *                         al final.
         * critical: para pruebas, sólo 1 hilo a la vez trabaja en la región definida con critical. No es lo más eficiente
         */
        #pragma omp parallel num_threads(TH)
        {
            if(CRIT)
            {
                #pragma omp for schedule(runtime)
                for(int i = 0; i < N; i++)
                {
                    double sigma = 0.0;
                    int j;

                    //Desenrollamiento de lazos
                    for(j = 0; j < N-4; j+=4){
                        if(i != j)  sigma += a[i * N + j] * x[j];
                        if(i != (j+1))  sigma += a[i * N + j + 1] * x[j + 1];
                        if(i != (j+2))  sigma += a[i * N + j + 2] * x[j + 2];
                        if(i != (j+3))  sigma += a[i * N + j + 3] * x[j + 3];               
                    }

                    for(; j < N; j++){
                        if(i != j)   sigma += a[i * N + j] * x[j];
                    }

                    x_new[i] = (b[i] - sigma) / a[i * N + i];

                    if(CRIT > 1)
                    {
                        #pragma omp atomic
                        norm2 += (x_new[i] - x[i]) * (x_new[i] - x[i]);
                    }
                    else
                    {
                        #pragma omp critical
                        {
                            norm2 += (x_new[i] - x[i]) * (x_new[i] - x[i]);
                        }
                    }

                }
            }
            else
            {
                #pragma omp for schedule(runtime) reduction(+ : norm2)
                for(int i = 0; i < N; i++)
                {
                    double sigma = 0.0;
                    int j;

                    //Desenrollamiento de lazos
                    for(j = 0; j < N-4; j+=4){
                        if(i != j)  sigma += a[i * N + j] * x[j];
                        if(i != (j+1))  sigma += a[i * N + j + 1] * x[j + 1];
                        if(i != (j+2))  sigma += a[i * N + j + 2] * x[j + 2];
                        if(i != (j+3))  sigma += a[i * N + j + 3] * x[j + 3];               
                    }

                    for(; j < N; j++){
                        if(i != j)   sigma += a[i * N + j] * x[j];
                    }

                    x_new[i] = (b[i] - sigma) / a[i * N + i];

                    norm2 += (x_new[i] - x[i]) * (x_new[i] - x[i]);
                }
            }
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