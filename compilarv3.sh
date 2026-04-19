#!/bin/bash
# Solicitamos un nodo con 64 cores y 16 GB de memoria durante 15 minutos
# Cambiar no sbatch para menos tempo e memoria para que nos den acceso máis fácilmente. Facer unha prueba con pouco
#SBATCH -n 1 -c 64 -t 00:30:00 --mem=16G
# Ponemos nombre a nuestro trabajo para poder identificarlo.
# ATENCIÓN - Debes sustituir el NN por el número de equipo.
#SBATCH --job-name p1acg1013
#SBATCH -o RES_v3.txt
#SBATCH -e ERR_v3.txt

gcc v3.c -o v3 -O0 -lm -fopenmp

for N in {1250,2000,3200}
do
    for TH in {2,4,8,12,32,64}
    do
        for SC in {0..3}
        do
            for CRIT in {0..2}
            do
                ./v3 $N $TH $SC $CRIT
            done
        done
    done
done