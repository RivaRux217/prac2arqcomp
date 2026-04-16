#!/bin/bash
# Solicitamos un nodo con 64 cores y 16 GB de memoria durante 15 minutos
# Cambiar no sbatch para menos tempo e memoria para que nos den acceso máis fácilmente. Facer unha prueba con pouco
#SBATCH -n 1 -c 64 -t 00:15:00 --mem=16G
# Ponemos nombre a nuestro trabajo para poder identificarlo.
# ATENCIÓN - Debes sustituir el NN por el número de equipo.
#SBATCH --job-name p1acg1013
#SBATCH -o RES_v1.txt
#SBATCH -e ERR_v1.txt

# Sustituir los valores de Di y Li por los calculados para la realización de la práctica.

gcc v1.c -o v1 -O0 -lm

for N in {1250,2000,3200}
do
    echo -e "\nVALOR DE N $N\n"
    for i in {1..10}
    do
        echo -e "\n\tEXPERIMENTO $i\n"
        ./v1 $N
    done
done

