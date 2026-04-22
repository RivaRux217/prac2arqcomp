#!/bin/bash
# Solicitamos un nodo con 64 cores y 16 GB de memoria durante 15 minutos
# Cambiar no sbatch para menos tempo e memoria para que nos den acceso máis fácilmente. Facer unha prueba con pouco
#SBATCH -n 1 -c 64 -t 00:30:00 --mem=16G
# Ponemos nombre a nuestro trabajo para poder identificarlo.
# ATENCIÓN - Debes sustituir el NN por el número de equipo.
#SBATCH --job-name p1acg1013
#SBATCH -o RES_v21.txt
#SBATCH -e ERR_v21.txt

# Sustituir los valores de Di y Li por los calculados para la realización de la práctica.

gcc v21.c -o v21 -O3 -lm

echo -e "1,0,0,0\n"

for N in {1250,2000,3200}
do
    echo -e "\nVALOR DE N $N\n"
    for i in {1..10}
    do
        echo -e "\n\tEXPERIMENTO $i\n"
        ./v21 $N 1
    done
done

echo -e "\n0,1,0,0\n"
for N in {1250,2000,3200}
do
    echo -e "\nVALOR DE N $N\n"
    for i in {1..10}
    do
        echo -e "\n\tEXPERIMENTO $i\n"
        ./v21 $N 0 1
    done
done

echo -e "\n0,0,1,0\n"
for N in {1250,2000,3200}
do
    echo -e "\nVALOR DE N $N\n"
    for i in {1..10}
    do
        echo -e "\n\tEXPERIMENTO $i\n"
        ./v21 $N 0 0 1
    done
done

echo -e "\n0,0,0,1\n"
for N in {1250,2000,3200}
do
    echo -e "\nVALOR DE N $N\n"
    for i in {1..10}
    do
        echo -e "\n\tEXPERIMENTO $i\n"
        ./v21 $N 0 0 0 1
    done
done

echo -e "\n1,1,0,0\n"
for N in {1250,2000,3200}
do
    echo -e "\nVALOR DE N $N\n"
    for i in {1..10}
    do
        echo -e "\n\tEXPERIMENTO $i\n"
        ./v21 $N 1 1
    done
done

echo -e "\n1,0,1,0\n"
for N in {1250,2000,3200}
do
    echo -e "\nVALOR DE N $N\n"
    for i in {1..10}
    do
        echo -e "\n\tEXPERIMENTO $i\n"
        ./v21 $N 1 0 1
    done
done

echo -e "\n1,0,0,1\n"
for N in {1250,2000,3200}
do
    echo -e "\nVALOR DE N $N\n"
    for i in {1..10}
    do
        echo -e "\n\tEXPERIMENTO $i\n"
        ./v21 $N 1 0 0 1
    done
done

echo -e "\n1,0,1,1\n"
for N in {1250,2000,3200}
do
    echo -e "\nVALOR DE N $N\n"
    for i in {1..10}
    do
        echo -e "\n\tEXPERIMENTO $i\n"
        ./v21 $N 1 0 1 1
    done
done

echo -e "\n0,0,1,1\n"
for N in {1250,2000,3200}
do
    echo -e "\nVALOR DE N $N\n"
    for i in {1..10}
    do
        echo -e "\n\tEXPERIMENTO $i\n"
        ./v21 $N 0 0 1 1
    done
done



