#!/bin/bash

for b1 in 0 1; do
  for b2 in 0 1; do
    for b3 in 0 1; do
      for b4 in 0 1; do

        # condición de exclusión:
        if [[ "$b4" -eq 1 && ( "$b2" -eq 1 || "$b3" -eq 1 ) ]]; then
          continue
        fi

        echo -e "\n$b1,$b2,$b3,$b4\n"

        for N in 1250 2000 3200; do
          echo -e "\nVALOR DE N $N\n"
          for i in {1..10}; do
            echo -e "\n\tEXPERIMENTO $i\n"
            ./v2 $N $b1 $b2 $b3 $b4
          done
        done

      done
    done
  done
done
