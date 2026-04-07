#!/bin/bash

make
mpirun -np 1 ./lbm -e 0
./gen_animate_gif.sh