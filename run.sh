#!/bin/bash

make

mpirun -np 1 ./lbm -e 0 --config cases/config-wing.txt
./gen_animate_gif.sh
mv output.gif ref.gif

for i in {1..6}
do
    mpirun -np 4 ./lbm -e $i --config cases/config-wing.txt
    ./gen_animate_gif.sh
    mv output.gif output_$i.gif
done