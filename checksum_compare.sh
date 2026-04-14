#!/bin/bash

make
mpirun -np 1 ./lbm -e 0
./display --gnuplot output.raw 3 | md5sum > ref.md5

for ex in 1 2 3 4 5 6
do
    for i in {1..10}
    do
        mpirun -np 4 ./lbm -e ${ex}
        ./display --gnuplot output.raw 3 | md5sum > ex_${ex}.md5
        diff ref.md5 ex_${ex}.md5 || echo "Error in example ${ex}"
    done
done

#remove all generated files
rm -f output.raw ref.md5 ex_*.md5