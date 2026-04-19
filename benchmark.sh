#!/bin/bash

make
echo "Running reference case..."
for i in {1..5}
do
    mpirun -np 1 ./lbm -n -e 0 --config cases/config-wing.txt
done
for process in 1 2 4 8 16
do
    echo "Running with ${process} processes..."
    for ex in 1 2 3 4 5 6
    do
        echo "Running example ${ex}..."
        for i in {1..5}
        do
            mpirun -np ${process} ./lbm -n -e ${ex} --config cases/config-wing.txt
        done
    done
done
#remove all generated files
rm -f output.raw ref.md5 ex_*.md5