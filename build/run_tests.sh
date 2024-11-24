#!/usr/bin/bash

file=tests.csv

echo nthreads, size, time > $file

trials=5
size=32000000
nthread=(1 2 4 8)

for t in ${nthread[@]}; do
    echo 'with ' $t ' threads'
    for i in $(seq 1 $trials); do
        echo -ne $i 'out of' $trials 'runs\r'
        echo -ne $t ', ' $size ', ' >> $file
        ./main $t >> $file
    done
    echo -ne '\n------------\n'
done
echo ''
echo 'Done! Test results written to ' $file
