#! /bin/bash
#

gcc -fPIC -shared -o libparallel_for.so parallel_for.c -lpthread 
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH 

gcc -O3 -o heated_plate_parallel heated_plate_parallel.c -L. -lparallel_for -lpthread

echo "heated_plate_parallel compiled successfully!"
echo "./heated_plate_parallel M N THREAD_COUNT"