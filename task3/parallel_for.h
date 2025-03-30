#ifndef PARALLEL_FOR
#define PARALLEL_FOR

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct {
    int start;
    int end;
    int increment;
    void *(*functor)(void*, int);
    void *arg;
    int thread_id;
    int num_threads;
} ThreadData;


void* thread_function(void* arg);

void parallel_for(int start, int end, int increment, void *(*functor)(void*, int), void *arg, int num_threads);

#endif