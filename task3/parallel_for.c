#include "parallel_for.h"

void* thread_function(void* arg) {
    ThreadData* data = (ThreadData*) arg;
    for (int i = data->start + data->thread_id * data->increment; i < data->end; i += data->increment * data->num_threads) {
        data->functor(data->arg, i);
    }
    return NULL;
}

void parallel_for(int start, int end, int increment, void *(*functor)(void*, int), void *arg, int num_threads) {
    pthread_t* threads = (pthread_t*) malloc(num_threads * sizeof(pthread_t));
    ThreadData* thread_data = (ThreadData*) malloc(num_threads * sizeof(ThreadData));

    // 初始化线程
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].start = start;
        thread_data[i].end = end;
        thread_data[i].increment = increment;
        thread_data[i].functor = functor;
        thread_data[i].arg = arg;
        thread_data[i].thread_id = i;
        thread_data[i].num_threads = num_threads;
        pthread_create(&threads[i], NULL, thread_function, &thread_data[i]);
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    free(threads);
    free(thread_data);
}
