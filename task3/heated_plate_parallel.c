#define _POSIX_C_SOURCE 199309L
#include "parallel_for.h"
#include <math.h>
#include <time.h>
#include <string.h>

#define max(a, b) ((a) > (b) ? (a) : (b))

typedef struct {
    double *u;
    double *w;
    double diff;            //记录最大差值（不断上升，初始化为0）
    int M, N;
} MatrixArg;

pthread_mutex_t mutex;

void* update_temperature(void *arg, int row);

int main(int argc, char* argv[]) {
    // u是上一轮的副本，w为当前迭代的温度。
    int M, N, THREAD_COUNT;
    if (argc != 4) {
        printf("Please input M, N and THREAD_COUNT!\n");
        exit(1);
    }
    M = atoi(argv[1]);
    N = atoi(argv[2]);
    THREAD_COUNT = atoi(argv[3]);
    double* u = (double*)malloc(M * N * sizeof(double));
    double* w = (double*)malloc(M * N * sizeof(double));
    memset(u, 0, M * N * sizeof(double));
    memset(w, 0, M * N * sizeof(double));
    double epsilon = 0.001, diff = epsilon, mean = 0.0;
    int iterations = 0;
    
    pthread_mutex_init(&mutex, NULL);

    // 初始化边界值，并计算均值
    for (int i = 1; i < M - 1; i++) {
        u[i * N + 0] = w[i * N + 0] = 100.0;
        u[i * N + N - 1] = w[i * N + N - 1] = 100.0;
        mean += w[i * N + 0] + w[i * N + N - 1];
    }
    for (int j = 0; j < N; j++) {
        u[0 * N + j] = w[0 * N + j] = 0.0;
        u[(M - 1) * N + j] = w[(M - 1) * N + j] = 100.0;
        mean += w[0 * N + j] + w[(M - 1) * N + j];
    }

    mean /= (2*N + 2*M - 4);
    for (int i = 1; i < M - 1; i++)
        for (int j = 1; j < N - 1; j++) 
            u[i*N + j] = w[i*N + j] = mean;

    // 记录开始时间
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    MatrixArg data;
    // 迭代计算温度分布
    while (diff >= epsilon) {
        //print_matrix(w);
        if (iterations % 2) 
            data = (MatrixArg){w, u, 0.0, M, N};
        else 
            data = (MatrixArg){u, w, 0.0, M, N};

        // 更新温度值
        parallel_for(1, M - 1, 1, update_temperature, &data, THREAD_COUNT);
        
        diff = data.diff;

        iterations++;

    }
    // 记录结束时间
    clock_gettime(CLOCK_MONOTONIC, &end);
    double wtime = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    
    printf ( "  Wallclock time = %f\n", wtime );
    //print_matrix(w);
    
    pthread_mutex_destroy(&mutex);
    free(u);
    free(w);
    return 0;
}


void* update_temperature(void *arg, int row) {
    MatrixArg *data = (MatrixArg *)arg;
    double* u = data->u, *w = data->w;
    double my_diff = 0.0;
    int M = data->M, N = data->N;
    for (int j = 1; j < N - 1; j++) {
        w[row*N + j] = (u[(row-1)*N + j] + u[(row+1)*N + j] + u[row*N + j - 1] + u[row*N + j + 1]) / 4.0;
        my_diff = max(my_diff, fabs(w[row*N + j] - u[row*N + j]));
    }
    pthread_mutex_lock(&mutex);
    data->diff = max(data->diff, my_diff);
    pthread_mutex_unlock(&mutex);
    return NULL;
}
