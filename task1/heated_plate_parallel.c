#define _POSIX_C_SOURCE 199309L
#include "parallel_for.h"
#include <math.h>
#include <time.h>

#define M 500
#define N 500
#define THREAD_COUNT 4
#define max(a, b) ((a) > (b) ? (a) : (b))

typedef struct {
    double *u;
    double *w;
    double diff;            //记录最大差值（不断上升，初始化为0）
} MatrixArg;

pthread_mutex_t mutex;

void* update_temperature(void *arg, int row);
void print_matrix(double* u);           // for debug，如果要查看矩阵输出请修改M和N为20（方便打印），关闭第89行注释

int main() {
    // u是上一轮的副本，w为当前迭代的温度。
    double* u = (double*)malloc(M * N * sizeof(double));
    double* w = (double*)malloc(M * N * sizeof(double));
    double epsilon = 0.001, diff = epsilon, mean = 0.0;
    int iterations = 0, iterations_print = 1;
    
    pthread_mutex_init(&mutex, NULL);

    printf ( "\n" );
    printf ( "HEATED_PLATE_PTHREAD\n" );
    printf ( "  C/Pthread version\n" );
    printf ( "  A program to solve for the steady state temperature distribution\n" );
    printf ( "  over a rectangular plate.\n" );
    printf ( "\n" );
    printf ( "  Spatial grid of %d by %d points.\n", M, N );
    printf ( "  The iteration will be repeated until the change is <= %e\n", epsilon ); 
    printf ( "  Number of processors available = %d\n", THREAD_COUNT);
    printf ( "  Number of threads =              %d\n", THREAD_COUNT);

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
    printf ( "\n" );
    printf ( "  MEAN = %f\n", mean );

    for (int i = 1; i < M - 1; i++)
        for (int j = 1; j < N - 1; j++) 
            u[i*N + j] = w[i*N + j] = mean;

    printf ( "\n" );
    printf ( " Iteration  Change\n" );
    printf ( "\n" );

    // 记录开始时间
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    MatrixArg data;
    // 迭代计算温度分布
    while (diff >= epsilon) {
        //print_matrix(w);
        if (iterations % 2) 
            data = (MatrixArg){w, u, 0.0};
        else 
            data = (MatrixArg){u, w, 0.0};

        // 更新温度值
        parallel_for(1, M - 1, 1, update_temperature, &data, THREAD_COUNT);
        
        diff = data.diff;

        iterations++;
        if ( iterations == iterations_print )
        {
            printf ( "  %8d  %f\n", iterations, diff );
            iterations_print = 2 * iterations_print;
            //print_matrix(w);
        }
        

    }
    // 记录结束时间
    clock_gettime(CLOCK_MONOTONIC, &end);
    double wtime = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf ( "\n" );
    printf ( "  %8d  %f\n", iterations, diff );
    printf ( "\n" );
    printf ( "  Error tolerance achieved.\n" );
    printf ( "  Wallclock time = %f\n", wtime );
    //print_matrix(w);

    printf ( "\n" );
    printf ( "HEATED_PLATE_Pthread:\n" );
    printf ( "  Normal end of execution.\n" );
    
    pthread_mutex_destroy(&mutex);
    free(u);
    free(w);
    return 0;
}


void* update_temperature(void *arg, int row) {
    MatrixArg *data = (MatrixArg *)arg;
    double* u = data->u, *w = data->w;
    double my_diff = 0.0;
    for (int j = 1; j < N - 1; j++) {
        w[row*N + j] = (u[(row-1)*N + j] + u[(row+1)*N + j] + u[row*N + j - 1] + u[row*N + j + 1]) / 4.0;
        my_diff = max(my_diff, fabs(w[row*N + j] - u[row*N + j]));
    }
    pthread_mutex_lock(&mutex);
    data->diff = max(data->diff, my_diff);
    pthread_mutex_unlock(&mutex);
    return NULL;
}

void print_matrix(double* w) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("%5.1lf ", w[i*N + j]);
        }
        printf("\n");
    }
}