#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define M 500
#define N 500
#define max(a, b) ((a) > (b) ? (a) : (b))

double update_temperature(double *u, double *w, int rows_per_proc, int rank, int size, MPI_Comm comm);
double calculate_diff(double *u, double *w, int rows_per_proc);
void print_matrix(double* u);           // for debug，如果要查看矩阵输出请修改M和N为20（方便打印）

int main(int argc, char *argv[]) {
    int rank, size, rows_per_proc, start_row, end_row;
    double *u, *w, *local_u, *local_w;
    double epsilon = 0.001, diff = epsilon, mean = 0.0;
    int iterations = 0, iterations_print = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Barrier(MPI_COMM_WORLD); 
    double start_time = MPI_Wtime(); 

    if (rank == 0) {
        printf ( "\n" );
        printf ( "HEATED_PLATE_OPENMP\n" );
        printf ( "  C/OpenMP version\n" );
        printf ( "  A program to solve for the steady state temperature distribution\n" );
        printf ( "  over a rectangular plate.\n" );
        printf ( "\n" );
        printf ( "  Spatial grid of %d by %d points.\n", M, N );
        printf ( "  The iteration will be repeated until the change is <= %e\n", epsilon ); 
        printf ( "  Number of processors available = %d\n", size );
        printf ( "  Number of process =              %d\n", size );
    }
    

    rows_per_proc = M / size;
    start_row = rank * rows_per_proc;
    end_row = (rank == size - 1) ? M : start_row + rows_per_proc;

    // 分配全局和局部矩阵
    if (rank == 0) {
        u = (double *)malloc(M * N * sizeof(double));
        w = (double *)malloc(M * N * sizeof(double));

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
        mean /= (2 * N + 2 * M - 4);

        if (rank == 0)
            printf ( "\n  MEAN = %f\n", mean );

        for (int i = 1; i < M - 1; i++)
            for (int j = 1; j < N - 1; j++)
                u[i*N + j] = w[i*N + j] = mean;
    }

    // 为每个进程分配局部矩阵（包含上下邻居的额外一行）
    local_u = (double *)malloc((rows_per_proc + 2) * N * sizeof(double));
    local_w = (double *)malloc((rows_per_proc + 2) * N * sizeof(double));
    memset(local_u, 0, (rows_per_proc + 2) * N * sizeof(double));
    memset(local_w, 0, (rows_per_proc + 2) * N * sizeof(double));

    // 广播初始条件
    MPI_Scatter(w, rows_per_proc * N, MPI_DOUBLE, &local_w[N], rows_per_proc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(u, rows_per_proc * N, MPI_DOUBLE, &local_u[N], rows_per_proc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) 
        printf ( "\n Iteration  Change\n" );
    
    // 开始迭代
    double local_diff = 0;
    while (diff >= epsilon) {

        // 更新温度值
        if (iterations % 2)
            local_diff = update_temperature(local_w, local_u, rows_per_proc, rank, size, MPI_COMM_WORLD);
        else 
            local_diff = update_temperature(local_u, local_w, rows_per_proc, rank, size, MPI_COMM_WORLD);

        MPI_Allreduce(&local_diff, &diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        iterations++;
        
        if (rank == 0 && iterations == iterations_print) {
            printf("  %8d  %f\n", iterations, diff);
            iterations_print *= 2;
            //print_matrix(w);
        }
    } 

    // 汇总结果到根进程
    MPI_Gather(&local_w[N], rows_per_proc * N, MPI_DOUBLE, w, rows_per_proc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD); 
    double end_time = MPI_Wtime(); 

    // 输出最终结果
    if (rank == 0) {
        printf ( "\n" );
        printf ( "  %8d  %f\n", iterations, diff );
        printf ( "\n" );
        printf ( "  Error tolerance achieved.\n" );
        printf ( "  Wallclock time = %f\n", end_time - start_time );
        //print_matrix(w);

        printf ( "\n" );
        printf ( "HEATED_PLATE_Pthread:\n" );
        printf ( "  Normal end of execution.\n" );
    }

    // 释放资源
    free(local_u);
    free(local_w);
    if (rank == 0) {
        free(u);
        free(w);
    }

    MPI_Finalize();
    return 0;
}

double update_temperature(double *u, double *w, int rows_per_proc, int rank, int size, MPI_Comm comm) {

    // 与上/下邻居交换边界行
    if (rank > 0) 
        MPI_Sendrecv(&u[N], N, MPI_DOUBLE, rank - 1, 0, &u[0], N, MPI_DOUBLE, rank - 1, 0, comm, MPI_STATUS_IGNORE);

    if (rank < size - 1) 
        MPI_Sendrecv(&u[rows_per_proc * N], N, MPI_DOUBLE, rank + 1, 0, &u[(rows_per_proc + 1) * N], N, MPI_DOUBLE, rank + 1, 0, comm, MPI_STATUS_IGNORE);
    
    double diff = 0.0;
    for (int i = 1; i < rows_per_proc + 1; i++) {
        if (i == 1 && rank == 0 || i == rows_per_proc && rank == size - 1)
            continue;
        
        for (int j = 1; j < N - 1; j++) {
            w[i * N + j] = (u[(i - 1) * N + j] + u[(i + 1) * N + j] + u[i * N + j - 1] + u[i * N + j + 1]) / 4.0;
            diff = max(diff, fabs(w[i * N + j] - u[i * N + j]));
        }
    }
    return diff;
}

double calculate_diff(double *u, double *w, int rows_per_proc) {
    double diff = 0.0;
    for (int i = 1; i < rows_per_proc + 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            diff = max(diff, fabs(w[i * N + j] - u[i * N + j]));
        }
    }
    return diff;
}

void print_matrix(double* w) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("%5.1lf ", w[i*N + j]);
        }
        printf("\n");
    }
}