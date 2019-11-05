#pragma once
#define N0 6
#define NN0 (N0*N0)
#define L0 8
#define beta0 1.
#define dt0 (beta/L0)
#define u0 4.
#define lmbd0 acosh(exp(dt*U0/2)) 
#define MTH0 200
#define M0 1000
#define nt0 2

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif
#define pi M_PI

#define rank omp_get_thread_num()
#define NT omp_get_num_threads()

#define PRINT_NOW {time_t rawtime; time(&rawtime); printf("\n%s", asctime(localtime(&rawtime)));}
#define PRINT_TIME {double tt = omp_get_wtime()-t0; printf("  %.1lf min", tt/60);}
#define SHOW_TIME(i,M,t0) {double tt = omp_get_wtime()-t0, pc = 1.*(i+1)/M;\
        printf("  %.1lf min  \t%.0lf%%\test: %.1f min\tleft: %.1f min     \r", tt/60, 100.*pc, tt/pc/60, (1-pc)*(tt/pc)/60); fflush(stdout);}


