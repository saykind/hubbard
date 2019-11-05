#include "qmc.h"
#include "sqrmtx.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <cstdio>
#include <cstring>


using namespace std;

int main(int argc, char* argv[]) {
    int nt = nt0, verbose = 0;
    double u = u0, beta = beta0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i],"-v")) {
            verbose = 1;
            continue;
        }
        sscanf(argv[i], "beta=%lf", &beta);
        sscanf(argv[i], "u=%lf", &u);
        sscanf(argv[i], "nt=%d", &nt);
    }
    omp_set_num_threads(nt);
    if (verbose) PRINT_NOW;
    double t0 = omp_get_wtime();
    SimQMC hubbard(L0, N0, beta, u);
    if (verbose) PRINT_TIME;
    if (verbose) cout << endl << "initialized" << endl;
    if (verbose) PRINT_NOW;

    t0 = omp_get_wtime();
    for (int m=0; m<MTH0; m++) {
        for (int l=0; l<L0; l++) {
            for (int i=0; i<NN0; i++)
                ++hubbard;
            hubbard.l_();
        }
        hubbard.g_();
        if (verbose) SHOW_TIME(m,MTH0,t0)
    }
    if (verbose) cout << endl << "thermalized" << endl;
    if (verbose) PRINT_NOW;

    t0 = omp_get_wtime();
    for (int m=0; m<M0; m++) {
        for (int l=0; l<L0; l++) {
            for (int i=0; i<NN0; i++)
                ++hubbard;
            hubbard.l_();
        }
        hubbard.g_();
        hubbard.n_();
        if (verbose) SHOW_TIME(m,M0,t0)
    }
    if (verbose) cout << endl << "simulated" << endl;
    if (verbose) PRINT_NOW;
    cout << hubbard << endl;
    return 0;
}

