#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <CRSMatrix.h>
#include <Matrix2D.h>
#include <Vector.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double getClock();

// Compute sparse matrix-vector multiplication
void atmux(double *val, double *x, double *y, int *col_ind, int *row_ptr, int n) {
    for (int t = 0; t < n; t++)
        y[t] = 0;

    // y = A^T x
    for (int i = 0; i < n; i++) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
            y[col_ind[k]] = y[col_ind[k]] + x[i] * val[k];
        }
    }
}

int main(int argc, char *argv[]) {
    double param_sparsity = 0.66;
    int param_iters = 10;

    if (argc != 2) {
        printf("Usage: %s <n>\n", argv[0]);
        printf("  <n> is the desired test size.\n");
        return 0;
    }

    // Reads the test parameters from the command line
    unsigned long param_n = 0;
    sscanf(argv[1], "%lu", &param_n);
    printf("- Input parameters\n");
    printf("size\t= %lu\n", param_n);

    // Allocates input/output resources and initializes data (if needed)
    Vector *out_vec = Vector_new(param_n);
    Vector *in_vec = Vector_new(param_n);
    Vector_rand(in_vec);
    Matrix2D *denseMat = Matrix2D_new(param_n, param_n);
    Matrix2D_randSparse(denseMat, param_sparsity);
    CRSMatrix *in_sparseMat = CRSMatrix_from(denseMat);

    if (!in_vec || !out_vec || !denseMat || !in_sparseMat) {
        printf("Error: not enough memory to run the test using n = %lu\n", param_n);
        return 0;
    }

    // Calls the corresponding function to perform the computation
    printf("- Executing test...\n");
    double time_start = getClock();
    // ================================================

    for (int iters = 0; iters < param_iters; iters++) {
        atmux(CRSMatrix_getData(in_sparseMat), Vector_getData(in_vec), Vector_getData(out_vec),
              CRSMatrix_colRef(in_sparseMat), CRSMatrix_rowRef(in_sparseMat), param_n);
    }

    // ================================================
    double time_finish = getClock();

    // Prints execution report
    double checksum = Vector_checksum(out_vec);
    printf("time (s)= %.6f\n", time_finish - time_start);
    printf("size\t= %lu\n", param_n);
    printf("sparsity= %g\n", param_sparsity);
    printf("chksum\t= %.0f\n", checksum);
    printf("iters\t= %i\n", param_iters);

    // Release allocated resources
    Matrix2D_delete(denseMat);
    CRSMatrix_delete(in_sparseMat);
    Vector_delete(in_vec);
    Vector_delete(out_vec);

    return 0;
}

double getClock() {
#ifdef _OPENMP
    return omp_get_wtime();
#elif __linux__ || __APPLE__
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1.0e9;
#else
    // Warning: this clock is invalid for parallel applications
    return (double)clock() / CLOCKS_PER_SEC;
#endif
}
