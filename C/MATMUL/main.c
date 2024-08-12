#include <stdio.h>
#include <matrix.h>
#include <clock.h>

// C (m x n) = A (m x p) * B (p x n)
void matmul(size_t m, size_t n, size_t p, double **A, double **B, double **C) {
    // Initialization
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            C[i][j] = 0;
        }
    }

    // Accumulation
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < p; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main(int argc, char *argv[]) {
    int param_iters = 1;

    if (argc != 2) {
        printf("Usage: %s <n> \n", argv[0]);
        printf("  <n> is the desired test size.\n");
        return 1;
    }

    // Reads the test parameters from the command line
    int param_n = 0;
    sscanf(argv[1], "%d", &param_n);
    printf("- Input parameters\n");
    printf("n\t= %i\n", param_n);
    size_t rows = param_n, cols = param_n;

    // Allocates input/output resources
    double **in1_mat = new_matrix(rows, cols);
    double **in2_mat = new_matrix(rows, cols);
    double **out_mat = new_matrix(rows, cols);
    if (!in1_mat || !in2_mat || !out_mat) {
        printf("Error: not enough memory to run the test using n = %i\n", param_n);
        return 1;
    }

    // Initializes data
    rand_matrix(in1_mat, rows, cols);
    rand_matrix(in2_mat, rows, cols);

    // Calls to the corresponding function to perform the computation
    printf("- Executing test...\n");
    double time_start = getClock();
    // ================================================

    for (int iters = 0; iters < param_iters; iters++) {
        matmul(rows, cols, cols, in1_mat, in2_mat, out_mat);
    }

    // ================================================
    double time_finish = getClock();

    // Prints an execution report
    double checksum = checksum_matrix(out_mat, rows, cols);
    printf("time (s)= %.6f\n", time_finish - time_start);
    printf("size\t= %i\n", param_n);
    printf("chksum\t= %.0f\n", checksum);
    if (param_iters > 1)
        printf("iters\t= %i\n", param_iters);

    // Release allocated resources
    delete_matrix(in1_mat);
    delete_matrix(in2_mat);
    delete_matrix(out_mat);

    return 0;
}
