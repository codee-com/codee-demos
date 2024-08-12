#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double getClock();

#define LEFTVAL 1.0
#define RIGHTVAL 10.0

double calculate_checksum(double *A, int n);


void compute(double *out_error, double *uk, double *T, int *T_mapping, int nx) {
    double maxdiff = 0;

    // Compute new values
    double dx = 1.0 / nx;
    double dt = 0.5 * dx * dx;
    for (int i = 1; i < nx - 1; ++i) {
        T[i] = 0;
    }
    for (int i = 1; i < nx - 1; ++i) {
        T[T_mapping[i]] += uk[i] + (dt / (dx * dx)) * (uk[i + 1] - 2 * uk[i] + uk[i - 1]);
    }

    // Check for convergence
    for (int i = 1; i < nx - 1; ++i) {
        double diff = fabs(uk[i] - T[i]);
        maxdiff = fmax(maxdiff, diff);
    }

    // Copy ukp1 to uk
    for (int i = 0; i < nx; i++) {
        uk[i] = T[i];
    }

    *out_error = maxdiff;
}

void initialize_uk(double uk[], int nx) {
    uk[0] = LEFTVAL;
    uk[nx - 1] = RIGHTVAL;
    for (int i = 1; i < nx - 1; ++i)
        uk[i] = 0.0;
}

void initialize_T(double ukp1[], int ukp1_mapping[], double uk[], int nx) {
    for (int i = 0; i < nx; ++i)
        ukp1[i] = uk[i];
    for (int i = 0; i < nx; ++i)
        ukp1_mapping[i] = i;
}

void cfd_heat_diffusion(double *uk, double threshold, int maxsteps, int nx, int *out_step, double *out_maxdiff) {
    int step;
    double maxdiff = threshold;

    double *T = malloc(sizeof(*T) * nx);
    int *T_mapping = malloc(sizeof(*T_mapping) * nx);
    if (!T || !T_mapping) {
        printf("Error: not enough memory for temporaries to run the test using n = %i\n", nx);
        return;
    }
    initialize_T(T, T_mapping, uk, nx);

    for (step = 0; (step < maxsteps) && (maxdiff >= threshold); ++step) {
        compute(&maxdiff, uk, T, T_mapping, nx);
    }

    *out_maxdiff = maxdiff;
    *out_step = step;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s <n>\n", argv[0]);
        printf("  <n> is the desired test size.\n");
        return 1;
    }

    // Reads the test parameters from the command line
    int param_n = atoi(argv[1]);
    int param_iters = 5000;
    double threshold = 0.001;

    // Allocates resources
    double *uk = malloc(sizeof(*uk) * param_n);
    if (!uk) {
        printf("Error: not enough memory to run the test using n = %i\n", param_n);
        return 1;
    }
    initialize_uk(uk, param_n);

    int iters;
    double maxdiff;

    printf("- Executing test...\n");
    double time_start = getClock();
    // ================================================

    cfd_heat_diffusion(uk, threshold, param_iters, param_n, &iters, &maxdiff);

    // ================================================
    double time_finish = getClock();
    double checksum = calculate_checksum(uk, param_n);

    printf("size = %d, iters = %d, maxiters = %d, maxdiff = %g, threshold = %g\n", param_n, iters, param_iters, maxdiff,
           threshold);
    if (maxdiff < threshold) {
        printf("converged in %d iterations\n", iters);
    } else {
        printf("failed to converge in %d iterations, maxdiff = %g\n", iters, maxdiff);
    }
    printf("time (s) = %.6f\n", time_finish - time_start);
    printf("checksum = %.6f\n", checksum);

    return 0;
}

// Generates a checksum based on the output data
double calculate_checksum(double *A, int n) {
    if (!A)
        return 0;

    double checkSum = 0.0;
    for (int row = 0; row < n; row++)
        checkSum += A[row];
    return checkSum;
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
