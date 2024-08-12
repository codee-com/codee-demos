#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>

#include "Vector.h"
#include "Matrix2D.h"

#ifdef _OPENMP
#include <omp.h>
#endif

double getClock();

void coulomb(double* vec, int size,  // List of charged particles (xyzq)
	double* mat, int rows, int cols, // Output matrix
	double x0, double y0, double z0, // Initial point
	double x1, double y1)            // Final point
{
	const double PI = 3.14159265358979324;
	const double e0 = 8.854187817e-12;
	double scaleX = (x1 - x0) / cols;
	double scaleY = (y1 - y0) / rows;
	
	for(int i=0; i<rows; i++) {
		for(int j=0; j<cols; j++) {
			double mat_ij = 0;
			for(int k = 0; k < size; k+=4) {
				double dx = vec[k+0] - (scaleX * j + x0);
				double dy = vec[k+1] - (scaleY * i + y0);
				double dz = vec[k+2] - (z0);
				double charge = 1e-9 * vec[k+3];
				double dist = sqrt(dx*dx + dy*dy + dz*dz);
				mat_ij += charge / dist;
			}
        	mat[j + i * cols] = mat_ij / (4 * PI * e0);
		}
	}
}

void chargePrinter(int pos, double* charge) {
	printf("%i> Charge at (x=%.0f, y=%.0f, z=%.0f) is %.0f nC\n",
		pos, charge[0], charge[1], charge[2], charge[3]);
}

int main(int argc, char* argv[]) {
	// Reads the test parameters from the command line
	double arg_n = 0.0, arg_density = 0.1;
	int param_iters = 1;
	if(argc >= 2) sscanf(argv[1], "%lf", &arg_n);
	if(argc >= 3) param_iters = atoi(argv[2]);
	if(argc >= 4) sscanf(argv[3], "%lf", &arg_density);
	
	if(arg_n < 1 || arg_n >= INT_MAX || param_iters < 1 || arg_density < 0.0 || arg_density > 1.0) {
		printf("This test computes the electric potential created\n");
		printf("by a set of charges in an n x n 2D plane.\n");
		printf("  The first parameter <n> is the desired test size.\n");
		printf("  The optional parameter [iters] is used to repeat the test several times.\n");
		printf("  The optional parameter [density] is the ratio of charges in the plane.\n");
		printf("Usage: %s <n> [iters] [density] \n", argv[0]);
		exit(0);
	}

	// Allocates input/output resources
	int param_n = (int)arg_n;
	int numCharges = arg_n * arg_n * arg_density + 0.5;
	Matrix2D* out_mat = Matrix2D_new(param_n, param_n);
	Vector* in_vec = Vector_new(4 * numCharges);
	
	if(!in_vec || !out_mat) {
		if(numCharges == 0) printf("Error: There are no charges in the domain\n");
		else printf("Error: Not enough memory to run the test using n = %i\n", param_n);
		exit(0);
	}

	// Initializes data if needed
	Vector_rand(in_vec);
		
	// Calls the function that performs the actual computation
	printf("- Executing test...\n");
	double time_start = getClock();
	for(int iters = 0; iters < param_iters; iters++) {
		coulomb(
			Vector_getData(in_vec), Vector_getSize(in_vec),
			Matrix2D_getData(out_mat)[0], param_n, param_n,
			0, 0, 0, param_n, param_n
		);
	}
	double time_finish = getClock();

	// Prints an execution report
	double checksum = Matrix2D_checksum(out_mat);
	printf("time (s)= %.6f\n", time_finish - time_start);
	printf("size\t= %i\n", param_n);
	printf("chksum\t= %.0f\n", checksum);
	if(param_iters > 1) printf("iters\t= %i\n", param_iters);

	if(param_n < 9) { // Show example for small problems
		printf("\n- Input vector b:\n");
		for(int i = 0; i < Vector_getSize(in_vec); i+=4)
			chargePrinter(i>>2, Vector_getData(in_vec)+i);
		printf("\n- Output matrix A:\n");
		Matrix2D_print(out_mat);
	}
	printf("\n");
	// Release allocated resources
	Vector_delete(in_vec);
	Matrix2D_delete(out_mat);
	
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
