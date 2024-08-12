/*
This is an oversimplification of the Livermore Unstructured Lagrangian
Explicit Shock Hydrodynamics (LULESH) program structure without the
workload functionality.

More information on LULESH:
  - https://computing.llnl.gov/projects/co-design/lulesh
  - https://github.com/LLNL/LULESH

                 Copyright (c) 2010-2013.
      Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
                  LLNL-CODE-461231
                All rights reserved.

                 Copyright (c) 2020
               Appentra Solutions S.L.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define WORKLOAD_CalcElemFBHourglassForce 12             // Adjusted to 58% of total runtime
#define WORKLOAD_ApplyMaterialPropertiesForElems 1000000 // Adjusted to 20% of total runtime
#define WORKLOAD_CalcElemVelocityGradient 25             // Adjusted to 15% of total runtime

#define MAX_TIME_STEPS 932
#define NUM_NODES 27000
#define NUM_ELEMS 30000
#define MAX_NODELIST 240000 // LULESH: 8 * NUM_ELEMS = 8 * 30000
#define NUM_REGIONS 1

#define MAX_TOLERANCE 1.0e-5

double getClock();
double calculate_checksum(double *A, int n);


// Precision specification
typedef double real8;

typedef int Index_t;  // array subscript and loop index
typedef real8 Real_t; // floating point representation
typedef int Int_t;    // integer representation

enum { VolumeError = -1, QStopError = -2 };

typedef struct Parameters {
    int param_iters;
    double param_tol;

    int nvertices;
    int nelements;

    double *A;
    double *T;
    int *nodes;
    double *energy;
} Parameters;

Parameters Parameters_create();
void Parameters_free(Parameters p);

void VerifyAndWriteFinalOutput(Real_t elapsed_time, Int_t nx, Int_t numRanks, Int_t locDom_cycle, Real_t *locDom_e,
                               Int_t locDom_numNode, Int_t locDom_numElem, Int_t locDom_numReg,
                               Int_t locDom_regElemSize[], Real_t *locDom_f);


// NOTE: We use the computation of PI only to represent an integration method
// and add workload artificially;
double CalcElemFBHourglassForce_workload(int N) {

    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double x = (i + 0.5) / N;
        sum += sqrt(1 - x * x);
    }
    return 4.0 / N * sum;
}

void CalcElemFBHourglassForce(Index_t i2, Real_t gamma[4][8], Real_t *hgfx, Real_t *hgfy, Real_t *hgfz) {
    for (Index_t i = 0; i < 8; i++) {
        double T = CalcElemFBHourglassForce_workload(WORKLOAD_CalcElemFBHourglassForce);
        hgfx[i] = gamma[0][0] * T;
        hgfy[i] = gamma[0][1] * T;
        hgfz[i] = gamma[0][2] * T;
    }
}

void CalcFBHourglassForceForElems(Index_t numElem, Index_t *domain_m_nodelist, Real_t *domain_m_fx, Real_t *domain_m_fy,
                                  Real_t *domain_m_fz) {
    /*************************************************
     *
     *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
     *               force.
     *
     *************************************************/
    Real_t gamma[4][8];

    gamma[0][0] = (1.);
    gamma[0][1] = (1.);
    gamma[0][2] = (-1.);
    gamma[0][3] = (-1.);
    gamma[0][4] = (-1.);
    gamma[0][5] = (-1.);
    gamma[0][6] = (1.);
    gamma[0][7] = (1.);
    gamma[1][0] = (1.);
    gamma[1][1] = (-1.);
    gamma[1][2] = (-1.);
    gamma[1][3] = (1.);
    gamma[1][4] = (-1.);
    gamma[1][5] = (1.);
    gamma[1][6] = (1.);
    gamma[1][7] = (-1.);
    gamma[2][0] = (1.);
    gamma[2][1] = (-1.);
    gamma[2][2] = (1.);
    gamma[2][3] = (-1.);
    gamma[2][4] = (1.);
    gamma[2][5] = (-1.);
    gamma[2][6] = (1.);
    gamma[2][7] = (-1.);
    gamma[3][0] = (-1.);
    gamma[3][1] = (1.);
    gamma[3][2] = (-1.);
    gamma[3][3] = (1.);
    gamma[3][4] = (1.);
    gamma[3][5] = (-1.);
    gamma[3][6] = (1.);
    gamma[3][7] = (-1.);

    /*************************************************/
    /*    compute the hourglass modes */

    for (Index_t i2 = 0; i2 < numElem; ++i2) {
        Real_t hgfx[8], hgfy[8], hgfz[8];

        CalcElemFBHourglassForce(i2, gamma, hgfx, hgfy, hgfz);

        // With the threaded version, we write into local arrays per elem
        // so we don't have to worry about race conditions
        Index_t n0si2 = domain_m_nodelist[(8) * i2 + 0];
        Index_t n1si2 = domain_m_nodelist[(8) * i2 + 1];
        Index_t n2si2 = domain_m_nodelist[(8) * i2 + 2];
        Index_t n3si2 = domain_m_nodelist[(8) * i2 + 3];
        Index_t n4si2 = domain_m_nodelist[(8) * i2 + 4];
        Index_t n5si2 = domain_m_nodelist[(8) * i2 + 5];
        Index_t n6si2 = domain_m_nodelist[(8) * i2 + 6];
        Index_t n7si2 = domain_m_nodelist[(8) * i2 + 7];

        domain_m_fx[n0si2] += hgfx[0];
        domain_m_fy[n0si2] += hgfy[0];
        domain_m_fz[n0si2] += hgfz[0];

        domain_m_fx[n1si2] += hgfx[1];
        domain_m_fy[n1si2] += hgfy[1];
        domain_m_fz[n1si2] += hgfz[1];

        domain_m_fx[n2si2] += hgfx[2];
        domain_m_fy[n2si2] += hgfy[2];
        domain_m_fz[n2si2] += hgfz[2];

        domain_m_fx[n3si2] += hgfx[3];
        domain_m_fy[n3si2] += hgfy[3];
        domain_m_fz[n3si2] += hgfz[3];

        domain_m_fx[n4si2] += hgfx[4];
        domain_m_fy[n4si2] += hgfy[4];
        domain_m_fz[n4si2] += hgfz[4];

        domain_m_fx[n5si2] += hgfx[5];
        domain_m_fy[n5si2] += hgfy[5];
        domain_m_fz[n5si2] += hgfz[5];

        domain_m_fx[n6si2] += hgfx[6];
        domain_m_fy[n6si2] += hgfy[6];
        domain_m_fz[n6si2] += hgfz[6];

        domain_m_fx[n7si2] += hgfx[7];
        domain_m_fy[n7si2] += hgfy[7];
        domain_m_fz[n7si2] += hgfz[7];
    }
}


// NOTE: We use the computation of PI only to represent an integration method
// and add workload artificially;
double ApplyMaterialPropertiesForElems_workload(int N) {

    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double x = (i + 0.5) / N;
        sum += sqrt(1 - x * x);
    }
    return 4.0 / N * sum;
}

void ApplyMaterialPropertiesForElems(Index_t domain_m_numElem, Real_t domain_m_eosvmin, Real_t domain_m_eosvmax,
                                     Real_t vnew[]) {
    Index_t numElem = domain_m_numElem;

    if (numElem != 0) {
        /* Expose all of the variables needed for material evaluation */
        Real_t eosvmin = domain_m_eosvmin;
        Real_t eosvmax = domain_m_eosvmax;

        double T = ApplyMaterialPropertiesForElems_workload(WORKLOAD_ApplyMaterialPropertiesForElems);

        // Bound the updated relative volumes with eosvmin/max
        if (eosvmin != (0.)) {
            for (Index_t i = 0; i < numElem; ++i) {
                if (vnew[i] < eosvmin)
                    vnew[i] = eosvmin;
            }
        }

        if (eosvmax != (0.)) {
            for (Index_t i = 0; i < numElem; ++i) {
                if (vnew[i] > eosvmax)
                    vnew[i] = eosvmax;
            }
        }
    }
}


// NOTE: We use the computation of PI only to represent an integration method
// and add workload artificially;
double CalcElemVelocityGradient_workload(int N) {

    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double x = (i + 0.5) / N;
        sum += sqrt(1 - x * x);
    }
    return 4.0 / N * sum;
}


void CalcElemVelocityGradient(Index_t k, Real_t *const d) {
    double T = CalcElemVelocityGradient_workload(WORKLOAD_CalcElemVelocityGradient);
    d[0] = k * T;
    d[1] = k * T;
    d[2] = k * T;
    d[5] = k * T;
    d[4] = k * T;
    d[3] = k * T;
}


void CalcKinematicsForElems(Index_t numElem, Real_t *domain_m_dxx, Real_t *domain_m_dyy, Real_t *domain_m_dzz) {
    // loop over all elements
    for (Index_t k = 0; k < numElem; ++k) {
        Real_t D[6];

        CalcElemVelocityGradient(k, D);

        // put velocity gradient quantities into their global arrays.
        domain_m_dxx[k] = D[0];
        domain_m_dyy[k] = D[1];
        domain_m_dzz[k] = D[2];
    }
}

int luleshmk(double param_tol, Int_t maxsteps, Int_t numElem, Int_t numNodes, Index_t *domain_m_nodelist,
             Real_t *domain_m_fx, Real_t *domain_m_fy, Real_t *domain_m_fz, Real_t *domain_m_dxx, Real_t *domain_m_dyy,
             Real_t *domain_m_dzz, Int_t domain_m_eosvmin, Int_t domain_m_eosvmax, Real_t *vnew, Real_t *locDom_f,
             Real_t *locDom_e) {
    // Timestep loop of the simulation
    int step = 0;
    double err = param_tol;

    for (step = 0, err = param_tol; err >= param_tol && step < maxsteps; step++) {
        // ---------- Call graph structure of function main() ------------
        // In each time-step iteration, first invoke TimeIncrement()

        // In each time-step iteration, second LagrangeLeapFrog()
        // Step 1: LagrangeNodal()
        // 1.1 CalcForceForNodes(): output domain_m_fx, domain_m_fy, domain_m_fz
        // 1.1.1 CalcVolumeForceForElems()
        // 1.1.1.1 CalcHourglassControlForElems()
        // 1.1.1.1.1 CalcFBHourglassForceForElems() -HOTSPOT 30%-
        // 1.2 CalcAccelerationForNodes()
        // 1.3 ApplyAccelerationBoundaryConditionsForNodes()
        // 1.4 CalcVelocityForNodes()
        // 1.5 CalcPositionForNodes()

        // Step 2: LagrangeElements()
        // 2.1 malloc vnew[numElem]
        // 2.2 CalcLagrangeElements()
        // 2.2.1 CalcKinematicsForElems() -HOTSPOT 8%-
        // 2.3 CalcQForElems()
        // 2.4 ApplyMaterialPropertiesForElems() -HOTSPOT 20%-
        // 2.5 UpdateVolumesForElems()
        // 2.6 free(vnew)

        // Step 3: CalcTimeConstraintsForElems()
        // ---------- Call graph structure of function main() ------------

        // 1.1.1.1.1 CalcFBHourglassForceForElems() -HOTSPOT 30%-
        CalcFBHourglassForceForElems(numElem, domain_m_nodelist, domain_m_fx, domain_m_fy, domain_m_fz);
        // output domain_m_fx, domain_m_fy, domain_m_fz

        // 2.2.1 CalcKinematicsForElems() -HOTSPOT 8%-
        CalcKinematicsForElems(numElem, domain_m_dxx, domain_m_dyy, domain_m_dzz);
        // output domain_m_dxx, domain_m_dyy, domain_m_dzz

        // 2.4 ApplyMaterialPropertiesForElems() -HOTSPOT 20%-
        ApplyMaterialPropertiesForElems(numElem, domain_m_eosvmin, domain_m_eosvmax, vnew);
        // output vnew

        // Finally update the output "locDom_e" in order to allow for the
        // verification of the numerical results at VerifyAndWriteFinalOutput()
        for (int i = 0; i < numNodes; ++i) {
            locDom_f[i] += domain_m_fx[i] + domain_m_fy[i] + domain_m_fz[i];
        }
        for (int i = 0; i < numElem; ++i) {
            locDom_e[i] += domain_m_dxx[i] + domain_m_dyy[i] + domain_m_dzz[i] + vnew[i];
        }

#if DEBUG_ENABLED
        printf("step: %d , maxsteps: %d , maxdiff: %g ,  threshold: %g , checksum_f: %g , checksum_e: %g\n", step,
               maxsteps, err, param_tol, calculate_checksum(locDom_f, numNodes), calculate_checksum(locDom_e, numElem));
#endif
    }

    return step;
}

int main(int argc, char *argv[]) {
    if (argc != 1) {
        printf("Usage: %s\n", argv[0]);
        exit(0);
    }

    printf("- Configuring the test...\n");

    Parameters p = Parameters_create();

    Int_t locDom_cycle = MAX_TIME_STEPS;
    Int_t locDom_numNode = NUM_NODES;
    Int_t locDom_numElem = NUM_ELEMS;
    Int_t locDom_numReg = NUM_REGIONS; // LULESH processes each region independently, so 1 emulates this.
    Int_t locDom_regElemSize[NUM_REGIONS] = {locDom_numElem};

    Int_t locDom_m_eosvmin = 0;
    Int_t locDom_m_eosvmax = 0;

    Real_t locDom_e[NUM_ELEMS]; // Check if 900 is enouch = nx*nx with opt_nx=30 below
    Real_t locDom_m_dxx[NUM_ELEMS];
    Real_t locDom_m_dyy[NUM_ELEMS];
    Real_t locDom_m_dzz[NUM_ELEMS];
    Real_t vnew[NUM_ELEMS];
    for (int i = 0; i < NUM_ELEMS; ++i) {
        locDom_e[i] = i + 1.0;
        locDom_m_dxx[i] = i + 1.0;
        locDom_m_dyy[i] = i + 1.0;
        locDom_m_dzz[i] = i + 1.0;
        vnew[i] = i + 1.0;
    }

    Index_t locDom_m_nodelist[MAX_NODELIST];
    for (int i = 0; i < MAX_NODELIST; ++i) {
        locDom_m_nodelist[i] = i % NUM_NODES; // Indirections are bounded to NUM_NODES
    }

    Real_t locDom_f[NUM_NODES];    // Added in the miniapp for verification purposes
    Real_t locDom_m_fx[NUM_NODES]; // LULESH compute force at each node mesh, so 27000 on x
    Real_t locDom_m_fy[NUM_NODES]; // and 27000 on y
    Real_t locDom_m_fz[NUM_NODES]; // and 27000 on z
    for (int i = 0; i < NUM_NODES; ++i) {
        locDom_f[i] = 2.0;
        locDom_m_fx[i] = 2.0;
        locDom_m_fy[i] = 3.0;
        locDom_m_fz[i] = 4.0;
    }

    Int_t numRanks = 1;
    Int_t myRank = 0;

    /* Set defaults that can be overridden by command line opts */
    Int_t opts_its = 9999999;
    Int_t opts_nx = 30;
    Int_t opts_numReg = 11;
    Int_t opts_numFiles = (int)(numRanks + 10) / 9;
    Int_t opts_showProg = 0;
    Int_t opts_quiet = 0;
    Int_t opts_viz = 0;
    Int_t opts_balance = 1;
    Int_t opts_cost = 1;

    printf("- Executing the test...\n");
    double time_start = getClock();
    // ================================================

    p.param_iters = 932; // Max number of iteration of LULESH default use case
    luleshmk(p.param_tol, locDom_cycle, locDom_numElem, locDom_numNode, locDom_m_nodelist, locDom_m_fx, locDom_m_fy,
             locDom_m_fz, locDom_m_dxx, locDom_m_dyy, locDom_m_dzz, locDom_m_eosvmin, locDom_m_eosvmax, vnew, locDom_f,
             locDom_e);

    // ================================================
    double time_finish = getClock();
    double elapsed_timeG = time_finish - time_start;

    printf("- Verifying the test...\n");
    VerifyAndWriteFinalOutput(elapsed_timeG,
                              /* *locDom ,*/
                              opts_nx, numRanks, locDom_cycle, locDom_e, locDom_numNode, locDom_numElem, locDom_numReg,
                              locDom_regElemSize, locDom_f);

    Parameters_free(p);

    return 0;
}

Parameters Parameters_create() {

    Parameters p;

    p.param_iters = MAX_TIME_STEPS;
    p.param_tol = MAX_TOLERANCE;
    p.nelements = NUM_ELEMS;
    p.nvertices = NUM_NODES;

    p.A = (double *)calloc(p.nvertices, sizeof(double));
    if (!p.A) {
        printf("Error: not enough memory to allocate array <A> ");
        printf("of size nvertices=%i\n", p.nvertices);
        exit(0);
    }
    p.T = (double *)calloc(p.nvertices, sizeof(double));
    if (!p.T) {
        printf("Error: not enough memory to allocate array <T> ");
        printf("of size nvertices=%i\n", p.nvertices);
        exit(0);
    }
    p.nodes = (int *)calloc(p.nelements, sizeof(int));
    if (!p.nodes) {
        printf("Error: not enough memory to allocate array <nodes> ");
        printf("of size nelements=%i\n", p.nelements);
        exit(0);
    }
    p.energy = (double *)calloc(p.nelements, sizeof(double));
    if (!p.energy) {
        printf("Error: not enough memory to allocate array <energy> ");
        printf("of size nelements=%i\n", p.nelements);
        exit(0);
    }

    return p;
}

void Parameters_free(Parameters p) {
    free(p.A);
    free(p.T);
    free(p.nodes);
    free(p.energy);
}

void VerifyAndWriteFinalOutput(Real_t elapsed_time, Int_t nx, Int_t numRanks, Int_t locDom_cycle, Real_t *locDom_e,
                               Int_t locDom_numNode, Int_t locDom_numElem, Int_t locDom_numReg,
                               Int_t locDom_regElemSize[], Real_t *locDom_f) {
    // GrindTime1 only takes a single domain into account, and is thus a good way to measure
    // processor speed indepdendent of MPI parallelism.
    // GrindTime2 takes into account speedups from MPI parallelism
    Real_t grindTime1 = ((elapsed_time * 1e6) / locDom_cycle) / (nx * nx * nx);
    Real_t grindTime2 = ((elapsed_time * 1e6) / locDom_cycle) / (nx * nx * nx * numRanks);

    Index_t ElemId = 0;
    printf("Run completed:  \n");
    printf("   Problem size        =  %i \n", nx);
    printf("   MPI tasks           =  %i \n", numRanks);
    printf("   Iteration count     =  %i \n", locDom_cycle);
    printf("   Final Origin Energy = %12.6e \n", locDom_e[ElemId]);

    Real_t MaxAbsDiff = (0.0);
    Real_t TotalAbsDiff = (0.0);
    Real_t MaxRelDiff = (0.0);

    for (Index_t j = 0; j < nx; ++j) {
        for (Index_t k = j + 1; k < nx; ++k) {
            Real_t AbsDiff = fabs(locDom_e[j * nx + k] - locDom_e[k * nx + j]);
            TotalAbsDiff += AbsDiff;

            if (MaxAbsDiff < AbsDiff)
                MaxAbsDiff = AbsDiff;

            Real_t RelDiff = AbsDiff / locDom_e[k * nx + j];

            if (MaxRelDiff < RelDiff)
                MaxRelDiff = RelDiff;
        }
    }

#if DEBUG_ENABLED
    // Computational workload check
    printf("   Number of nodes    = %d \n", locDom_numNode);
    printf("   Number of elements = %d \n", locDom_numElem);
    printf("   Number of regions  = %d \n", locDom_numReg);
    for (Int_t i = 0; i < locDom_numReg; i++)
        printf("      Region %d of size %d \n", (i + 1), locDom_regElemSize[i]);
#endif

    // Quick symmetry check
    printf("   Testing Plane 0 of Energy Array on rank 0:\n");
    printf("        MaxAbsDiff   = %12.6e\n", MaxAbsDiff);
    printf("        TotalAbsDiff = %12.6e\n", TotalAbsDiff);
    printf("        MaxRelDiff   = %12.6e\n\n", MaxRelDiff);

    // Timing information
    printf("\nElapsed time         = %10.2f (s)\n", elapsed_time);
    printf("Grind time (us/z/c)  = %10.8g (per dom)  (%10.8g overall)\n", grindTime1, grindTime2);
    printf("FOM                  = %10.8g (z/s)\n\n", 1000.0 / grindTime2); // zones per second

    // Added in the miniapp to verify the correctness of the code
    // Verification of the numerical results at VerifyAndWriteFinalOutput()
    Real_t checksum_f = (0.0);
    for (int i = 0; i < locDom_numNode; ++i) {
        checksum_f += locDom_f[i];
    }
    Real_t checksum_e = (0.0);
    for (int i = 0; i < locDom_numElem; ++i) {
        checksum_e += locDom_e[i];
    }

    printf("\n");
    printf("numNodes   = %d\n", locDom_numNode);
    printf("numElems   = %d\n", locDom_numElem);
    printf("checksum_f = %g\n", calculate_checksum(locDom_f, locDom_numNode));
    printf("checksum_e = %g\n", calculate_checksum(locDom_e, locDom_numElem));

    return;
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
