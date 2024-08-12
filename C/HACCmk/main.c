#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
//#include <mpi.h>

extern void Step10_orig( int count1, float xxi, float yyi, float zzi, float fsrrmax2, float mp_rsm2, float *xx1, float *yy1, float *zz1, float *mass1, float *dxi, float *dyi, float *dzi );

#ifdef TIMEBASE
extern unsigned long long timebase();
#else
extern double mysecond();
#endif

// Frequency in MHz, L1 cache size in bytes, Peak MFlops per node */
// BG/Q 
#define MHz 1600e6
#define NC 16777216
#define PEAK (16.*12800.)

// 4 core Intel(R) Xeon(R) CPU E5-2643 0 @ 3.30GHz - Rick's machine celero 
// Peak flop rate: AVX: 8 Flops/cycle * 4 cores * 3291.838 Mcycles/second 
/*
#define MHz 3291.838e6     
#define NC (64*1024*1024)
#define PEAK 105338.816
*/

// 4 core Intel(R) Xeon(R) CPU           E5430  @ 2.66GHz - crush
// Peak flop rate: SSE: 4 Flops/cycle * 4 cores * 2666.666 Mcycles/second
/*
#define MHz 2666.7e6
#define NC ( 256*1024*1024 )
#define PEAK 42667.2
*/

// 4 core Intel(R) Core(TM) i7-3820 CPU @ 3.60GHz - strength
// Peak flop rate: AVX: 8 Flops/cycle * 4 cores * 3600 MCycles/second
/*
#define MHz 3600.e6
#define NC ( 32*1024*1024 )
#define PEAK 115200.
*/

#define N 15000      /* Vector length, must be divisible by 4  15000 */

#define ETOL 1.e-4  /* Tolerance for correctness */

char  M1[NC], M2[NC];

int main( int argc, char *argv[] )
{
  static float xx[N], yy[N], zz[N], mass[N], vx1[N], vy1[N], vz1[N];
  float fsrrmax2, mp_rsm2, fcoeff, dx1, dy1, dz1;

  int n, count, i, rank, nprocs;
  unsigned long long tm1, tm2, tm3, tm4, total = 0;
  double t3, elapsed = 0.0, validation, final;
  double t1, t2;

  //MPI_Init( &argc, &argv );
  //MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  //MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  
  rank = 0;
  nprocs = 1;

  count = 327;

  if ( rank == 0 ) 
  {  
     printf( "count is set %d\n", count );
     printf( "Total MPI ranks %d\n", nprocs );
  } 

//#pragma omp parallel
//{
//  if ( (rank == 0) && (omp_get_thread_num() == 0) )
//  {
//     printf( "Number of OMP threads %d\n\n", omp_get_num_threads() );
     //printf( "      N         Time,us        Validation result\n" );
//  }   
//}

#ifdef TIMEBASE
  tm3 = timebase();
#endif

  final = 0.;

  for ( n = 400; n < N; n = n + 20 ) 
  {
      /* Initial data preparation */
      fcoeff = 0.23f;  
      fsrrmax2 = 0.5f; 
      mp_rsm2 = 0.03f;
      dx1 = 1.0f/(float)n;
      dy1 = 2.0f/(float)n;
      dz1 = 3.0f/(float)n;
      xx[0] = 0.f;
      yy[0] = 0.f;
      zz[0] = 0.f;
      mass[0] = 2.f;
      
      for ( i = 1; i < n; i++ )
      {
          xx[i] = xx[i-1] + dx1;
          yy[i] = yy[i-1] + dy1;
          zz[i] = zz[i-1] + dz1;
          mass[i] = (float)i * 0.01f + xx[i];
      }
    
      for ( i = 0; i < n; i++ )
      {
          vx1[i] = 0.f;
          vy1[i] = 0.f;
          vz1[i] = 0.f;
      }
    
      /* Data preparation done */
    
    
      /* Clean L1 cache */
      for ( i = 0; i < NC; i++ ) M1[i] = 4;
      for ( i = 0; i < NC; i++ ) M2[i] = M1[i];

#ifdef TIMEBASE
      tm1 = timebase();
#else
      t1 = mysecond();
#endif

      for ( i = 0; i < count; ++i)
      {
        Step10_orig( n, xx[i], yy[i], zz[i], fsrrmax2, mp_rsm2, xx, yy, zz, mass, &dx1, &dy1, &dz1 );
    
        vx1[i] = vx1[i] + dx1 * fcoeff;
        vy1[i] = vy1[i] + dy1 * fcoeff;
        vz1[i] = vz1[i] + dz1 * fcoeff;
      }

#ifdef TIMEBASE
      tm2 = timebase();
#else
      t2 = mysecond();
#endif

      validation = 0.;
      for ( i = 0; i < n; i++ )
      {
         validation = validation + ( vx1[i] + vy1[i] + vz1[i] );
      }

      final = final + validation;
    
#ifdef TIMEBASE      
      t3 = 1e6 * (double)(tm2 - tm1) / MHz; // time in us
#else
      t3 = (t2 - t1) * 1e6;
#endif 

      elapsed = elapsed + t3;
    
  }

#ifdef TIMEBASE
  tm4 = timebase();
  if ( rank == 0 )
  {
      printf( "\nKernel elapsed time, s: %18.8lf\n", elapsed*1e-6 );
      printf(   "Total  elapsed time, s: %18.8lf\n", (double)(tm4 - tm3) / MHz ); 
      printf(   "Result validation: %18.8lf\n", final );
      printf(   "Result expected  : 6636045675.12190628\n" );
  }
#else
  if ( rank == 0 )
  {
      printf( "\nKernel elapsed time, s: %18.8lf\n", elapsed*1e-6 );
  }    
#endif
  
  //MPI_Finalize();

  return 0;
}
