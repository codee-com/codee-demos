# NuCCOR

The Nuclear Coupled Cluster Oak Ridge (NuCCOR) algorithm calculates properties
of atomic nuclei from “first principles”.

The code of this repository (Microscopic Transport Code (MTC)) is just the
kernel of NuCCOR.

Both `mtc.F90` and `mtc_openmp.F90` are the same kernel algorithm. The idea is
to optimize `mtc_openmp.F90` with Codee and leave `mtc.F90` without changes, so
that the benchmark (`mtc_main.F90`) can measure the speedup achieved with the
help of Codee.

## How to run NuCCOR

The NuCCOR application setup is as follows. The usage of the executable binary:
`$ ./mtc.x nparticles basissize nblocks nonzero_fraction run_benchmark[yes/no] kernel_name`

### Suggested values for the parameters:
- nparticles: Number of particles in the system, typically values in the region from 16-132. A good testcase is 50.
- basissize: Basis size for tensors. Production sizes go up to 10k, but depending on total memory available, a good testcase is around 500.
- nblocks: Number of tensor symmetry blocks, typically values around 70-100. Depending on memory, a good testcase is 10.
- nonzero_fraction: Fraction of tensor with non-zero elements. A good testcase is any number <= 1.
- run_benchmark \[yes/no\]: if set to 'no', the program will parse and print the MTC parameters without running the benchmark.
- kernel_name: Choose which of the included kernels to run. If not specified, the 'contract' kernel is run.

### Example run of the NuCCOR executable binary used in the experiments:
```
$ ./mtc.x 30 70 10 0.1 yes contract_simple
 nh:           30
 np:           40
 nab:           16
 nc:            4
 nij:            9
 nk:            3
 Allocated cmap:  T
 Allocated kmap:  T
Memory usage for simple contract:     12.900 Gb
 Time spent in simple contraction reference:    60.642983010000535     
 Time spent in contraction Openmp version:    8.6296607419999418     
 Test simple contraction: OK
```
