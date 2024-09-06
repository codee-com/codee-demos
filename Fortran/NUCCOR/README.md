# NuCCOR

The Nuclear Coupled Cluster Oak Ridge (NuCCOR) algorithm calculates properties
of atomic nuclei from “first principles”.

The code of this repository (Microscopic Transport Code (MTC)) is just the
kernel of NuCCOR.

Both `mtc.F90` and `mtc_openmp.F90` are the same kernel algorithm. The idea is
to optimize `mtc_openmp.F90` with Codee and leave `mtc.F90` without changes, so
that the benchmark (`mtc_main.F90`) can measure the speedup achieved with the
help of Codee.

## How to run NUCCOR

Example of invoking the benchmark binary:
`./mtc.x 8 100 8 0.5`
