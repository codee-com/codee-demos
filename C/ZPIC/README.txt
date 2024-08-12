## For ZPIC ew provide multiple versions of the source code, so that starting on the original code
## you can produce CPU and GPU code using OpenMP and OpenACC directives.
## At each step the corresponding ZPIC version can be build and run on Perlmutter.
## We provide two scripts that build and run all the versions:
##    nvidia_ZPIC.sh: run the projects using the nvidia compiler.
##    gnu_ZPIC.sh   : run the project using the gcc 12.1 compiler.

ZPIC EM2D versions
------------------

- em2d-0-serial
- em2d-1-outlining
- em2d-2-inlining
- em2d-3-aos
- em2d-4-loopfission
- em2d-5a-omp-atomic-cpu
- em2d-6a-omp-atomic-gpu
- em2d-7a-acc-atomic-gpu
