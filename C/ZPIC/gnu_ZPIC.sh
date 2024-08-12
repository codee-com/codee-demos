#!/bin/bash

module use /global/cfs/cdirs/m1759/yunhe/Modules/perlmutter/modulefiles
module load gcc/12.1.0

echo ""
echo "em2d-0-serial"
echo ""
make -C em2d-0-serial/ CC=gcc CFLAGS="-Ofast" clean zpic
em2d-0-serial/zpic

echo ""
echo "em2d-1-outlining"
echo ""
make -C em2d-1-outlining/ CC=gcc CFLAGS="-Ofast" clean zpic
em2d-1-outlining/zpic

echo ""
echo "em2d-2-inlining"
echo ""
make -C em2d-2-inlining/ CC=gcc CFLAGS="-Ofast" clean zpic
em2d-2-inlining/zpic

echo ""
echo "em2d-3-aos"
echo ""
make -C em2d-3-aos/ CC=gcc CFLAGS="-Ofast" clean zpic
em2d-3-aos/zpic

echo ""
echo "em2d-4-loopfission"
echo ""
make -C em2d-4-loopfission/ CC=gcc CFLAGS="-Ofast" clean zpic
em2d-4-loopfission/zpic

echo ""
echo "em2d-5a-omp-atomic-cpu"
echo ""
make -C em2d-5a-omp-atomic-cpu/ CC=gcc CFLAGS="-Ofast -fopenmp" clean zpic
em2d-5a-omp-atomic-cpu/zpic

echo ""
echo "em2d-6a-omp-atomic-gpu"
echo ""
make -C em2d-6a-omp-atomic-gpu/ CC=gcc CFLAGS="-Ofast -fopenmp -foffload=nvptx-none=\"-Ofast -misa=sm_80\"" clean zpic
em2d-6a-omp-atomic-gpu/zpic

echo ""
echo "em2d-7a-acc-atomic-gpu"
echo ""
make -C em2d-7a-acc-atomic-gpu/ CC=gcc CFLAGS="-Ofast -openacc -foffload=nvptx-none=\"-Ofast -misa=sm_80\"" clean zpic
em2d-7a-acc-atomic-gpu/zpic
