#!/bin/bash

module load PrgEnv-nvidia

echo ""
echo "em2d-0-serial"
echo ""
make -C em2d-0-serial/ CC=nvc CFLAGS="-fast" clean zpic
em2d-0-serial/zpic

echo ""
echo "em2d-1-outlining"
echo ""
make -C em2d-1-outlining/ CC=nvc CFLAGS="-fast" clean zpic
em2d-1-outlining/zpic

echo ""
echo "em2d-2-inlining"
echo ""
make -C em2d-2-inlining/ CC=nvc CFLAGS="-fast" clean zpic
em2d-2-inlining/zpic

echo ""
echo "em2d-3-aos"
echo ""
make -C em2d-3-aos/ CC=nvc CFLAGS="-fast" clean zpic
em2d-3-aos/zpic

echo ""
echo "em2d-4-loopfission"
echo ""
make -C em2d-4-loopfission/ CC=nvc CFLAGS="-fast" clean zpic
em2d-4-loopfission/zpic

echo ""
echo "em2d-5a-omp-atomic-cpu"
echo ""
make -C em2d-5a-omp-atomic-cpu/ CC=nvc CFLAGS="-fast -mp -target=multicore -Minfo=mp" clean zpic
em2d-5a-omp-atomic-cpu/zpic

echo ""
echo "em2d-6a-omp-atomic-gpu"
echo ""
make -C em2d-6a-omp-atomic-gpu/ CC=nvc CFLAGS="-fast -mp -target=gpu -Minfo=mp" clean zpic
em2d-6a-omp-atomic-gpu/zpic

echo ""
echo "em2d-7a-acc-atomic-gpu"
echo ""
make -C em2d-7a-acc-atomic-gpu/ CC=nvc CFLAGS="-fast -acc -target=gpu -Minfo=acc" clean zpic
em2d-7a-acc-atomic-gpu/zpic
