#!/bin/bash -x

export BG_SHAREDMEMSIZE=32
export L1P_POLICY=std

export PROG=HACCmk
export NODES=1

export RANKS_PER_NODE=1
export OMP_NUM_THREADS=8
export BG_THREADLAYOUT=1   # 1 - default next core first; 2 - my core first


export NPROCS=$((NODES*RANKS_PER_NODE)) 
export OUTPUT=HACCmk_${NPROCS}_${OMP_NUM_THREADS}
export VARS="FAST_WAKEUP=TRUE:BG_SHAREDMEMSIZE=${BG_SHAREDMEMSIZE}:OMP_NUM_THREADS=${OMP_NUM_THREADS}:L1P_POLICY=${L1P_POLICY}:BG_THREADLAYOUT=${BG_THREADLAYOUT}:XLSMPOPTS=stack=8000000"
rm -f core.* ${OUTPUT}.cobaltlog ${OUTPUT}.error ${OUTPUT}.output
qsub -A Performance -n ${NODES} --proccount ${RANKS_PER_NODE} --mode c${RANKS_PER_NODE} -t 0:10:00 -O $OUTPUT --env ${VARS} $PROG
