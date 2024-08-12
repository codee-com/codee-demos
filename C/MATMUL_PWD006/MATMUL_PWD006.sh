#!/bin/bash

function printRunCommand(){
    ## Print the command
    printf "\n$ $@\n"
    ## Run the command
    $@
}

module load PrgEnv-nvidia

rm -f matmul_pwd006 *.o

echo ""
echo ""
echo "Matmul"
printRunCommand "nvc -fast -mp -target=gpu -Minfo=mp -I include matrix.c clock.c main_pwd006.c -o matmul_pwd006"
printRunCommand "./matmul_pwd006 3000"
