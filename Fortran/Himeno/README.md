# Himeno Benchmark

Fortran 90 version source https://i.riken.jp/en/supercom/documents/himenobmt/

# HowTo

## Compile with CMAKE
```
$ cmake . \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -G "Ninja" \
  -B build
```
`$ ninja -C build/ run`

## Compile with Makefile
```
$ make build
```
`$ make run`
