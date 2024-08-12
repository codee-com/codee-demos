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
$ FC=gfortran make build
```
`$ FC=gfortran make run`
