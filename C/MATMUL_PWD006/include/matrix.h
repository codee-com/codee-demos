#include <stdlib.h>

/*   Matrices stored as a linearized array plus a row pointer array
                m
                |
                v
              { a,       b,       c }
a -> 1 0 0      |        |        |
b -> 0 2 0      v        v        v
b -> 0 0 3 => { 1, 0, 0, 0, 2, 0, 0, 0, 3 } => m[1][1] == 2
                                                   ^   == a[1*3+1] == *(m[0]+3+1)
                                                   |   == b[1]     == *(m[1]+1)
         The row pointer array enables this syntax +
*/

double **new_matrix(size_t rows, size_t cols);
void delete_matrix(double **mat);
double **rand_matrix(double **mat, size_t rows, size_t cols);
double checksum_matrix(double **mat, size_t rows, size_t cols);
