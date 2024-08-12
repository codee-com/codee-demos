#include <matrix.h>

// Creates a new dense matrix with the specified rows and columns
double **new_matrix(size_t rows, size_t cols) {
    if (rows < 1 || cols < 1)
        return NULL;

    // Allocate a dynamic array of doubles to store the matrix data linearized
    size_t matBytes = cols * rows * sizeof(double);
    double *memPtr = (double *)malloc(matBytes);
    if (!memPtr) {
        return NULL;
    }

    // Allocate an array of pointers to store the beginning of each row
    double **mat = (double **)calloc(rows, sizeof(double *));
    if (!mat) {
        free(memPtr);
        return NULL;
    }

    // Set the row pointers (eg. mat[2] points to the first double of row 3)
    for (size_t i = 0; i < rows; i++)
        mat[i] = memPtr + i * cols;

    return mat;
}

// Deletes the matrix and the resources allocated by it
void delete_matrix(double **mat) {
    if (mat) {
        free(mat[0]);
        free(mat);
    }
}

// Generates a random dense matrix
double **rand_matrix(double **mat, size_t rows, size_t cols) {
    if (!mat)
        return NULL;

    for (size_t row = 0; row < rows; row++)
        for (size_t col = 0; col < cols; col++)
            mat[row][col] = rand() % 10;
    return mat;
}

// Generates a checksum based on the matrix data
double checksum_matrix(double **mat, size_t rows, size_t cols) {
    if (!mat)
        return 0.0;

    double checkSum = 0.0;
    for (size_t row = 0; row < rows; row++)
        for (size_t col = 0; col < cols; col++)
            checkSum += mat[row][col];
    return checkSum;
}
