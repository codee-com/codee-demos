#pragma once
#ifndef _MATRIX2D_H_
#define _MATRIX2D_H_

// Matrix structure
typedef struct Matrix2D {
	int rows;
	int cols;
	long long size;
	double** data;
} Matrix2D;

// Creates a new dense matrix with the specified rows and columns
Matrix2D* Matrix2D_new(int rows, int cols);

// Deletes the matrix and the resources allocated by it
void Matrix2D_delete(Matrix2D* mat);

// Generates a random matrix
Matrix2D* Matrix2D_rand(Matrix2D* mat);

// Generates a sequential matrix
Matrix2D* Matrix2D_seq(Matrix2D* mat);

// Generates an empty matrix
Matrix2D* Matrix2D_zero(Matrix2D* mat);

// Generates an identity matrix
Matrix2D* Matrix2D_identity(Matrix2D* mat);

// Generates a random matrix with sparse data
Matrix2D* Matrix2D_randSparse(Matrix2D* mat, float sparsity);

// Generates a checksum based on the matrix data
double Matrix2D_checksum(Matrix2D* mat);

// Returns the matrix value at the specified row and column
double Matrix2D_getVal(const Matrix2D* mat, int row, int col);

// Prints the selected matrix
void Matrix2D_print(const Matrix2D* mat);

// Get number of rows
#define Matrix2D_getRows(matPtr) matPtr->rows

// Get number of columns
#define Matrix2D_getCols(matPtr) matPtr->cols

// Get total number of elements
#define Matrix2D_getSize(matPtr) matPtr->size

// Get raw pointer to matrix data (as list of pointers to rows)
#define Matrix2D_getData(matPtr) matPtr->data

// Get raw pointer to matrix data (as consecutive memory)
#define Matrix2D_get1D(matPtr) matPtr->data[0]

#endif
