// Include module header
#include "Matrix2D.h"

// Include ohter headers
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// Creates a new dense matrix with the specified rows and columns
Matrix2D* Matrix2D_new(int rows, int cols) {
	if(rows < 1 || cols < 1) return 0;
	Matrix2D* _this = (Matrix2D*)malloc(sizeof(Matrix2D));
	if(!_this) return 0;
	
	_this->rows = rows;
	_this->cols = cols;
	_this->size = rows * cols;
	// Last row pointer is always null
	_this->data = (double**)calloc(rows + 1, sizeof(double));
	if(_this->data == 0) return 0;
	
	size_t matBytes = cols * rows * sizeof(double);
	double *memPtr = (double*)malloc(matBytes);
	if(memPtr) {
		for(int i = 0; i < rows; i++)
			_this->data[i] = memPtr + i * cols;
		return _this;
	}

	// on memory allocation error
	if(_this->data) free(_this->data);
	if(_this) free(_this);
	return 0;
}

// Deletes the matrix and the resources allocated by it
void Matrix2D_delete(Matrix2D* mat) {
	if(mat && mat->data) {
		if(mat->data[0])
			free(mat->data[0]);
		free(mat->data);
	}
	if(mat) free(mat);
}

// Generates a random dense matrix
Matrix2D* Matrix2D_rand(Matrix2D* mat) {
	if(!mat) return 0;
	for(int row = 0; row < mat->rows; row++)
		for(int col = 0; col < mat->cols; col++)
			mat->data[row][col] = rand() % 10;
	return mat;
}

// Generates a sequential matrix
Matrix2D* Matrix2D_seq(Matrix2D* mat) {
	if(!mat) return 0;
	int cont = 1;
	for(int row = 0; row < mat->rows; row++)
		for(int col = 0; col < mat->cols; col++)
			mat->data[row][col] = cont++;
	return mat;
}

// Generates an identity matrix
Matrix2D* Matrix2D_identity(Matrix2D* mat) {
	if(!mat) return 0;
	for(int row = 0; row < mat->rows; row++)
		for(int col = 0; col < mat->cols; col++)
			mat->data[row][col] = (row == col) ? 1 : 0;
	return mat;
}

// Generates an empty matrix
Matrix2D* Matrix2D_zero(Matrix2D* mat) {
	if(!mat) return 0;
	size_t nBytes = sizeof(double) * mat->size;
	memset(mat->data[0], 0, nBytes);
	return mat;
}

// Generates a random sparse matrix
Matrix2D* Matrix2D_randSparse(Matrix2D* mat, float sparsity) {
	if(!mat) return 0;
	float prob = 1.0f - sparsity;
	for(int row = 0; row < mat->rows; row++) {
		for(int col = 0; col < mat->cols; col++) {
			int value = rand() % 99 + 1;
			double mask = (double)rand() / RAND_MAX < prob;
			mat->data[row][col] = value * mask;
		}
	}
	return mat;
}

// Generates a checksum based on the matrix data
double Matrix2D_checksum(Matrix2D* mat) {
	if(!mat) return 0;
	
	double checkSum = 0.0;
	for(int row = 0; row < mat->rows; row++)
		for(int col = 0; col < mat->cols; col++)
			checkSum += mat->data[row][col];
	return checkSum;
}

// Returns the matrix value at the specified row and column
double Matrix2D_getVal(const Matrix2D* mat, int row, int col) {
	assert(mat && row >= 0 && col >= 0);
	assert(row < mat->rows && col < mat->cols);
	return mat->data[row][col];
}

// Prints the selected dense matrix
void Matrix2D_print(const Matrix2D* mat) {
	if(!mat) { printf("Matrix2D: <NULL>\n"); return; }

	for(int row = 0; row < mat->rows; row++) {
		printf("%c[ ", row == 0 ? '[' : ' ');
		for(int col = 0; col < mat->cols; col++)
			printf("%.3e ", mat->data[row][col]);
		printf("]%c\n", row == mat->cols - 1 ? ']' : ' ');
	}
}
