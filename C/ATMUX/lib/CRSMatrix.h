#pragma once
#ifndef _CRSMATRIX_H_
#define _CRSMATRIX_H_

#include "Matrix2D.h"

typedef struct CRSMatrix {
    int rows;
    int cols;
    long long size;
    double *data;
    int *colRef;
    int *rowRef;
} CRSMatrix;

// Creates a new CRS sparse matrix (size = non-zero elements)
CRSMatrix *CRSMatrix_new(int rows, int cols, int size);

// Creates a new CRS sparse matrix from a dense matrix
CRSMatrix *CRSMatrix_from(const Matrix2D *mat);

// Deletes the CRS sparse matrix and the resources allocated by it
void CRSMatrix_delete(CRSMatrix *mat);

// Prints the selected sparse matrix in dense form
void CRSMatrix_print(const CRSMatrix *mat);

// Prints the selected sparse matrix in CRS form
void CRSMatrix_debug(const CRSMatrix *mat);

// Obtains CRS matrix value at the specified location
double CRSMatrix_getVal(const CRSMatrix *mat, int row, int col);

// Get number of rows
#define CRSMatrix_getRows(matPtr) matPtr->rows

// Get number of columns
#define CRSMatrix_getCols(matPtr) matPtr->cols

// Get total number of elements
#define CRSMatrix_getSize(matPtr) matPtr->size

// Get raw pointer to matrix data
#define CRSMatrix_getData(matPtr) matPtr->data

// Get raw pointer to CRS column references
#define CRSMatrix_colRef(matPtr) matPtr->colRef

// Get raw pointer to CRS row references
#define CRSMatrix_rowRef(matPtr) matPtr->rowRef

#endif
