// Include module header
#include "CRSMatrix.h"

// Include other headers
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


// Creates a new CRS sparse matrix (size = non-zero elements)
CRSMatrix *CRSMatrix_new(int rows, int cols, int size) {
    if (rows < 1 || cols < 1 || rows * cols < size)
        return 0;
    CRSMatrix *_this = (CRSMatrix *)malloc(sizeof(CRSMatrix));
    if (!_this)
        return 0;

    _this->rows = rows;
    _this->cols = cols;
    _this->size = size;
    _this->data = (double *)malloc(size * sizeof(double));
    _this->colRef = (int *)malloc(size * sizeof(int));
    _this->rowRef = (int *)malloc((rows + 1) * sizeof(int));
    if (_this->data && _this->colRef && _this->rowRef)
        return _this;

    // on memory allocation error
    if (_this->data)
        free(_this->data);
    if (_this->colRef)
        free(_this->colRef);
    if (_this->rowRef)
        free(_this->rowRef);
    if (_this)
        free(_this);
    return 0;
}


// Creates a new CRS sparse matrix from a dense matrix
CRSMatrix *CRSMatrix_from(const Matrix2D *mat) {
    if (!mat)
        return 0;
    size_t nonZero = 0;
    const double EPS = 1e-9;
    for (int row = 0; row < mat->rows; row++)
        for (int col = 0; col < mat->cols; col++)
            if (fabs(mat->data[row][col]) > EPS)
                nonZero++;

    CRSMatrix *CRSMat = CRSMatrix_new(mat->rows, mat->cols, nonZero);
    if (!CRSMat)
        return 0;

    long long sparsePos = 0;
    for (int row = 0; row < mat->rows; row++) {
        CRSMat->rowRef[row] = sparsePos;
        for (int col = 0; col < mat->cols; col++) {
            double denseValue = mat->data[row][col];
            if (fabs(denseValue) < EPS)
                continue;

            CRSMat->data[sparsePos] = denseValue;
            CRSMat->colRef[sparsePos] = col;
            sparsePos++;
        }
    }

    CRSMat->rowRef[mat->rows] = sparsePos;

    return CRSMat;
}


// Deletes the CRS sparse matrix and the resources allocated by it
void CRSMatrix_delete(CRSMatrix *mat) {
    if (!mat)
        return;
    if (mat->data)
        free(mat->data);
    if (mat->colRef)
        free(mat->colRef);
    if (mat->rowRef)
        free(mat->rowRef);
    free(mat);
}


// Obtains CRS matrix value at the specified location
double CRSMatrix_getVal(const CRSMatrix *mat, int row, int col) {
    assert(mat && row >= 0 && col >= 0);
    assert(row < mat->rows && col < mat->cols);
    long long initRow = mat->rowRef[row];
    long long endRow = mat->rowRef[row + 1];

    for (long long sparsePos = initRow; sparsePos < endRow; sparsePos++) {
        int sparseCol = mat->colRef[sparsePos];
        if (sparseCol > col)
            break;
        if (col == sparseCol)
            return mat->data[sparsePos];
    }

    return 0.0;
}


// Prints the selected sparse matrix in dense form
void CRSMatrix_print(const CRSMatrix *mat) {
    if (!mat) {
        printf("CRSMatrix: <NULL>\n");
        return;
    }

    for (int row = 0; row < mat->rows; row++) {
        printf("%c[ ", row == 0 ? '[' : ' ');
        for (int col = 0; col < mat->cols; col++)
            printf("%5.1f ", CRSMatrix_getVal(mat, row, col));
        printf("]%c\n", row == mat->cols - 1 ? ']' : ' ');
    }
}


// Prints the selected sparse matrix in CRS form
void CRSMatrix_debug(const CRSMatrix *mat) {
    if (!mat) {
        printf("CRSMatrix: <NULL>\n");
        return;
    }
    printf("CRSMatrix: %i x %i\n", mat->rows, mat->cols);

    long long size = mat->size;
    printf("  data  [%2Li] = {", size);
    for (long long i = 0; i < size; i++)
        printf("%5.1f%c", mat->data[i], i == size - 1 ? ' ' : ',');
    printf("}\n");

    printf("  colRef[%2Li] = {", size);
    for (long long i = 0; i < size; i++)
        printf("%5i%c", mat->colRef[i], i == size - 1 ? ' ' : ',');
    printf("}\n");

    long long sparseRows = mat->rows + 1;
    printf("  rowRef[%2Li] = {", sparseRows);
    for (long long i = 0; i < sparseRows; i++)
        printf("%5i%c", mat->rowRef[i], i == sparseRows - 1 ? ' ' : ',');
    printf("}\n");
}
