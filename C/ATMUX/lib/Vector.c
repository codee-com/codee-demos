// Include module header
#include "Vector.h"

// Include other headers
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Creates a new vector with the specified size
Vector *Vector_new(long long size) {
    if (size < 1)
        return 0;
    Vector *_this = (Vector *)malloc(sizeof(Vector));
    if (!_this)
        return 0;

    _this->size = size;
    size_t nBytes = size * sizeof(double);
    _this->data = (double *)malloc(nBytes);
    if (_this->data)
        return _this;

    // on memory allocation error
    if (_this)
        free(_this);
    return 0;
}

// Deletes the vector and the resources allocated by it
void Vector_delete(Vector *vec) {
    if (vec->data)
        free(vec->data);
    if (vec)
        free(vec);
}

// Generates a random vector
Vector *Vector_rand(Vector *vec) {
    if (!vec)
        return 0;
    for (size_t i = 0; i < vec->size; i++)
        vec->data[i] = rand() % 100;
    return vec;
}

// Generates a sequential vector
Vector *Vector_seq(Vector *vec) {
    if (!vec)
        return 0;
    for (size_t i = 0; i < vec->size; i++)
        vec->data[i] = (double)i;
    return vec;
}

// Generates an empty vector
Vector *Vector_zero(Vector *vec) {
    if (!vec)
        return 0;
    size_t nBytes = sizeof(double) * vec->size;
    memset(vec->data, 0, nBytes);
    return vec;
}

double Vector_checksum(Vector *vec) {
    if (!vec)
        return 0;

    double checkSum = 0.0;
    for (int i = 0; i < vec->size; i++)
        checkSum += vec->data[i];
    return checkSum;
}

// Returns the vector value at the specified location
double Vector_getVal(const Vector *vec, long long pos) {
    assert(vec && pos >= 0);
    assert(pos < vec->size);
    return vec->data[pos];
}

// Prints the selected vector
void Vector_print(const Vector *vec) {
    if (!vec) {
        printf("Vector: <NULL>\n");
        return;
    }
    printf("{ ");
    for (long long pos = 0; pos < vec->size; pos++) {
        printf("%5.1f", vec->data[pos]);
        printf("%c", pos == vec->size - 1 ? ' ' : ',');
    }
    printf("}\n");
}
