#pragma once
#ifndef _VECTOR_H_
#define _VECTOR_H_

#ifndef vType
	#define vType double
#endif

// Vector structure
typedef struct Vector {
	long long size;
	vType* data;
} Vector;

// Creates a new vector with the specified size
Vector* Vector_new(long long size);

// Deletes the vector and the resources allocated by it
void Vector_delete(Vector* vec);

// Generates a random vector
Vector* Vector_rand(Vector* vec);

// Generates a sequential vector
Vector* Vector_seq(Vector* vec);

// Generates an empty vector
Vector* Vector_zero(Vector* vec);

// Generates a checksum based on the vector data
vType Vector_checksum(Vector* vec);

// Returns the matrix value at the specified location
vType Vector_getVal(const Vector* vec, long long pos);

// Prints the selected vector
void Vector_print(const Vector* mat);

// Get total number of elements
#define Vector_getSize(vecPtr) vecPtr->size

// Get raw pointer to vector data
#define Vector_getData(vecPtr) vecPtr->data

#endif
