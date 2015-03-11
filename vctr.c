#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

#include "sclr.h"
#include "mtrx_util.h"
#include "vctr.h"


/*===========================================================================*
 *                           VECTOR FUNCTIONS                                *
 *===========================================================================*/

/*--------------------------- Initialization --------------------------------*/

vector_t *vctr_init(const char *str) {

  // Count the number of elements the array will have.
  size_t count = 0;
  if (str[0] != ' ')
    ++count;

  for (size_t i = 1; str[i] != '\0'; ++i) {
    if (str[i - 1] == ' ' && str[i] != ' ')
      ++count;
  }

  // Create the array.
  vector_t *vector = vctr_empty(count);

  // Parse each scalar value out of the array.
  count = 0;
  size_t start = 0;
  size_t i = 1;
  for (; str[i] != '\0'; ++i) {

    // 'Rising edge'
    if (str[i - 1] == ' ' && str[i] != ' ')
      start = i;

    // 'Falling edge'
    if (str[i - 1] != ' ' && str[i] == ' ')
      vector->values[count++] = atos(str, start, i);
  }
  if (count < vector->length)
    vector->values[count] = atos(str, start, i);
  return vector;
}


vector_t *vctr_empty(size_t length) {
	vector_t *vector = (vector_t *)malloc(sizeof(vector_t));

	vector->values = (scalar_t *)malloc(length * sizeof(scalar_t));
	vector->length = length;

	return vector;
}


vector_t *vctr_zeros(size_t length) {
	vector_t *vector = (vector_t *)malloc(sizeof(vector_t));

	vector->values = (scalar_t *)calloc(length, sizeof(scalar_t));
	vector->length = length;

	return vector;
}


vector_t *vctr_ones(size_t length) {
	vector_t *vector = (vector_t *)malloc(sizeof(vector_t));

	vector->values = (scalar_t *)malloc(length * sizeof(scalar_t));
	vector->length = length;

	for (size_t i = 0; i < length; ++i)
		vector->values[i] = 1;

	return vector;
}


vector_t *vctr_rnd(size_t length, uint32_t max) {
  vector_t *vector = vctr_empty(length);

  for (size_t i = 0; i < length; ++i)
    vector->values[i] = rand() % max;

  return vector;
}


vector_t *vctr_copy(vector_t *vector) {
	vector_t *copy = vctr_empty(vector->length);
	for (size_t i = 0; i < vector->length; ++i)
		copy->values[i] = vector->values[i];
	return copy;
}


/*---------------------------- Deallocation ---------------------------------*/

void vctr_destroy(vector_t *vector) {
	free(vector->values);
	free(vector);
	vector = NULL;
}


/*-------------------------------- Display ----------------------------------*/

void vctr_print(vector_t *vector) {
	for (size_t i = 0; i < vector->length; ++i)
		printf("%f\n", vector->values[i]);
}


/*------------------------------ Comparisons --------------------------------*/

bool vctr_eq_len(vector_t *A, vector_t *B) {
	return (A->length == B->length);
}


bool vctr_eq(vector_t *A, vector_t *B) {
  if (!vctr_eq_len(A, B))
    return false;
  for (size_t i = 0; i < A->length; ++i) {
    if (A->values[i] != B->values[i])
      return false;
  }
  return true;
}


/*------------------------------ Max and Min --------------------------------*/

scalar_t vctr_max(vector_t *vector) {
  scalar_t max = vector->values[0];
  for (size_t i = 1; i < vector->length; ++i) {
    if (vector->values[i] > max)
      max = vector->values[i];
  }
  return max;
}


scalar_t vctr_min(vector_t *vector) {
  scalar_t min = vector->values[0];
  for (size_t i = 1; i < vector->length; ++i) {
    if (vector->values[i] < min)
      min = vector->values[i];
  }
  return min;
}


/*------------------------------- Operations --------------------------------*/

scalar_t vctr_dot_prod(vector_t *A, vector_t *B) {

	// Check that the vectors are of equal length.
	if (!vctr_eq_len(A, B))
		return 0;

	scalar_t prod = 0;
	for (size_t i = 0; i < A->length; ++i)
		prod += A->values[i] * B->values[i];

	return prod;
}


vector_t *vctr_cross_prod(vector_t *A, vector_t *B) {

	// Vectors must be of length 3 to perform the cross product.
	if (A->length != 3 || B->length != 3)
		return NULL;

	vector_t *C = vctr_empty(3);

	C->values[0] = A->values[1] * B->values[2] - A->values[2] * B->values[1];
	C->values[1] = A->values[2] * B->values[0] - A->values[0] * B->values[2];
	C->values[2] = A->values[0] * B->values[1] - A->values[1] * B->values[0];

	return C;
}


scalar_t vctr_mag(vector_t *vector) {

	scalar_t mag = 0;

	for (size_t i = 0; i < vector->length; ++i)
		mag += vector->values[i] * vector->values[i];

	return sqrt(mag);
}

