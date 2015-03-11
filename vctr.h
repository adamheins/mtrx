#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#include "sclr.h"


#ifndef VCTR_H_
#define VCTR_H_


/*===========================================================================*
 *                               STRUCTURES                                  *
 *===========================================================================*/

/*****************************************************************************
 * Defines a 1-dimensional vector.                                           *
 *                                                                           *
 * Fields:                                                                   *
 *     values - The array of double precision values stored in the vector.   *
 *     length - The length of the vector.                                    *
 *****************************************************************************/
typedef struct vector {
  scalar_t *values;
  size_t length;
} vector_t;


/*===========================================================================*
 *                              FUNCTIONS                                    *
 *===========================================================================*/

/*--------------------------- Initialization --------------------------------*/

/*****************************************************************************
 * Creates a vector from the string of values.                               *
 *                                                                           *
 * Fields:                                                                   *
 *     str - The string to use to create the vector.                         *
 *                                                                           *
 * Returns: A pointer to the new vector.                                     *
 *****************************************************************************/
vector_t *vctr_init(const char *str);

/*****************************************************************************
 * Creates a vector of specified length with uninitialized values.           *
 *                                                                           *
 * Fields:                                                                   *
 *     length - The length of the vector being created.                      *
 *                                                                           *
 * Returns: A pointer to the empty vector.                                   *
 *****************************************************************************/
vector_t *vctr_empty(size_t length);


/*****************************************************************************
 * Creates a vector of specified length with all values set to zero.         *
 *                                                                           *
 * Fields:                                                                   *
 *     length - The length of the vector being created.                      *
 *                                                                           *
 * Returns: A pointer to the created vector.                                 *
 *****************************************************************************/
vector_t *vctr_zeros(size_t length);


/*****************************************************************************
 * Creates a vector of specified length with all values set to one.          *
 *                                                                           *
 * Fields:                                                                   *
 *     length - The length of the vector being created.                      *
 *                                                                           *
 * Returns: A pointer to the created vector.                                 *
 *****************************************************************************/
vector_t *vctr_ones(size_t length);


/*****************************************************************************
 * Creates a vector of specified length containing random intergers between  *
 * 0 and the given maximum value.                                            *
 *                                                                           *
 * Fields:                                                                   *
 *     length - The length of the vector being created.                      *
 *     max - The maximum possible value in the vector, exclusive.            *
 *                                                                           *
 * Returns: A pointer to the created vector.                                 *
 *****************************************************************************/
vector_t *vctr_rnd(size_t length, uint32_t max);


/*****************************************************************************
 * Creates a new vector that is a copy of the existing one.                  *
 *                                                                           *
 * Fields:                                                                   *
 *     vector - The vector to be copied.                                     *
 *                                                                           *
 * Returns: A pointer to the new copied vector.                              *
 *****************************************************************************/
vector_t *vctr_copy(vector_t *vector);


/*---------------------------- Deallocation ---------------------------------*/

/*****************************************************************************
 * Frees all of the memory associated with this vector and sets it equal to  *
 * NULL.                                                                     *
 *                                                                           *
 * Fields:                                                                   *
 *     vector - The vector to be freed.                                      *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void vctr_destroy(vector_t *vector);


/*-------------------------------- Display ----------------------------------*/

/*****************************************************************************
 * Prints the values contained in the vector to standard output.             *
 *                                                                           *
 * Fields:                                                                   *
 *     vector - The vector to print.                                         *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void vctr_print(vector_t * vector);


/*------------------------------ Comparisons --------------------------------*/

/*****************************************************************************
 * Compares the length of two vectors, returning true if the length of the   *
 * vectors are equals, false otherwise.                                      *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first vector being compared.                                  *
 *     B - The second vector being compared.                                 *
 *                                                                           *
 * Returns: True if the vectors have equal length, false otherwise.          *
 *****************************************************************************/
bool vctr_eq_len(vector_t *A, vector_t *B);


/*****************************************************************************
 * Compares two vectors. Returns true if the values contained by each vector *
 * is the same, false otherwise.                                             *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first vector being compared.                                  *
 *     B - The second vector being compared.                                 *
 *                                                                           *
 * Returns: True if the vectors are equal, false otherwise.                  *
 *****************************************************************************/
bool vctr_eq(vector_t *A, vector_t *B);


/*------------------------------ Max and Min --------------------------------*/

/*****************************************************************************
 * Returns the maximum value in the vector.                                  *
 *                                                                           *
 * Fields:                                                                   *
 *     vector - The vector in which to find the maximum value.               *
 *                                                                           *
 * Returns: The maximum value in the vector.                                 *
 *****************************************************************************/
scalar_t vctr_max(vector_t *vector);


/*****************************************************************************
 * Returns the minimum value in the vector.                                  *
 *                                                                           *
 * Fields:                                                                   *
 *     vector - The vector in which to find the minimum value.               *
 *                                                                           *
 * Returns: The minimum value in the vector.                                 *
 *****************************************************************************/
scalar_t vctr_min(vector_t *vector);


/*------------------------------- Operations --------------------------------*/

/*****************************************************************************
 * Calculates the dot product of two vectors.                                *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first vector being dotted.                                    *
 *     B - The second vector being dotted.                                   *
 *                                                                           *
 * Returns: The dot (scalar) product of the two vectors.                     *
 *****************************************************************************/
scalar_t vctr_dot_prod(vector_t *A, vector_t *B);


/*****************************************************************************
 * Calculates the cross product of two vectors. The vectors must be of       *
 * length 3.                                                                 *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first vector being crossed.                                   *
 *     B - The second vector being crossed.                                  *
 *                                                                           *
 * Returns: The cross (vector) product of the two vectors.                   *
 *****************************************************************************/
vector_t *vctr_cross_prod(vector_t *A, vector_t *B);


/*****************************************************************************
 * Calculates the magnitude of a vector.                                     *
 *                                                                           *
 * Fields:                                                                   *
 *     vector - The vector for which the magnitude is being calculated.      *
 *                                                                           *
 * Returns: The magnitude of the vector.                                     *
 *****************************************************************************/
scalar_t vctr_mag(vector_t *vector);

#endif
