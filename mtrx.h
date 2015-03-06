/*****************************************************************************
 * A simple linear algebra library containing functions for vector and       *
 * matrix operations.                                                        *
 *                                                                           *
 * Author: Adam Heins                                                        *
 *                                                                           *
 * Last Updated: 2014-11-29                                                  *
 *****************************************************************************/

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#ifndef MTRX_H_
#define MTRX_H_


/*===========================================================================*
 *                               STRUCTURES                                  *
 *===========================================================================*/

/*****************************************************************************
 * Defines a scalar type.                                                    *
 *****************************************************************************/
typedef double scalar_t;


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


/*****************************************************************************
 * Defines a 2-dimensional matrix.                                           *
 *                                                                           *
 * Fields:                                                                   *
 *     values  - The 2D array of double precision values stored in the       *
 *               matrix.                                                     *
 *	   rows    - The number of rows this matrix has.                         *
 *     columns - The number of columns this matrix has.                      *
 *****************************************************************************/
typedef struct matrix {
	scalar_t **values;
	size_t rows;
	size_t columns;
} matrix_t;


/*****************************************************************************
 * Defines an indexing array. This structure contains a series of integer    *
 * indexes for accessing elements in vectors and matrices.                   *
 *                                                                           *
 * Fields:                                                                   *
 *     values - The indexes contained by the indexer.                        *
 *     length - The length of the indexer.                                   *
 *****************************************************************************/
typedef struct indexer {
	size_t *values;
	size_t length;
} indexer_t;


/*===========================================================================*
 *                           VECTOR FUNCTIONS                                *
 *===========================================================================*/

/*--------------------------- Initialization --------------------------------*/

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


/*===========================================================================*
 *                           MATRIX FUNCTIONS                                *
 *===========================================================================*/

/*--------------------------- Initialization --------------------------------*/

/*****************************************************************************
 * Creates a matrix with the specified number of rows and columns.           *
 *                                                                           *
 * Fields:                                                                   *
 *     rows - The number of rows in the matrix.                              *
 *     cols - The number of columns in the matrix.                           *
 *                                                                           *
 * Returns: A pointer to the created matrix.                                 *
 *****************************************************************************/
matrix_t *mtrx_empty(size_t rows, size_t cols);


/*****************************************************************************
 * Creates a matrix with the specified number of rows and columns and all    *
 * values initialized to zero.                                               *
 *                                                                           *
 * Fields:                                                                   *
 *     rows - The number of rows in the matrix.                              *
 *     cols - The number of columns in the matrix.                           *
 *                                                                           *
 * Returns: A pointer to the created matrix.                                 *
 *****************************************************************************/
matrix_t *mtrx_zeros(size_t rows, size_t cols);


/*****************************************************************************
 * Creates a matrix with the specified number of rows and columns and all    *
 * values initialized to one.                                                *
 *                                                                           *
 * Fields:                                                                   *
 *     rows - The number of rows in the matrix.                              *
 *     cols - The number of columns in the matrix.                           *
 *                                                                           *
 * Returns: A pointer to the created matrix.                                 *
 *****************************************************************************/
matrix_t *mtrx_ones(size_t rows, size_t cols);


/*****************************************************************************
 * Creates an identity matrix of size n x n.                                 *
 *                                                                           *
 * Fields:                                                                   *
 *     n - The number of rows and columns in the matrix.                     *
 *                                                                           *
 * Returns: A pointer to the created identity matrix.                        *
 *****************************************************************************/
matrix_t *mtrx_id(size_t n);


/*****************************************************************************
 * Creates a diagonal matrix using a vector to specify the diagonal entries. *
 *                                                                           *
 * Fields:                                                                   *
 *     vector - The vector to place on the diagonal of the matrix.           *
 *                                                                           *
 * Returns: A pointer to the created diagonal matrix.                        *
 *****************************************************************************/
matrix_t *mtrx_diag(vector_t *vector);


/*****************************************************************************
 * Creates a new a matrix filled with random positive integers, up to the    *
 * given maximum value.                                                      *
 *                                                                           *
 * Fields:                                                                   *
 *     rows - The number of rows in the matrix.                              *
 *     cols - The number of columns in the matrix.                           *
 *     max - The maximum possible value in the matrix.                       *
 *                                                                           *
 * Returns: A pointer to the new matrix copy.                                *
 *****************************************************************************/
matrix_t *mtrx_rnd(size_t rows, size_t cols, uint32_t max);


/*****************************************************************************
 * Creates a new a matrix that is a copy of the existing one.                *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix of which to make a copy.                          *
 *                                                                           *
 * Returns: A pointer to the new matrix copy.                                *
 *****************************************************************************/
matrix_t *mtrx_copy(matrix_t *matrix);


/*---------------------------- Deallocation ---------------------------------*/

/*****************************************************************************
 * Frees the memory associated with the matrix.                              *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix to deallocate.                                    *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void mtrx_destroy(matrix_t *matrix);


/*-------------------------------- Display ----------------------------------*/

/*****************************************************************************
 * Prints a matrix to standard output.                                       *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix to print.                                         *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void mtrx_print(matrix_t *matrix);


/*------------------------------ Comparisons --------------------------------*/

/*****************************************************************************
 * Checks if the matrix is square. In other words, if the number of rows is  *
 * the same as the number columns in the matrix.                             *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix to check for squareness.                          * *                                                                           *
 * Returns: True if the matrix is square, false otherwise.                   *
 *****************************************************************************/
bool mtrx_is_sqr(matrix_t *matrix);


/*****************************************************************************
 * Checks if two matrices have equal dimensions. In other words, checks if   *
 * the matrices have the same number of rows and the same number of columns. *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first matrix.                                                 *
 *     B - The second matrix.                                                *
 *                                                                           *
 * Returns: True if the matrix have equal dimensions, false otherwise.       *
 *****************************************************************************/
bool mtrx_eq_dim(matrix_t *A, matrix_t *B);


/*****************************************************************************
 * Checks if two matrices are equal. To be equal, the matrices must have the *
 * same dimensions and contain all of the same values in the same positions. *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first matrix.                                                 *
 *     B - The second matrix.                                                *
 *                                                                           *
 * Returns: True if the matrix have equal dimensions, false otherwise.       *
 *****************************************************************************/
bool mtrx_eq(matrix_t *A, matrix_t *B);


/*------------------------------ Max and Min --------------------------------*/

/*****************************************************************************
 * Returns the maximum value in the matrix.                                  *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix in which to find the maximum value.               *
 *                                                                           *
 * Returns: The maximum value in the matrix.                                 *
 *****************************************************************************/
scalar_t mtrx_max(matrix_t *matrix);


/*****************************************************************************
 * Returns the minimum value in the matrix.                                  *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix in which to find the minimum value.               *
 *                                                                           *
 * Returns: The minimum value in the matrix.                                 *
 *****************************************************************************/
scalar_t mtrx_min(matrix_t *matrix);


/*--------------------------- Matrix Arithmetic -----------------------------*/

/*****************************************************************************
 * Adds to matrices together.                                                *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first matrix.                                                 *
 *     B - The second matrix.                                                *
 *                                                                           *
 * Returns: A new matrix that is the result of B added to A.                 *
 *****************************************************************************/
matrix_t *mtrx_add(matrix_t *A, matrix_t *B);


/*****************************************************************************
 * Subtracts one matrix from another.                                        *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first matrix.                                                 *
 *     B - The second matrix, which is subtracted from the first.            *
 *                                                                           *
 * Returns: A new matrix which is equal to B subtracted from A.              *
 *****************************************************************************/
matrix_t *mtrx_subtract(matrix_t *A, matrix_t *B);


/*****************************************************************************
 * Scales every element in the matrix by a scalar value.                     *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix to be scaled.                                     *
 *     scalar - The value by which to scale.                                 *
 *                                                                           *
 * Returns: A new scaled matrix.                                             *
 *****************************************************************************/
matrix_t *mtrx_scale(matrix_t *matrix, scalar_t scalar);


/*****************************************************************************
 * Multiplies two matrices together.                                         *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first matrix, of size m x n.                                  *
 *     B - The second matrix, of size n x p.                                 *
 *                                                                           *
 * Returns: A m x p matrix that is the result of multiplying A and B.        *
 *****************************************************************************/
matrix_t *mtrx_mult(matrix_t *A, matrix_t *B);


/*****************************************************************************
 * Multiplies a matrix by a vector.                                          *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix, of size m x n.                                   *
 *     vector - The vector, of length n.                                     *
 *                                                                           *
 * Returns: A new vector that is the result of the matrix multipled by the   *
 *     vector.                                                               *
 *****************************************************************************/
vector_t *mtrx_mult_vctr(matrix_t *matrix, vector_t *vector);


/*------------------------- Pointwise Arithmetic ----------------------------*/

/*****************************************************************************
 * Performs pointwise multiplication of two matrices. The resulting matrix   *
 * consists of elements that are the product of the elements in the same     *
 * positions in the original two matrices. The two matrices being multiplied *
 * must be of the same size.                                                 *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first matrix.                                                 *
 *     B - The second matrix.                                                *
 *                                                                           *
 * Returns: A new matrix consisting of elements that are the result of       *
 *     multplying the elements in the same positions of the original two     *
 *     matrices.                                                             *
 *****************************************************************************/
matrix_t *mtrx_pw_mult(matrix_t *A, matrix_t *B);


/*****************************************************************************
 * Performs pointwise division of two matrices. The resulting matrix         *
 * consists of elements that are the quotient of the element in the same     *
 * position in the first matrix divided by the element in the same position  *
 * in the second matrix. The two matrices being divided must be of the same  *
 * size.                                                                     *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first matrix. The elements of this matrix are the dividends.  *
 *     B - The second matrix. The elements of this matrix are the divisors.  *
 *                                                                           *
 * Returns: A new matrix consisting of elements that are the result of       *
 *     multplying the elements in the same positions of the original two     *
 *     matrices.                                                             *
 *****************************************************************************/
matrix_t *mtrx_pw_div(matrix_t *A, matrix_t *B);


/*****************************************************************************
 * Performs pointwise exponentiation of two matrices. The resulting matrix   *
 * consists of elements that are the result of the element in the same       *
 * position in the first matrix raised to the power of the element in the    *
 * same position in the second matrix.The two matrices being multiplied must *
 * be of the same size.                                                      *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The first matrix. The elements of this matrix are the bases.      *
 *     B - The second matrix. The elements of this matrix are the exponents. *
 *                                                                           *
 * Returns: A new matrix consisting of elements that are the result of       *
 *     multplying the elements in the same positions of the original two     *
 *     matrices.                                                             *
 *****************************************************************************/
matrix_t *mtrx_pw_pow(matrix_t *A, matrix_t *B);


/*-------------------------- Matrix Manipulation ----------------------------*/

/*****************************************************************************
 * Swaps two rows in the matrix.                                             *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix in which rows are to be swapped.                  *
 *     r1 - The index of the first row to be swapped.                        *
 *     r2 - The index of the row with which r1 is to be swapped.             *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void mtrx_row_swap(matrix_t *matrix, size_t r1, size_t r2);


/*****************************************************************************
 * Swaps two columns the matrix.                                             *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix in which rows are to be swapped.                  *
 *     c1 - The index of the first column to be swapped.                     *
 *     c2 - The index of the column with which c1 is to be swapped.          *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void mtrx_col_swap(matrix_t *matrix, size_t c1, size_t c2);


/*****************************************************************************
 * Multiplies each element in a single row of the matrix by a scalar.        *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix in which rows are to be swapped.                  *
 *     row - The index of the row to be scaled.                              *
 *     scalar - The scalar by which to multiply each element in the row at   *
 *              index row.                                                   *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void mtrx_scale_row(matrix_t *matrix, size_t row, scalar_t scalar);


/*****************************************************************************
 * Multiplies each element in a single column of the matrix by a scalar.     *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix in which columnss are to be swapped.              *
 *     col - The index of the column to be scaled.                           *
 *     scalar - The scalar by which to multiply each element in the column   *
 *              at index col.                                                *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void mtrx_scale_col(matrix_t *matrix, size_t col, scalar_t scalar);


/*----------------- Row and Column Accessors and Mutators -------------------*/

/*****************************************************************************
 * Returns the row of a matrix as a vector.                                  *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix from which to get the column.                     *
 *     row - The index of the row to return.                                 *
 *                                                                           *
 * Returns: A vector containing the values of the row.                       *
 *****************************************************************************/
vector_t *mtrx_get_row(matrix_t *matrix, size_t row);


/*****************************************************************************
 * Returns the column of a matrix as a vector.                               *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix from which to get the column.                     *
 *     col - The index of the column to return.                              *
 *                                                                           *
 * Returns: A vector containing the values of the column.                    *
 *****************************************************************************/
vector_t *mtrx_get_col(matrix_t *matrix, size_t col);


/*****************************************************************************
 * Sets the values in a row of the matrix to those in a given vector.        *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix in which to set the row.                          *
 *     vector - The vector containing the values to copy into the specified  *
 *              row.                                                         *
 *     row - The index of the row to set.                                    *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void mtrx_set_row(matrix_t *matrix, vector_t *vector, size_t row);


/*****************************************************************************
 * Sets the values in a column of the matrix to those in a given vector.     *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix in which to set the column.                       *
 *     vector - The vector containing the values to copy into the specified  *
 *              column.                                                      *
 *     col - The index of the column to set.                                 *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void mtrx_set_col(matrix_t *matrix, vector_t *vector, size_t col);


/*----------------------------- Sub-matrices --------------------------------*/

/*****************************************************************************
 * Copies the values in a block of an existing matrix into a new matrix.     *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix from which to get the submatrix.                  *
 *     rows - The indexes of the rows in the original matrix to include in   *
 *            the submatrix.                                                 *
 *     cols - The indexes of the columns in the original matrix to include   *
 *            in the submatrix.                                              *
 *                                                                           *
 * Returns: A new matrix representing the submatrix.                         *
 *****************************************************************************/
matrix_t *mtrx_sub_matrix(matrix_t *matrix, indexer_t *rows, indexer_t *cols);


/*****************************************************************************
 * Copies the values in a block of an existing matrix into a new matrix.     *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix from which to get the block.                      *
 *     r1 - The first row contained in the block, inclusive.                 *
 *     r2 - The row immediately after the last row contained in the block.   *
 *     c1 - The first column contained in the block, inclusive.              *
 *     c2 - The column immediately after the last column contained in the    *
 *          block.                                                           *
 *                                                                           *
 * Returns: A new matrix representing the block.                             *
 *****************************************************************************/
matrix_t *mtrx_sub_block(matrix_t *matrix, size_t r1, size_t r2, size_t c1,
    size_t c2);


/*------------------------------- Operations --------------------------------*/

/*****************************************************************************
 * Creates a new matrix that is the transpose of the original matrix.        *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix of which to take the transpose.                   *
 *                                                                           *
 * Returns: A new matrix that is the transpose of the original matrix. If    *
 *     the original matrix was of size m x n, this new matrix will be n x m. *
 *****************************************************************************/
matrix_t *mtrx_transpose(matrix_t *matrix);


/*****************************************************************************
 * Calculates the determinant of the matrix.                                 *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix of which to calculate the determinant.            *
 *                                                                           *
 * Returns: The value of the determinant of the matrix.                      *
 *****************************************************************************/
scalar_t mtrx_det(matrix_t *matrix);


/*****************************************************************************
 * Finds the inverse of a matrix, if it is invertible. The matrix must be    *
 * square.                                                                   *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix of which to find the inverse.                     *
 *                                                                           *
 * Returns: The inverse matrix. If the matrix is not square or is not        *
 *     invertible, returns NULL instead.                                     *
 *****************************************************************************/
matrix_t *mtrx_inv(matrix_t *matrix);


/*----------------------------- System Solving ------------------------------*/

/*****************************************************************************
 * Solves a system of linear equations of the from Ax = b, where A is a      *
 * matrix and b is a vector, and x is the unknown solution vector.           *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The coefficient matrix.                                           *
 *     B - The constant vector.                                              *
 *                                                                           *
 * Returns: A vector representing the solution to the linear system of       *
 *     equations.                                                            *
 *****************************************************************************/
vector_t *mtrx_solve(matrix_t *matrix, vector_t *vector);


/*****************************************************************************
 * Checks if a matrix is diagonally dominant.                                *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix to check for diagonal dominance.                  *
 *                                                                           *
 * Returns: True if the matrix is diagonally dominant, false otherwise.      *
 *****************************************************************************/
bool mtrx_is_diag_dom(matrix_t *matrix);


/*****************************************************************************
 * Performs matrix operations to make the given matrix diagonally dominant.  *
 * Some matrices cannot be made diagonally dominant, in which case the       *
 * function will return false.                                               *
 *                                                                           *
 * Fields:                                                                   *
 *     matrix - The matrix to make diagonally dominant.                      *
 *                                                                           *
 * Returns: True if the matrix was successfully made diagonally dominant,    *
 *     false otherwise.                                                      *
 *****************************************************************************/
bool mtrx_make_diag_dom(matrix_t *matrix);


/*****************************************************************************
 * Solves a system of linear equations of the from Ax = b, where A is a      *
 * matrix and b is a vector, and x is the unknown solution vector,           *
 * numerically.                                                              *
 *                                                                           *
 * Fields:                                                                   *
 *     A - The coefficient matrix.                                           *
 *     B - The constant vector.                                              *
 *     X - The solution vector. Initially, this vector should contain an     *
 *         initial guess for the solution.                                   *
 *     tolerance - The maximum difference in the true value of each element  *
 *                 of the solution and the calulated value.                  *
 *                                                                           *
 * Returns: A vector representing the solution to the linear system of       *
 *     equations, with each element being within the acceptable tolerance.   *
 *****************************************************************************/
void mtrx_solve_gs(matrix_t *A, vector_t *B, vector_t *X, scalar_t tolerance);

#endif
