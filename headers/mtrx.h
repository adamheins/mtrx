/*****************************************************************************
 * A simple linear algebra library containing functions for vector and       *
 * matrix operations.                                                        *
 *                                                                           *
 * Author: Adam Heins                                                        *
 *                                                                           *
 * Last Updated: 2014-11-29                                                  *
 *****************************************************************************/

#ifndef MTRX_H_
#define MTRX_H_


/*===========================================================================*
 *                               STRUCTURES                                  *
 *===========================================================================*/


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
	double **values;
	unsigned int rows;
	unsigned int columns;
} matrix_t;


/*****************************************************************************
 * Defines a 1-dimensional vector.                                           *
 *                                                                           *
 * Fields:                                                                   *
 *     values - The array of double precision values stored in the vector.   *
 *     length - The length of the vector.                                    *
 *****************************************************************************/
typedef struct vector {
	double *values;
	unsigned int length;
} vector_t;


/*****************************************************************************
 * Defines an indexing array. This structure contains a series of integer    *
 * indexes for accessing elements in vectors and matrices.                   *
 *                                                                           *
 * Fields:                                                                   *
 *     values - The indexes contained by the indexer.                        *
 *     length - The length of the indexer.                                   *
 *****************************************************************************/             
typedef struct indexer {
	unsigned int *values;
	unsigned int length;
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
vector_t *vctr_empty(unsigned int);


/*****************************************************************************
 * Creates a vector of specified length with all values set to zero.         *
 *                                                                           *
 * Fields:                                                                   *
 *     length - The length of the vector being created.                      *
 *                                                                           *
 * Returns: A pointer to the created vector.                                 *
 *****************************************************************************/
vector_t *vctr_zeros(unsigned int);


/*****************************************************************************
 * Creates a vector of specified length with all values set to one.          *
 *                                                                           *
 * Fields:                                                                   *
 *     length - The length of the vector being created.                      *
 *                                                                           *
 * Returns: A pointer to the created vector.                                 *
 *****************************************************************************/
vector_t *vctr_ones(unsigned int);


/*****************************************************************************
 * Creates a new vector that is a copy of the existing one.                  *
 *                                                                           *
 * Fields:                                                                   *
 *     vector - The vector to be copied.                                     *
 *                                                                           *
 * Returns: A pointer to the new copied vector.                              *
 *****************************************************************************/
vector_t *vctr_copy(vector_t *);


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
void vctr_destroy(vector_t *);


/*-------------------------------- Display ----------------------------------*/

/*****************************************************************************
 * Prints the values contained in the vector to standard output.             *
 *                                                                           *
 * Fields:                                                                   *
 *     vector - The vector to print.                                         *
 *                                                                           *
 * Returns: void                                                             *
 *****************************************************************************/
void vctr_print(vector_t *);


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
int vctr_eq_len(vector_t *, vector_t *);


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
int vctr_eq(vector_t *, vector_t *);


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
double vctr_dot_prod(vector_t *, vector_t *);


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
vector_t *vctr_cross_prod(vector_t *, vector_t *);


/*****************************************************************************
 * Calculates the magnitude of a vector.                                     *
 *                                                                           *
 * Fields:                                                                   *
 *     vector - The vector for which the magnitude is being calculated.      * 
 *                                                                           *
 * Returns: The magnitude of the vector.                                     *
 *****************************************************************************/
double vctr_mag(vector_t *);


/*===========================================================================*
 *                           MATRIX FUNCTIONS                                *
 *===========================================================================*/


/*--------------------------- Initialization --------------------------------*/

matrix_t *mtrx_zeros(unsigned int, unsigned int);
matrix_t *mtrx_ones(unsigned int, unsigned int);
matrix_t *mtrx_id(unsigned int);
matrix_t *mtrx_diag(vector_t *);
matrix_t *mtrx_copy(matrix_t *);


/*---------------------------- Deallocation ---------------------------------*/

void mtrx_destroy(matrix_t *);


/*-------------------------------- Display ----------------------------------*/

void mtrx_print(matrix_t *);


/*------------------------------ Comparisons --------------------------------*/

int mtrx_is_sqr(matrix_t *);
int mtrx_eq_dim(matrix_t *, matrix_t *);
int mtrx_eq(matrix_t *, matrix_t *);

double mtrx_max(matrix_t *);
double mtrx_min(matrix_t *);


/*--------------------------- Matrix Arithmetic -----------------------------*/

matrix_t *mtrx_mult(matrix_t *, matrix_t *);
matrix_t *mtrx_scale(matrix_t *, double);
matrix_t *mtrx_add(matrix_t *, matrix_t *);
matrix_t *mtrx_subtract(matrix_t *, matrix_t *);
vector_t *mtrx_mult_vctr(matrix_t *, vector_t *);


/*------------------------- Pointwise Arithmetic ----------------------------*/

matrix_t *mtrx_pw_mult(matrix_t *, matrix_t *);
matrix_t *mtrx_pw_div(matrix_t *, matrix_t *);
matrix_t *mtrx_pw_pow(matrix_t *, matrix_t *);


/*-------------------------- Matrix Manipulation ----------------------------*/

void mtrx_row_swap(matrix_t *, unsigned int, unsigned int);
void mtrx_col_swap(matrix_t *, unsigned int, unsigned int);
void mtrx_scale_row(matrix_t *, unsigned int, double);
void mtrx_scale_col(matrix_t *, unsigned int, double);


/*----------------- Row and Column Accessors and Mutators -------------------*/

vector_t *mtrx_get_row(matrix_t *, unsigned int);
vector_t *mtrx_get_col(matrix_t *, unsigned int);
void mtrx_set_col(matrix_t *, vector_t *, unsigned int);
void mtrx_set_row(matrix_t *, vector_t *, unsigned int);


/*----------------------------- Sub-matrices --------------------------------*/

matrix_t *mtrx_sub_block(matrix_t *, unsigned int, unsigned int, unsigned int, unsigned int);
matrix_t *mtrx_submatrix(matrix_t *, indexer_t *, indexer_t *);


/*------------------------------- Operations --------------------------------*/

matrix_t *mtrx_transpose(matrix_t *);
double mtrx_det(matrix_t *);
matrix_t *mtrx_inv(matrix_t *);

/*----------------------------- System Solving ------------------------------*/
vector_t *mtrx_solve(matrix_t *, vector_t *);
void mtrx_solve_gs(matrix_t *, vector_t *, vector_t *, double);

#endif
