#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "mtrx.h"


//#define NULL (void *)0;
typedef unsigned int U32;
typedef int S32;

/*=========================================================================================*
 *                             HELPER FUNCTIONS & STRUCTS                                  *
 *=========================================================================================*/

/* Stores the results of an LU decomposition. */
typedef struct {
	matrix_t *L;
	matrix_t *U;
	matrix_t *P;
	unsigned int swaps;
} lu_factors_t;


/* Frees the memory associated with an lu_factors_t object. */
void lu_factors_destroy(lu_factors_t *lu_factors) {
	mtrx_destroy(lu_factors->L);
	mtrx_destroy(lu_factors->U);
	mtrx_destroy(lu_factors->P);
	free(lu_factors);
}


/*=========================================================================================*
 *                                  VECTOR FUNCTIONS                                       *
 *=========================================================================================*/

vector_t *vctr_empty(unsigned int length) {

	vector_t *vector = (vector_t *)malloc(sizeof(vector_t));

	vector->values = (double *)malloc(length * sizeof(double));
	vector->length = length;

	return vector;
}


vector_t *vctr_zeros(unsigned int length) {

	vector_t *vector = (vector_t *)malloc(sizeof(vector_t));

	vector->values = (double *)calloc(length, sizeof(double));
	vector->length = length;

	return vector;
}


vector_t *vctr_ones(unsigned int length) {

	vector_t *vector = (vector_t *)malloc(sizeof(vector_t));

	vector->values = (double *)malloc(length * sizeof(double));
	vector->length = length;

	for (unsigned int i = 0; i < length; ++i)
		vector->values[i] = 1;

	return vector;
}


/* Free the memory of the vector object. */
void vctr_destroy(vector_t *vector) {
	free(vector->values);
	free(vector);
	vector = NULL;
}


vector_t *vctr_copy(vector_t *vector) {
	vector_t *copy = vctr_empty(vector->length);
	for (unsigned int i = 0; i < vector->length; ++i)
		copy->values[i] = vector->values[i];
	return copy;
}


/* Prints out the vector to standard output. */
void vctr_print(vector_t *vector) {
	for (unsigned int i = 0; i < vector->length; ++i)
		printf("%f\n", vector->values[i]);
}


/* Returns true if the two vectors are of equal length, false otherwise. */
int vctr_eq_len(vector_t *A, vector_t *B) {
	return (A->length == B->length);
}


/* Calculates the dot product of two vectors. */
double vctr_dot_prod(vector_t *A, vector_t *B) {

	// Check that the vectors are of equal length.
	if (!vctr_eq_len(A, B))
		return 0;

	double prod = 0;
	for (unsigned int i = 0; i < A->length; ++i)
		prod += A->values[i] * B->values[i];

	return prod;
}


/* Calculates the cross product of two 3-dimensional vectors. */
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


/* Calculates the magnitude of a vector. */
double vctr_mag(vector_t *vector) {

	double mag = 0;

	for (unsigned int i = 0; i < vector->length; ++i)
		mag += vector->values[i] * vector->values[i];

	return sqrt(mag);
}


/*=========================================================================================*
 *                                  MATRIX FUNCTIONS                                       *
 *=========================================================================================*/


matrix_t *mtrx_zeros(unsigned int rows, unsigned int cols) {

	matrix_t *matrix = (matrix_t *)malloc(sizeof(matrix_t));

	// Allocate the array of rows.
	matrix->values = (double **)calloc(rows, sizeof(double *));

	// Allocate the columns arrays.
	unsigned int i;
	for (i = 0; i < rows; ++i)
		matrix->values[i] = (double *)calloc(cols, sizeof(double));

	matrix->rows = rows;
	matrix->columns = cols;

	return matrix;
}


matrix_t *mtrx_ones(unsigned int rows, unsigned int cols) {

	matrix_t *matrix = (matrix_t *)malloc(sizeof(matrix_t));

	// Allocate the array of rows.
	matrix->values = (double **)malloc(rows * sizeof(double *));

	// Allocate the columns arrays.
	unsigned int i;
	for (i = 0; i < rows; ++i)
		matrix->values[i] = (double *)malloc(cols * sizeof(double));

	matrix->rows = rows;
	matrix->columns = cols;

	for (unsigned int r = 0; r < rows; ++r) {
		for (unsigned int c = 0; c < cols; ++c)
			matrix->values[r][c] = 1;
	}

	return matrix;
}


/*
* Generates an identity matrix (1's on the diagonal, 0's elsewhere) of size n.
*/
matrix_t *mtrx_id(unsigned int n) {

	// Initialize empty matrix.
	matrix_t *id = mtrx_zeros(n, n);

	// Fill diagonal with zeros.
	for (unsigned int i = 0; i < n; ++i)
		id->values[i][i] = 1;

	return id;
}


matrix_t *mtrx_diag(vector_t *vector) {

	matrix_t *matrix = mtrx_zeros(vector->length, vector->length);

	for (unsigned int i = 0; i < vector->length; ++i)
		matrix->values[i][i] = vector->values[i];

	return matrix;
}


matrix_t *mtrx_rnd(U32 rows, U32 cols, U32 max) {
	matrix_t *matrix = mtrx_zeros(rows, cols);

	//srand(time(0));

	for (U32 i = 0; i < rows; ++i) {
		for (U32 j = 0; j < cols; ++j) {
			matrix->values[i][j] = rand() % max;
		}
	}

	return matrix;
}


/* Free the memory associated with the matrix object. */
void mtrx_destroy(matrix_t *matrix) {
	unsigned int i;
	for (i = 0; i < matrix->rows; ++i)
		free(matrix->values[i]);
	free(matrix->values);
	free(matrix);
	matrix = NULL;
}


/* Copies an existing matrix into a new one. */
matrix_t *mtrx_copy(matrix_t *matrix) {

	matrix_t *copy = mtrx_zeros(matrix->rows, matrix->columns);

	for (unsigned int i = 0; i < matrix->rows; ++i) {
		for (unsigned int j = 0; j < matrix->columns; ++j)
			copy->values[i][j] = matrix->values[i][j];
	}

	return copy;
}


/* Print out the matrix to standard output. */
void mtrx_print(matrix_t *matrix) {
	for (unsigned int i = 0; i < matrix->rows; ++i) {
		for (unsigned int j = 0; j < matrix->columns; ++j)
			printf("%f ", matrix->values[i][j]);
		printf("\n");
	}
}


/* Returns true if the matrix is square, false otherwise. */
int mtrx_is_sqr(matrix_t *matrix) {
	return matrix->rows == matrix->columns;
}


/* Returns true if the two matrices have equal dimensions, false otherwise. */
int mtrx_eq_dim(matrix_t *A, matrix_t *B) {
	return (A->rows == B->rows && A->columns == B->columns);
}


/* Returns the maximum value in the matrix. */
double mtrx_max(matrix_t *matrix) {

	double largest = matrix->values[0][0];

	for (unsigned int i = 0; i < matrix->rows; ++i) {
		for (unsigned int j = 0; j < matrix->columns; ++j) {
			if (matrix->values[i][j] > largest)
				largest = matrix->values[i][j];
		}
	}

	return largest;
}


/* Returns the minimum value in the matrix. */
double mtrx_min(matrix_t *matrix) {

	double smallest = matrix->values[0][0];

	for (unsigned int i = 0; i < matrix->rows; ++i) {
		for (unsigned int j = 0; j < matrix->columns; ++j) {
			if (matrix->values[i][j] < smallest)
				smallest = matrix->values[i][j];
		}
	}

	return smallest;
}


/* Naive multiplication algorithm. */
matrix_t *mtrx_naive_mult(matrix_t *A, matrix_t *B) {

	// Create a new matrix to hold the result.
	matrix_t *C = mtrx_zeros(A->rows, B->columns);
	
	for (unsigned int i = 0; i < B->columns; ++i) {
		for (unsigned int j = 0; j < A->rows; ++j) {
			for (unsigned int k = 0; k < A->columns; ++k)
				C->values[j][i] += A->values[j][k] * B->values[k][i];
		}
	}
	return C;
}


/*
 * Splits a block from an existing matrix, and returns it as a new matrix.
 */
matrix_t *mtrx_split(matrix_t *matrix, U32 r1, U32 r2, U32 c1, U32 c2) {

	U32 num_rows = r2 - r1;
	U32 num_cols = c2 - c1;
	matrix_t *child_matrix = mtrx_zeros(num_rows, num_cols);

	for (U32 i = 0; i < num_rows; ++i) {
		for (U32 j = 0; j < num_cols; ++j)
			child_matrix->values[i][j] = matrix->values[r1 + i][c1 + j];
	}

	return child_matrix;
}


/*
 * Splits a block from an existing matrix, and returns it as a new matrix. Also pads
 * the new block matrix to a specified number of rows and columns.
 */
matrix_t *mtrx_padded_split(matrix_t *matrix, U32 r1, U32 r2, U32 c1, U32 c2, U32 rows, U32 cols) {

	U32 num_rows = r2 - r1;
	U32 num_cols = c2 - c1;
	matrix_t *child_matrix = mtrx_zeros(rows, cols);

	for (U32 i = 0; i < num_rows; ++i) {
		for (U32 j = 0; j < num_cols; ++j)
			child_matrix->values[i][j] = matrix->values[r1 + i][c1 + j];
	}

	return child_matrix;
}


/* Fast multiplication algorithm. */
matrix_t *mtrx_strassen_mult(matrix_t *A, matrix_t *B) {

	// Ending condition.
	if (A->rows <= 1 || A->columns <= 1)
		return mtrx_naive_mult(A, B);

	U32 Ahr = round(A->rows / 2.0);
	U32 Ahc = round(A->columns / 2.0);
	U32 Bhr = round(B->rows / 2.0);
	U32 Bhc = round(B->columns / 2.0);

	// Split A into four submatrices.
	matrix_t *A11 = mtrx_split(A, 0, Ahr, 0, Ahc);
	matrix_t *A12 = mtrx_padded_split(A, 0, Ahr, Ahc, A->columns, Ahr, Ahc);
	matrix_t *A21 = mtrx_padded_split(A, Ahr, A->rows, 0, Ahc, Ahr, Ahc);
	matrix_t *A22 = mtrx_padded_split(A, Ahr, A->rows, Ahc, A->columns, Ahr, Ahc);

	// Split B into four submatrices.
	matrix_t *B11 = mtrx_split(B, 0, Bhr, 0, Bhc);
	matrix_t *B12 = mtrx_padded_split(B, 0, Bhr, Bhc, B->columns, Bhr, Bhc);
	matrix_t *B21 = mtrx_padded_split(B, Bhr, B->rows, 0, Bhc, Bhr, Bhc);
	matrix_t *B22 = mtrx_padded_split(B, Bhr, B->rows, Bhc, B->columns, Bhr, Bhc);

	// Calculate M matrices. This is the recursive section of the algorithm.
	matrix_t *M1 = mtrx_strassen_mult(mtrx_add(A11, A22), mtrx_add(B11, B22));
	matrix_t *M2 = mtrx_strassen_mult(mtrx_add(A21, A22), B11);
	matrix_t *M3 = mtrx_strassen_mult(A11, mtrx_subtract(B12, B22));
	matrix_t *M4 = mtrx_strassen_mult(A22, mtrx_subtract(B21, B11));
	matrix_t *M5 = mtrx_strassen_mult(mtrx_add(A11, A12), B22);
	matrix_t *M6 = mtrx_strassen_mult(mtrx_subtract(A21, A11), mtrx_add(B11, B12));
	matrix_t *M7 = mtrx_strassen_mult(mtrx_subtract(A12, A22), mtrx_add(B21, B22));

	mtrx_destroy(A11);
	mtrx_destroy(A12);
	mtrx_destroy(A21);
	mtrx_destroy(A22);

	mtrx_destroy(B11);
	mtrx_destroy(B12);
	mtrx_destroy(B21);
	mtrx_destroy(B22);

	// Calculate each quadrant of C from M matrices.
	matrix_t *C11 = mtrx_add(mtrx_subtract(mtrx_add(M1, M4), M5), M7);
	matrix_t *C12 = mtrx_add(M3, M5);
	matrix_t *C21 = mtrx_add(M2, M4);
	matrix_t *C22 = mtrx_add(mtrx_add(mtrx_subtract(M1, M2), M3), M6);

	mtrx_destroy(M1);
	mtrx_destroy(M2);
	mtrx_destroy(M3);
	mtrx_destroy(M4);
	mtrx_destroy(M5);
	mtrx_destroy(M6);
	mtrx_destroy(M7);

	matrix_t *C = mtrx_zeros(A->rows, B->columns);

	for (U32 i = 0; i < C11->rows; ++i) {
		for (U32 j = 0; j < C11->columns; ++j)
			C->values[i][j] = C11->values[i][j];
	}

	for (U32 i = 0; i < C11->rows; ++i) {
		for (U32 j = 0; j < C->columns / 2; ++j)
			C->values[i][j + C11->columns] = C12->values[i][j];
	}

	for (U32 i = 0; i < C->rows / 2; ++i) {
		for (U32 j = 0; j < C11->columns; ++j)
			C->values[i + C11->rows][j] = C21->values[i][j];
	}

	for (U32 i = 0; i < C->rows / 2; ++i) {
		for (U32 j = 0; j < C->columns / 2; ++j)
			C->values[i + C11->rows][j + C11->columns] = C22->values[i][j];
	}

	mtrx_destroy(C11);
	mtrx_destroy(C12);
	mtrx_destroy(C21);
	mtrx_destroy(C22);

	return C;
}


/* Multiply two matrices together. */
matrix_t *mtrx_mult(matrix_t *A, matrix_t *B) {
	
	// Check for matrices inappropriately sized for multiplication.
	if (A->columns != B->rows)
		return NULL;

	// Currently defaulting to naive multiplication.
	return mtrx_naive_mult(A, B);
}


/* Scale the entire matrix by a scalar constant. */
matrix_t *mtrx_scale(matrix_t *matrix, double scalar) {

	matrix_t *multiple = mtrx_zeros(matrix->rows, matrix->columns);

	for (unsigned int i = 0; i < matrix->rows; ++i) {
		for (unsigned int j = 0; j < matrix->columns; ++j)
			multiple->values[i][j] = matrix->values[i][j] * scalar;
	}

	return multiple;
}


/*
* Add two matrices.
*/
matrix_t *mtrx_add(matrix_t *A, matrix_t *B) {

	// Check that matrices have the same dimensions.
	if (!mtrx_eq_dim(A, B))
		return NULL;

	// Create a new matrix to hold the result.
	matrix_t *C = mtrx_zeros(A->rows, A->columns);

	for (unsigned int i = 0; i < A->rows; ++i) {
		for (unsigned int j = 0; j < A->columns; j++)
			C->values[i][j] = A->values[i][j] + B->values[i][j];
	}

	return C;
}


/*
* Subtract a matrix from another one.
*/
matrix_t *mtrx_subtract(matrix_t *A, matrix_t *B) {

	// Check that matrices have the same dimensions.
	if (!mtrx_eq_dim(A, B))
		return NULL;

	// Create a new matrix to hold the result.
	matrix_t *C = mtrx_zeros(A->rows, A->columns);

	// Subtract each element in B from the corresponding element in A to
	// generate the resulting element in C.
	for (unsigned int i = 0; i < A->rows; ++i) {
		for (unsigned int j = 0; j < A->columns; j++)
			C->values[i][j] = A->values[i][j] - B->values[i][j];
	}

	return C;
}


/* Multiplies an (n x m) matrix with an (m x 1) vector. */
vector_t *mtrx_mult_vctr(matrix_t *A, vector_t *B) {

	if (A->columns != B->length)
		return NULL;

	vector_t *C = vctr_zeros(A->rows);

	for (unsigned int i = 0; i < A->rows; ++i) {
		for (unsigned int j = 0; j < A->columns; ++j)
			C->values[i] += A->values[i][j] * B->values[j];
	}

	return C;
}


/*
* Pointwise multiplication of a matrix.
* Multiplies each element of A by the corresponding element of B.
*/
matrix_t *mtrx_pw_mult(matrix_t *A, matrix_t *B) {

	// Check that matrices have the same dimensions.
	if (!mtrx_eq_dim(A, B))
		return NULL;

	matrix_t *C = mtrx_zeros(A->rows, A->columns);

	for (unsigned int r = 0; r < A->rows; ++r) {
		for (unsigned int c = 0; c < A->columns; ++c)
			C->values[r][c] = A->values[r][c] * B->values[r][c];
	}

	return C;
}


/*
* Pointwise division of a matrix.
* Divides each element of A by the corresponding element of B.
*/
matrix_t *mtrx_pw_div(matrix_t *A, matrix_t *B) {

	// Check that matrices have the same dimensions.
	if (!mtrx_eq_dim(A, B))
		return NULL;

	matrix_t *C = mtrx_zeros(A->rows, A->columns);

	for (unsigned int r = 0; r < A->rows; ++r) {
		for (unsigned int c = 0; c < A->columns; ++c)
			C->values[r][c] = A->values[r][c] / B->values[r][c];
	}

	return C;
}


/* 
 * Pointwise exponentiation of a matrix.
 * Raises each element of A to the power of the corresponding element of B.
 */
matrix_t *mtrx_pw_pow(matrix_t *A, matrix_t *B) {

	// Check that matrices have the same dimensions.
	if (!mtrx_eq_dim(A, B))
		return NULL;

	matrix_t *C = mtrx_zeros(A->rows, A->columns);

	for (unsigned int r = 0; r < A->rows; ++r) {
		for (unsigned int c = 0; c < A->columns; ++c)
			C->values[r][c] = pow(A->values[r][c], B->values[r][c]);
	}

	return C;
}


/* Swap two rows in a matrix. */
void mtrx_row_swap(matrix_t *matrix, unsigned int r1, unsigned int r2) {

	double temp;

	for (unsigned int i = 0; i < matrix->columns; ++i) {
		temp = matrix->values[r1][i];
		matrix->values[r1][i] = matrix->values[r2][i];
		matrix->values[r2][i] = temp;
	}	
}


/* Swap two columns in a matrix. */
void mtrx_col_swap(matrix_t *matrix, unsigned int c1, unsigned int c2) {
	
	double temp;

	for (unsigned int i = 0; i < matrix->rows; ++i) {
		temp = matrix->values[i][c1];
		matrix->values[i][c1] = matrix->values[i][c2];
		matrix->values[1][c2] = temp;
	}
}


/* Multiplies a row in the matrix by a scalar value. */
void mtrx_scale_row(matrix_t *matrix, unsigned int row, double scalar) {
	for (unsigned int i = 0; i < matrix->columns; ++i)
		matrix->values[row][i] *= scalar;
}


/* Multiplies a column in the matrix by a scalar value. */
void mtrx_scale_col(matrix_t *matrix, unsigned int col, double scalar) {
	for (unsigned int i = 0; i < matrix->rows; ++i)
		matrix->values[i][col] *= scalar;
}


/* Gets the row of a matrix as a vector. */
vector_t *mtrx_get_row(matrix_t *matrix, unsigned int row) {
	
	vector_t *vector = vctr_empty(matrix->columns);

	for (unsigned int i = 0; i < vector->length; ++i)
		vector->values[i] = matrix->values[row][i];

	return vector;
}


/* Gets the column of a matrix as a vector. */
vector_t *mtrx_get_col(matrix_t *matrix, unsigned int col) {
	
	vector_t *vector = vctr_empty(matrix->rows);

	for (unsigned int i = 0; i < vector->length; ++i)
		vector->values[i] = matrix->values[i][col];

	return vector;
}


/* Sets the values in the specified column of a matrix to those of a vector. */
void mtrx_set_col(matrix_t *matrix, vector_t *vector, unsigned int col) {

	// Check if vector is the same length as the number of rows in the matrix.
	if (matrix->rows != vector->length)
		return;

	// Replace values in matrix column with those of the vector.
	for (unsigned int i = 0; i < vector->length; ++i)
		matrix->values[i][col] = vector->values[i];
}


/* Sets the values in the specified row of a matrix to those of a vector. */
void mtrx_set_row(matrix_t *matrix, vector_t *vector, unsigned int row) {

	// Check if vector is the same length as the number of rows in the matrix.
	if (matrix->columns != vector->length)
		return;

	// Replace values in matrix row with those of the vector.
	for (unsigned int i = 0; i < vector->length; ++i)
		matrix->values[i][row] = vector->values[i];
}


/* Creates a submatrix from the given matrix. */
matrix_t *mtrx_submatrix(matrix_t *A, indexer_t *rows, indexer_t *columns) {

	matrix_t *sub = mtrx_zeros(rows->length, columns->length);

	for (unsigned int i = 0; i < rows->length; ++i) {
		for (unsigned int j = 0; j < columns->length; ++j)
			sub->values[i][j] = A->values[rows->values[i]][columns->values[j]];
	}

	return sub;
}


/* Performs pivoting to ensure a stable answer for LU decomposition. */
int lu_decomposition_pivot(matrix_t *L, matrix_t *P, unsigned int col) {

	unsigned int pivot = col;

	// Find largest absolute value in the column. That will be the pivot.
	for (unsigned int i = col + 1; i < L->rows; ++i) {
		if (fabs(L->values[i][col]) > fabs(L->values[pivot][col]))
			pivot = i;
	}

	// Swap rows such that the pivot is in the current position of [col, col].
	if (pivot != col) {
		mtrx_row_swap(L, col, pivot);
		mtrx_row_swap(P, col, pivot);
		return 1;
	}
	return 0;
}


// A is nxn and B is nx1
lu_factors_t *lu_decomposition(matrix_t *A) {

	unsigned int n = A->rows;
	unsigned int swaps = 0;

	matrix_t *L = mtrx_copy(A);
	matrix_t *U = mtrx_id(n);
	matrix_t *P = mtrx_id(n);

	swaps += lu_decomposition_pivot(L, P, 0);

	// First row of U.
	for (unsigned int j = 1; j < n; ++j)
		U->values[0][j] = L->values[0][j] / L->values[0][0];

	for (unsigned int j = 1; j < n - 1; ++j) {

		// Compute the next column of L.
		for (unsigned int i = j; i < n; ++i) {
			for (unsigned int k = 0; k < j; ++k)
				L->values[i][j] -= L->values[i][k] * U->values[k][j];
		}

		// Pivoting.
		swaps += lu_decomposition_pivot(L, P, j);

		// Compute the next row of U.
		for (unsigned int k = j + 1; k < n; ++k) {
			U->values[j][k] = L->values[j][k];
			for (unsigned int i = 0; i < j; ++i)
				U->values[j][k] -= L->values[j][i] * U->values[i][k];
			U->values[j][k] /= L->values[j][j];
		}
	}

	for (unsigned int k = 0; k < n - 1; ++k)
		L->values[n - 1][n - 1] -= L->values[n - 1][k] * U->values[k][n - 1];

	// Zero the upper triangle of L to make it a lower triangular matrix.
	for (unsigned int i = 0; i < n - 1; ++i) {
		for (unsigned int j = i + 1; j < n; j++)
			L->values[i][j] = 0;
	}

	lu_factors_t *lu_factors = (lu_factors_t *)malloc(sizeof(lu_factors_t));
	lu_factors->L = L;
	lu_factors->U = U;
	lu_factors->P = P;
	lu_factors->swaps = swaps;

	return lu_factors;
}


/* Backward subsitution. */
vector_t *back_sub(matrix_t *A, vector_t *B) {

	unsigned int n = A->rows;

	vector_t *X = vctr_empty(n);

	for (unsigned int i = n - 1; i >= 0; --i) {
		X->values[i] = B->values[i];
		for (unsigned int k = i + 1; k < n; ++k)
			X->values[i] -= A->values[i][k] * X->values[k];
		X->values[i] /= A->values[i][i];
	}

	return X;
}


/* Forward substitution. */
vector_t *forward_sub(matrix_t *A, vector_t *B) {

	unsigned int n = A->rows;

	vector_t *X = vctr_empty(n);

	// Forward substitution.
	for (unsigned int i = 0; i < n; ++i) {
		X->values[i] = B->values[i];
		for (unsigned int k = 0; k < i; ++k)
			X->values[i] -= A->values[i][k] * X->values[k];
		X->values[i] /= A->values[i][i];
	}

	return X;
}





/* Generates the transpose of a matrix. */
matrix_t *mtrx_transpose(matrix_t *matrix) {

	matrix_t *transpose = mtrx_zeros(matrix->columns, matrix->rows);

	for (unsigned int i = 0; i < matrix->rows; ++i) {
		for (unsigned int j = 0; j < matrix->columns; ++j)
			transpose->values[j][i] = matrix->values[i][j];
	}

	return transpose;
}


/* Calculates the determinant of the matrix. */
double mtrx_det(matrix_t *A) {

	// Check that the matrix is square.
	if (!mtrx_is_sqr(A))
		return 0;

	lu_factors_t *lu_factors = lu_decomposition(A);

	double det = 1;
	for (unsigned int i = 0; i < lu_factors->L->rows; ++i)
		det *= lu_factors->L->values[i][i];
	if (lu_factors->swaps % 2)
		det *= -1;

	lu_factors_destroy(lu_factors);

	return det;
}


/* Determines the inverse of a matrix, if the matrix is invertible. */
matrix_t *mtrx_inv(matrix_t *A) {

	if (!mtrx_is_sqr(A))
		return NULL;

	matrix_t *inverse = mtrx_zeros(A->rows, A->columns);
	vector_t *id_col = vctr_zeros(A->rows);
	id_col->values[0] = 1;

	// Use LU decomposition to decompose A into an upper and lower triangular matrix.
	lu_factors_t *lu_factors = lu_decomposition(A);

	// Solve the system.
	vector_t *C = mtrx_mult_vctr(lu_factors->P, id_col);
	vector_t *D = forward_sub(lu_factors->L, C);
	vector_t *X = back_sub(lu_factors->U, D);
	mtrx_set_col(inverse, X, 0);


	for (unsigned int i = 1; i < id_col->length; ++i) {
		id_col->values[i - 1] = 0;
		id_col->values[i] = 1;

		C = mtrx_mult_vctr(lu_factors->P, id_col);
		D = forward_sub(lu_factors->L, C);
		X = back_sub(lu_factors->U, D);

		mtrx_set_col(inverse, X, i);
	}

	// Release memory from temporary objects.
	vctr_destroy(C);
	vctr_destroy(D);
	vctr_destroy(X);
	lu_factors_destroy(lu_factors);

	return inverse;
}


/* 
 * Solves a linear system of equations defined by [A][x]=[B],
 * where [A] is a matrix, and [x] and [B] are vectors. 
 */
vector_t *mtrx_solve(matrix_t *A, vector_t *B) {

	// Use LU decomposition to decompose A into an upper and lower triangular matrix.
	lu_factors_t *lu_factors = lu_decomposition(A);

	// Solve the system.
	vector_t *C = mtrx_mult_vctr(lu_factors->P, B);
	vector_t *D = forward_sub(lu_factors->L, C);
	vector_t *X = back_sub(lu_factors->U, D);

	// Release memory from temporary objects.
	vctr_destroy(C);
	vctr_destroy(D);
	lu_factors_destroy(lu_factors);

	return X;
}


/*
 * Checks if a matrix is diagonally dominant. Returns true if the matrix is
 * diagonally dominant, false otherwise.
 */
int mtrx_is_diag_dom(matrix_t *matrix) {

	// Iterate through each row of the matrix and check if that row fits the
	// criteria for diagonal dominance.
	for (U32 i = 0; i < matrix->rows; ++i) {
		double sum = 0;

		// Sum the values of the elements in the row that are not on the
		// diagonal.
		for (U32 j = 0; j < matrix->columns; ++j) {
			if (j != i)
				sum += fabs(matrix->values[i][j]);
		}

		// Absolute value of diagonal entry was not greater than or equal to 
		// the sum of absolute values of the other entries. Matrix is not
		// diagonally dominant.
		if (fabs(matrix->values[i][i]) < sum)
			return 0;
	}

	return 0;
}


int mtrx_make_diag_dom(matrix_t *matrix) {

	if (!mtrx_is_sqr(matrix))
		return 0;

	U32 *mark = (U32 *)malloc(matrix->rows * sizeof(U32));

	for (U32 i = 0; i < matrix->rows; ++i) {
		U32 largest_col = 0;
		double sum = 0;
		for (U32 j = 1; j < matrix->columns; ++j) {
			if (fabs(matrix->values[i][j]) > fabs(matrix->values[i][largest_col]))
				largest_col = j;
			sum += fabs(matrix->values[i][j]);
		}

		// Remove largest absolute value from the sum.
		sum -= fabs(matrix->values[i][largest_col]);

		// Check if this row's largest absolute value was not greater than or equal to
		// the sum of the absolute values of the other elements. If so, the matrix cannot
		// be rearranged into a diagonally dominant form.
		if (fabs(matrix->values[i][largest_col]) < sum)
			return 0;

		// Another row already has its largest value in this column, so the
		// matrix cannot be made diagonally dominant.
		if (mark[largest_col] != -1)
			return 0;
		mark[largest_col] = i;
	}

	// Make the matrix diagonally dominant by swapping rows.
	for (U32 i = 0; i < matrix->rows; ++i) {
		if (mark[i] != i)
			mtrx_row_swap(matrix, i, mark[i]);
	}

	// Matrix was successfully made diagonally dominant.
	return 1;
}


/*
 * Solve a linear system using the Gauss-Siedel iterative technique.
 * The A matrix must be diagonally dominanant to ensure convergence.
 */
void mtrx_solve_gs(matrix_t *A, vector_t *B, vector_t *X, double tolerance) {

	unsigned int max_iter = 1000000;
	unsigned int iter = 0;

	int done = 0;

	while (!done && iter < max_iter) {
		done = 1;
		for (unsigned int i = 0; i < X->length; ++i) {
			double prev = X->values[i];
			X->values[i] = B->values[i];
			for (unsigned int j = 0; j < i; ++j)
				X->values[i] -= A->values[i][j] * X->values[j];
			for (unsigned int j = i + 1; j < X->length; ++j)
				X->values[i] -= A->values[i][j] * X->values[j];
			X->values[i] /= A->values[i][i];

			// Check if this element satisfies the tolerance requirement.
			if (fabs(X->values[i] - prev) > tolerance)
				done = 0;
		}
		++iter;
	}
}


int main() {

	matrix_t *A = mtrx_rnd(5, 7, 10);
	matrix_t *B = mtrx_rnd(7, 5, 10);

	srand(time(0));

	//time_t now = time(0);
	matrix_t *C = mtrx_mult(A, B);
	mtrx_print(C);
	//printf("Naive: %i\n", time(0) - now);
	//now = time(0);
	matrix_t *D = mtrx_strassen_mult(A, B);
	printf("\n");
	mtrx_print(D);
	//printf("'Fast': %i\n", time(0) - now);

	mtrx_destroy(A);
	mtrx_destroy(B);
	mtrx_destroy(C);
	mtrx_destroy(D);

	system("PAUSE");
	return 0;
}