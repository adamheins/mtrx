#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"


#define NULL (void *)0;

typedef struct {
	int *values;
	unsigned int length;
} array_t;


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

	matrix_t *copy = (matrix_t *)malloc(sizeof(matrix_t));
	matrix_init(copy, matrix->rows, matrix->columns);

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


matrix_t *mtrx_naive_mult(matrix_t *A, matrix_t *B) {

	// Create a new matrix to hold the result.
	matrix_t *C = (matrix_t *)malloc(sizeof(matrix_t));
	matrix_init(C, A->rows, B->columns);
	
	for (int i = 0; i < B->columns; ++i) {
		for (int j = 0; j < A->rows; ++j) {
			for (int k = 0; k < A->columns; ++k) {
				C->values[j][i] += A->values[j][k] * B->values[k][i];
			}
		}
	}
	return C;
}


matrix_t *fast_matrix_multiply(matrix_t *A, matrix_t *B) {
	//TODO
}


matrix_t *mtrx_mult(matrix_t *A, matrix_t *B) {
	
	// Check for matrices inappropriately sized for multiplication.
	if (A->columns != B->rows)
		return NULL;

	// Currently defaulting to naive multiplication.
	return mtrx_naive_mult(A, B);
}


matrix_t *mtrx_scale(matrix_t *matrix, double scalar) {

	matrix_t *multiple = (matrix_t *)malloc(sizeof(matrix_t));
	matrix_init(multiple, matrix->rows, matrix->columns);

	for (unsigned int i = 0; i < matrix->rows; ++i) {
		for (unsigned int j = 0; j < matrix->columns; ++j)
			multiple->values[i][j] = matrix->values[i][j] * scalar;
	}

	return multiple;
}

/* Returns true if the two matrices have equal dimensions, false otherwise. */
int mtrx_eq_dim(matrix_t *A, matrix_t *B) {
	return (A->rows == B->rows && A->columns == B->columns);
}


/*
* Pointwise multiplication of a matrix.
* Multiplies each element of A by the corresponding element of B.
*/
matrix_t *mtrx_pw_mult(matrix_t *A, matrix_t *B) {

	// Check that matrices have the same dimensions.
	if (!mtrx_eq_dim(A, B))
		return NULL;

	matrix_t *C = (matrix_t *)malloc(sizeof(matrix_t));
	matrix_init(C, A->rows, A->columns);

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

	matrix_t *C = (matrix_t *)malloc(sizeof(matrix_t));
	matrix_init(C, A->rows, A->columns);

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

	matrix_t *C = (matrix_t *)malloc(sizeof(matrix_t));
	matrix_init(C, A->rows, A->columns);

	for (unsigned int r = 0; r < A->rows; ++r) {
		for (unsigned int c = 0; c < A->columns; ++c)
			C->values[r][c] = pow(A->values[r][c], B->values[r][c]);
	}

	return C;
}


/*
* Add two matrices.
*/
matrix_t *mtrx_add(matrix_t *A, matrix_t *B) {

	// Check that matrices have the same dimensions.
	if (!mtrx_eq_dim(A, B))
		return NULL;

	// Create a new matrix to hold the result.
	matrix_t *C = (matrix_t *)malloc(sizeof(matrix_t));
	matrix_init(C, A->rows, B->columns);

	for (int i = 0; i < A->rows; ++i) {
		for (int j = 0; j < A->columns; j++)
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
	matrix_t *C = (matrix_t *)malloc(sizeof(matrix_t));
	matrix_init(C, A->rows, B->columns);

	// Subtract each element in B from the corresponding element in A to
	// generate the resulting element in C.
	for (int i = 0; i < A->rows; ++i) {
		for (int j = 0; j < A->columns; j++)
			C->values[i][j] = A->values[i][j] - B->values[i][j];
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


/* Generates the transpose of a matrix. */
matrix_t *mtrx_transpose(matrix_t *matrix) {

	matrix_t *transpose = mtrx_zeros(matrix->columns, matrix->rows);

	for (unsigned int i = 0; i < matrix->rows; ++i) {
		for (unsigned int j = 0; j < matrix->columns; ++j)
			transpose->values[j][i] = matrix->values[i][j];
	}

	return transpose;
}


/* Performs pivoting to ensure a stable answer for LU decomposition. */
void lu_decomposition_pivot(matrix_t *L, matrix_t *P, unsigned int col) {

	unsigned int pivot = col;

	// Find largest absolute value in the column. That will be the pivot.
	for (unsigned int i = col + 1; i < L->rows; ++i) {
		if (abs(L->values[i][col]) > abs(L->values[pivot][col]))
			pivot = i;
	}

	// Swap rows such that the pivot is in the current position of [col, col].
	if (pivot != col) {
		mtrx_row_swap(L, col, pivot);
		mtrx_row_swap(P, col, pivot);
	}
}


// A is nxn and B is nx1
lu_factors_t *lu_decomposition(matrix_t *A) {

	unsigned int n = A->rows;
	unsigned int swaps = 0;

	matrix_t *L = copy_matrix(A);
	matrix_t *U = mtrx_id(n, n);
	matrix_t *P = mtrx_id(n, n);

	lu_decomposition_pivot(L, P, 0);

	// First row of U.
	for (int j = 1; j < n; ++j)
		U->values[0][j] = L->values[0][j] / L->values[0][0];

	for (int j = 1; j < n - 1; ++j) {

		// Compute the next column of L.
		for (int i = j; i < n; ++i) {
			for (int k = 0; k < j; ++k)
				L->values[i][j] -= L->values[i][k] * U->values[k][j];
		}

		// Pivoting.
		lu_decomposition_pivot(L, P, j);

		// Compute the next row of U.
		for (int k = j + 1; k < n; ++k) {
			U->values[j][k] = L->values[j][k];
			for (int i = 0; i < j; ++i)
				U->values[j][k] -= L->values[j][i] * U->values[i][k];
			U->values[j][k] /= L->values[j][j];
		}
	}

	for (int k = 0; k < n - 1; ++k)
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

	// Testing.
	/*
	printf("L=\n");
	matrix_print(L);
	printf("\nU=\n");
	matrix_print(U);
	printf("\nP=\n");
	matrix_print(P);
	printf("\nA=\n");
	matrix_print(A);
	printf("\nprod=\n");
	matrix_print(matrix_multiply(matrix_multiply(mtrx_transpose(P), L), U));
	printf("\n");
	*/

	return lu_factors;
}


vector_t *back_sub(matrix_t *A, vector_t *B) {

	unsigned int n = A->rows;

	vector_t *X = vctr_empty(n);

	for (int i = n - 1; i >= 0; --i) {
		X->values[i] = B->values[i];
		for (int k = i + 1; k < n; ++k)
			X->values[i] -= A->values[i][k] * X->values[k];
		X->values[i] /= A->values[i][i];
	}

	return X;
}


vector_t *forward_sub(matrix_t *A, vector_t *B) {

	unsigned int n = A->rows;

	vector_t *X = vctr_empty(n);

	// Forward substitution.
	for (int i = 0; i < n; ++i) {
		X->values[i] = B->values[i];
		for (int k = 0; k < i; ++k)
			X->values[i] -= A->values[i][k] * X->values[k];
		X->values[i] /= A->values[i][i];
	}

	return X;
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


vector_t *vctr_cross_prod(vector_t *A, vector_t *B) {
	//TODO
}


/* Calculates the magnitude of a vector. */
double vctr_mag(vector_t *vector) {

	double mag = 0;

	for (unsigned int i = 0; i < vector->length; ++i)
		mag += vector->values[i] * vector->values[i];

	return sqrt(mag);
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


double mtrx_det(matrix_t *A) {

	// Check that the matrix is square.
	if (!mtrx_is_sqr(A))
		return 0;

	lu_factors_t *lu_factors = lu_decomposition(A);

	double det = 1;
	for (unsigned int i = 0; i < lu_factors->L->rows; ++i)
		det *= lu_factors->L->values[i][i];

	lu_factors_destroy(lu_factors);
		
	return det;
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


/* Creates a submatrix from the given matrix. */
matrix_t *mtrx_submatrix(matrix_t *A, array_t *rows, array_t *columns) {

	matrix_t *sub = (matrix_t *)malloc(sizeof(matrix_t));
	matrix_init(sub, rows->length, columns->length);

	for (int i = 0; i < rows->length; ++i) {
		for (int j = 0; j < columns->length; ++j)
			sub->values[i][j] = A->values[rows->values[i]][columns->values[j]];
	}

	return sub;
}







/*
Helpful built ins:

det()
cross product
gauss-siedel

*/

int main() {

	matrix_t *A = mtrx_zeros(4, 4);
	vector_t *B = vctr_empty(3);

	/*
	A->values[0][0] = 6;
	A->values[0][1] = -4;
	A->values[0][2] = 1;

	A->values[1][0] = -4;
	A->values[1][1] = 6;
	A->values[1][2] = -4;

	A->values[2][0] = 1;
	A->values[2][1] = -4;
	A->values[2][2] = 6;
	*/
	
	B->values[0] = -14;
	B->values[1] = 36;
	B->values[2] = 6;

	/*
	A->values[0][0] = -2;
	A->values[0][1] = 2;
	A->values[0][2] = -3;

	A->values[1][0] = -1;
	A->values[1][1] = 1;
	A->values[1][2] = 3;

	A->values[2][0] = 2;
	A->values[2][1] = 0;
	A->values[2][2] = -1;*/
	
	A->values[0][0] = 3;
	A->values[0][1] = -1;
	A->values[0][2] = 4;
	A->values[0][3] = 5;

	A->values[1][0] = -1;
	A->values[1][1] = 2;
	A->values[1][2] = 7;
	A->values[1][3] = 8;

	A->values[2][0] = -4;
	A->values[2][1] = 3;
	A->values[2][2] = 5;
	A->values[2][3] = -1;

	A->values[3][0] = 6;
	A->values[3][1] = -2;
	A->values[3][2] = -3;
	A->values[3][3] = 4;

	matrix_t *inv = mtrx_inv(A);
	matrix_print(inv);
	matrix_print(matrix_multiply(A, inv));
	//printf("%f\n", mtrx_det(A));

	mtrx_destroy(A);
	mtrx_destroy(inv);

	system("PAUSE");
	return 0;
}