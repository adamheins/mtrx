/*****************************************************************************
 * Copyright (c) 2015 Adam Heins                                             *
 *                                                                           *
 * This file is part of the mtrx project, which is distributed under the MIT *
 * license. For the full terms, see the included LICENSE file.               *
 *****************************************************************************/

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
#include "mtrx.h"


/*===========================================================================*
 *                       HELPER FUNCTIONS AND STRUCTS                        *
 *===========================================================================*/

// Stores the results of an LU decomposition.
// L is the lower triangular, U is the upper triangular, and P is the
// permutation matrix. swaps is the number of row swaps that were performed
// by pivoting.
typedef struct {
  matrix_t *L;
  matrix_t *U;
  matrix_t *P;
  uint32_t swaps;
} lu_factors_t;


// Frees the memory associated with an lu_factors_t object.
void lu_factors_destroy(lu_factors_t *lu_factors) {
  mtrx_destroy(lu_factors->L);
  mtrx_destroy(lu_factors->U);
  mtrx_destroy(lu_factors->P);
  free(lu_factors);
}


// Splits a block from an existing matrix, and returns it as a new matrix. Also
// pads the new block matrix to a specified number of rows and columns.
matrix_t *mtrx_padded_split(matrix_t *matrix, size_t r1, size_t r2, size_t c1,
    size_t c2, size_t rows, size_t cols) {

  size_t num_rows = r2 - r1;
  size_t num_cols = c2 - c1;
  matrix_t *child_matrix = mtrx_zeros(rows, cols);

  for (size_t i = 0; i < num_rows; ++i) {
    for (size_t j = 0; j < num_cols; ++j)
      child_matrix->values[i][j] = matrix->values[r1 + i][c1 + j];
  }

  return child_matrix;
}


// Performs pivoting to ensure a stable answer for LU decomposition.
bool lu_decomposition_pivot(matrix_t *L, matrix_t *P, size_t col) {
  size_t pivot = col;

  // Find largest absolute value in the column. That will be the pivot.
  for (size_t i = col + 1; i < L->rows; ++i) {
    if (fabs(L->values[i][col]) > fabs(L->values[pivot][col]))
      pivot = i;
  }

  // Swap rows such that the pivot is in the current position of [col, col].
  if (pivot != col) {
    mtrx_row_swap(L, col, pivot);
    mtrx_row_swap(P, col, pivot);
    return true;
  }
  return false;
}


// Performs LU decomposition on a matrix, M, producing a lower- and an upper-
// triangular matrix, L and U, such that [L][U] = [M].
lu_factors_t *lu_decomposition(matrix_t *matrix) {
  size_t n = matrix->rows;
  size_t swaps = 0;

  matrix_t *L = mtrx_copy(matrix);
  matrix_t *U = mtrx_id(n);
  matrix_t *P = mtrx_id(n);

  swaps += lu_decomposition_pivot(L, P, 0);

  // First row of U.
  for (size_t j = 1; j < n; ++j)
    U->values[0][j] = L->values[0][j] / L->values[0][0];

  for (size_t j = 1; j < n - 1; ++j) {

    // Compute the next column of L.
    for (size_t i = j; i < n; ++i) {
      for (size_t k = 0; k < j; ++k)
        L->values[i][j] -= L->values[i][k] * U->values[k][j];
    }

    // Pivoting.
    swaps += lu_decomposition_pivot(L, P, j);

    // Compute the next row of U.
    for (size_t k = j + 1; k < n; ++k) {
      U->values[j][k] = L->values[j][k];
      for (size_t i = 0; i < j; ++i)
        U->values[j][k] -= L->values[j][i] * U->values[i][k];
      U->values[j][k] /= L->values[j][j];
    }
  }

  for (size_t k = 0; k < n - 1; ++k)
    L->values[n - 1][n - 1] -= L->values[n - 1][k] * U->values[k][n - 1];

  // Zero the upper triangle of L to make it a lower triangular matrix.
  for (size_t i = 0; i < n - 1; ++i) {
    for (size_t j = i + 1; j < n; j++)
      L->values[i][j] = 0;
  }

  lu_factors_t *lu_factors = (lu_factors_t *)malloc(sizeof(lu_factors_t));
  lu_factors->L = L;
  lu_factors->U = U;
  lu_factors->P = P;
  lu_factors->swaps = swaps;

  return lu_factors;
}


// Perform the naive matrix multiplication algorithm.
matrix_t *mtrx_naive_mult(matrix_t *A, matrix_t *B) {
  matrix_t *C = mtrx_zeros(A->rows, B->columns);

  for (size_t i = 0; i < B->columns; ++i) {
    for (size_t j = 0; j < A->rows; ++j) {
      for (size_t k = 0; k < A->columns; ++k)
        C->values[j][i] += A->values[j][k] * B->values[k][i];
    }
  }
  return C;
}


// Perform the 'fast' Strassen matrix multiplication algorithm. Breaks down the
// matrices recursively and only uses seven multiplications rather than the
// usual eight used by the naive algorithm. Faster for larger matrices.
matrix_t *mtrx_strassen_mult(matrix_t *A, matrix_t *B) {

  // Ending condition.
  if (A->rows <= 128 || A->columns <= 128)
    return mtrx_naive_mult(A, B);

  size_t Ahr = round(A->rows / 2.0);
  size_t Ahc = round(A->columns / 2.0);
  size_t Bhr = round(B->rows / 2.0);
  size_t Bhc = round(B->columns / 2.0);

  // Split A into four submatrices.
  matrix_t *A11 = mtrx_sub_block(A, 0, Ahr, 0, Ahc);
  matrix_t *A12 = mtrx_padded_split(A, 0, Ahr, Ahc, A->columns, Ahr, Ahc);
  matrix_t *A21 = mtrx_padded_split(A, Ahr, A->rows, 0, Ahc, Ahr, Ahc);
  matrix_t *A22 = mtrx_padded_split(A, Ahr, A->rows, Ahc, A->columns, Ahr, Ahc);

  // Split B into four submatrices.
  matrix_t *B11 = mtrx_sub_block(B, 0, Bhr, 0, Bhc);
  matrix_t *B12 = mtrx_padded_split(B, 0, Bhr, Bhc, B->columns, Bhr, Bhc);
  matrix_t *B21 = mtrx_padded_split(B, Bhr, B->rows, 0, Bhc, Bhr, Bhc);
  matrix_t *B22 = mtrx_padded_split(B, Bhr, B->rows, Bhc, B->columns, Bhr, Bhc);

  // Calculate M matrices. This is the recursive section of the algorithm.
  matrix_t *M[7];
  M[0] = mtrx_strassen_mult(mtrx_add(A11, A22), mtrx_add(B11, B22));
  M[1] = mtrx_strassen_mult(mtrx_add(A21, A22), B11);
  M[2] = mtrx_strassen_mult(A11, mtrx_subtract(B12, B22));
  M[3] = mtrx_strassen_mult(A22, mtrx_subtract(B21, B11));
  M[4] = mtrx_strassen_mult(mtrx_add(A11, A12), B22);
  M[5] = mtrx_strassen_mult(mtrx_subtract(A21, A11), mtrx_add(B11, B12));
  M[6] = mtrx_strassen_mult(mtrx_subtract(A12, A22), mtrx_add(B21, B22));

  mtrx_destroy(A11);
  mtrx_destroy(A12);
  mtrx_destroy(A21);
  mtrx_destroy(A22);

  mtrx_destroy(B11);
  mtrx_destroy(B12);
  mtrx_destroy(B21);
  mtrx_destroy(B22);

  // Calculate each quadrant of C from M matrices.
  matrix_t *C11 = mtrx_add(mtrx_subtract(mtrx_add(M[0], M[3]), M[4]), M[6]);
  matrix_t *C12 = mtrx_add(M[2], M[4]);
  matrix_t *C21 = mtrx_add(M[1], M[3]);
  matrix_t *C22 = mtrx_add(mtrx_add(mtrx_subtract(M[0], M[1]), M[2]), M[5]);

  for (uint8_t i = 0; i < 7; ++i)
    mtrx_destroy(M[i]);

  matrix_t *C = mtrx_empty(A->rows, B->columns);

  for (size_t i = 0; i < C11->rows; ++i) {
    for (size_t j = 0; j < C11->columns; ++j)
      C->values[i][j] = C11->values[i][j];
  }

  for (size_t i = 0; i < C11->rows; ++i) {
    for (size_t j = 0; j < C->columns / 2; ++j)
      C->values[i][j + C11->columns] = C12->values[i][j];
  }

  for (size_t i = 0; i < C->rows / 2; ++i) {
    for (size_t j = 0; j < C11->columns; ++j)
      C->values[i + C11->rows][j] = C21->values[i][j];
  }

  for (size_t i = 0; i < C->rows / 2; ++i) {
    for (size_t j = 0; j < C->columns / 2; ++j)
      C->values[i + C11->rows][j + C11->columns] = C22->values[i][j];
  }

  mtrx_destroy(C11);
  mtrx_destroy(C12);
  mtrx_destroy(C21);
  mtrx_destroy(C22);

  return C;
}


// Performs backward substitution to solve a system.
vector_t *back_sub(matrix_t *A, vector_t *B) {
  size_t n = A->rows;

  vector_t *X = vctr_empty(n);

  // For loop syntax adjusted to properly decrement an unsigned value.
  for (size_t i = n; i-- > 0;) {
    X->values[i] = B->values[i];
    for (size_t k = i + 1; k < n; ++k)
      X->values[i] -= A->values[i][k] * X->values[k];
    X->values[i] /= A->values[i][i];
  }
  return X;
}


// Performs forward substitution to solve a system.
vector_t *forward_sub(matrix_t *A, vector_t *B) {
  size_t n = A->rows;

  vector_t *X = vctr_empty(n);

  for (size_t i = 0; i < n; ++i) {
    X->values[i] = B->values[i];
    for (size_t k = 0; k < i; ++k)
      X->values[i] -= A->values[i][k] * X->values[k];
    X->values[i] /= A->values[i][i];
  }
  return X;
}


// Determine the number of columns in the matrix.
size_t get_num_columns(const char *str) {

  size_t count = 0;
  if (str[0] != ' ')
    ++count;

  for (size_t i = 1; str[i] != '\0' && str[i] != ';'; ++i) {
    if (str[i - 1] == ' ' && str[i] != ' ')
      ++count;
  }

  return count;
}


// Parse a row from the matrix.
scalar_t *row_from_str(const char *str, size_t cols, size_t low, size_t high) {
  scalar_t *values = (scalar_t *)malloc(cols * sizeof(scalar_t));

  // Parse each scalar value out of the array.
  size_t count = 0;
  size_t start = 0;
  size_t i = low + 1;

  for (; i < high; ++i) {
    // 'Rising edge'
    if (str[i - 1] == ' ' && str[i] != ' ')
      start = i;

    // 'Falling edge'
    if (str[i - 1] != ' ' && str[i] == ' ')
      values[count++] = atos(str, start, i);
  }

  values[count] = atos(str, start, i);
  return values;
}


/*===========================================================================*
 *                           INDEXER FUNCTIONS                               *
 *===========================================================================*/

/*--------------------------- Initialization --------------------------------*/

indexer_t *indexer_init(size_t length) {
  indexer_t *indexer = (indexer_t *)malloc(sizeof(indexer_t));

  indexer->values = (size_t *)malloc(length * sizeof(size_t));
  indexer->length = length;

  return indexer;
}


/*===========================================================================*
 *                           MATRIX FUNCTIONS                                *
 *===========================================================================*/

/*--------------------------- Initialization --------------------------------*/

matrix_t *mtrx_init(const char *str) {

  // Figure out how many rows the matrix will have.
  size_t row_count = 1;
  for (size_t i = 0; str[i] != '\0'; ++i) {
    if (str[i] == ';')
      ++row_count;
  }

  // Create the matrix.
  matrix_t *matrix = (matrix_t *)malloc(sizeof(matrix_t));
  matrix->rows = row_count;
  matrix->columns = get_num_columns(str);
  matrix->values = (scalar_t **)malloc(row_count * sizeof(scalar_t *));

  size_t start = 0;
  size_t index = 0;
  size_t count = 0;
  for (; str[index] != '\0'; ++index) {
    if (str[index] == ';') {
      matrix->values[count++] = row_from_str(str, matrix->columns, start, index);
      start = index + 1;
    }
  }

  matrix->values[count] = row_from_str(str, matrix->columns, start, index);
  return matrix;
}


matrix_t *mtrx_empty(size_t rows, size_t cols) {
  matrix_t *matrix = (matrix_t *)malloc(sizeof(matrix_t));

  // Allocate the array of rows.
  matrix->values = (scalar_t **)malloc(rows * sizeof(scalar_t *));

  // Allocate the columns arrays.
  for (size_t i = 0; i < rows; ++i)
    matrix->values[i] = (scalar_t *)malloc(cols * sizeof(scalar_t));

  matrix->rows = rows;
  matrix->columns = cols;

  return matrix;
}


matrix_t *mtrx_zeros(size_t rows, size_t cols) {
  matrix_t *matrix = (matrix_t *)malloc(sizeof(matrix_t));

  // Allocate the array of rows.
  matrix->values = (scalar_t **)calloc(rows, sizeof(scalar_t *));

  // Allocate the columns arrays.
  for (size_t i = 0; i < rows; ++i)
    matrix->values[i] = (scalar_t *)calloc(cols, sizeof(scalar_t));

  matrix->rows = rows;
  matrix->columns = cols;

  return matrix;
}


matrix_t *mtrx_ones(size_t rows, size_t cols) {
  matrix_t *matrix = (matrix_t *)malloc(sizeof(matrix_t));

  // Allocate the array of rows.
  matrix->values = (scalar_t **)malloc(rows * sizeof(scalar_t *));

  // Allocate the columns arrays.
  for (size_t i = 0; i < rows; ++i)
    matrix->values[i] = (scalar_t *)malloc(cols * sizeof(scalar_t));

  matrix->rows = rows;
  matrix->columns = cols;

  for (size_t r = 0; r < rows; ++r) {
    for (size_t c = 0; c < cols; ++c)
      matrix->values[r][c] = 1;
  }

  return matrix;
}


matrix_t *mtrx_id(size_t n) {
  matrix_t *id = mtrx_zeros(n, n);

  // Fill diagonal with zeros.
  for (size_t i = 0; i < n; ++i)
    id->values[i][i] = 1;

  return id;
}


matrix_t *mtrx_diag(vector_t *vector) {
  matrix_t *matrix = mtrx_zeros(vector->length, vector->length);

  for (size_t i = 0; i < vector->length; ++i)
    matrix->values[i][i] = vector->values[i];

  return matrix;
}


matrix_t *mtrx_rnd(size_t rows, size_t cols, uint32_t max) {
  matrix_t *matrix = mtrx_empty(rows, cols);

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j)
      matrix->values[i][j] = rand() % max;
  }

  return matrix;
}


matrix_t *mtrx_copy(matrix_t *matrix) {
  matrix_t *copy = mtrx_empty(matrix->rows, matrix->columns);

  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j)
      copy->values[i][j] = matrix->values[i][j];
  }

  return copy;
}


/*---------------------------- Deallocation ---------------------------------*/

void mtrx_destroy(matrix_t *matrix) {
  for (size_t i = 0; i < matrix->rows; ++i)
    free(matrix->values[i]);
  free(matrix->values);
  free(matrix);
  matrix = NULL;
}


/*-------------------------------- Display ----------------------------------*/

void mtrx_print(matrix_t *matrix) {
  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j)
      printf("%f ", matrix->values[i][j]);
    printf("\n");
  }
}


/*------------------------------ Comparisons --------------------------------*/

bool mtrx_is_sqr(matrix_t *matrix) {
  return matrix->rows == matrix->columns;
}


bool mtrx_eq_dim(matrix_t *A, matrix_t *B) {
  return (A->rows == B->rows && A->columns == B->columns);
}


bool mtrx_eq(matrix_t *A, matrix_t *B) {
  if (!mtrx_eq_dim(A, B))
    return false;
  for (size_t i = 0; i < A->rows; ++i) {
    for (size_t j = 0; j < A->columns; ++j) {
      if (A->values[i][j] != B->values[i][j])
        return false;
    }
  }
  return true;
}


/*------------------------------ Max and Min --------------------------------*/

scalar_t mtrx_max(matrix_t *matrix) {
  scalar_t largest = matrix->values[0][0];

  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j) {
      if (matrix->values[i][j] > largest)
        largest = matrix->values[i][j];
    }
  }
  return largest;
}


scalar_t mtrx_min(matrix_t *matrix) {
  scalar_t smallest = matrix->values[0][0];

  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j) {
      if (matrix->values[i][j] < smallest)
        smallest = matrix->values[i][j];
    }
  }
  return smallest;
}


/*--------------------------- Matrix Arithmetic -----------------------------*/

matrix_t *mtrx_add(matrix_t *A, matrix_t *B) {

  // Check that matrices have the same dimensions.
  if (!mtrx_eq_dim(A, B))
    return NULL;

  // Create a new matrix to hold the result.
  matrix_t *C = mtrx_empty(A->rows, A->columns);

  for (size_t i = 0; i < A->rows; ++i) {
    for (size_t j = 0; j < A->columns; j++)
      C->values[i][j] = A->values[i][j] + B->values[i][j];
  }

  return C;
}


matrix_t *mtrx_subtract(matrix_t *A, matrix_t *B) {

  // Check that matrices have the same dimensions.
  if (!mtrx_eq_dim(A, B))
    return NULL;

  // Create a new matrix to hold the result.
  matrix_t *C = mtrx_empty(A->rows, A->columns);

  // Subtract each element in B from the corresponding element in A to
  // generate the resulting element in C.
  for (size_t i = 0; i < A->rows; ++i) {
    for (size_t j = 0; j < A->columns; ++j)
      C->values[i][j] = A->values[i][j] - B->values[i][j];
  }

  return C;
}


matrix_t *mtrx_scale(matrix_t *matrix, scalar_t scalar) {
  matrix_t *multiple = mtrx_empty(matrix->rows, matrix->columns);

  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j)
      multiple->values[i][j] = matrix->values[i][j] * scalar;
  }

  return multiple;
}


matrix_t *mtrx_mult(matrix_t *A, matrix_t *B) {

  // Check for matrices inappropriately sized for multiplication.
  if (A->columns != B->rows)
    return NULL;

  return mtrx_strassen_mult(A, B);
}


vector_t *mtrx_mult_vctr(matrix_t *A, vector_t *B) {
  if (A->columns != B->length)
    return NULL;

  vector_t *C = vctr_zeros(A->rows);

  for (size_t i = 0; i < A->rows; ++i) {
    for (size_t j = 0; j < A->columns; ++j)
      C->values[i] += A->values[i][j] * B->values[j];
  }

  return C;
}


/*------------------------- Pointwise Arithmetic ----------------------------*/

matrix_t *mtrx_pw_mult(matrix_t *A, matrix_t *B) {
  if (!mtrx_eq_dim(A, B))
    return NULL;

  matrix_t *C = mtrx_empty(A->rows, A->columns);

  for (size_t r = 0; r < A->rows; ++r) {
    for (size_t c = 0; c < A->columns; ++c)
      C->values[r][c] = A->values[r][c] * B->values[r][c];
  }

  return C;
}


matrix_t *mtrx_pw_div(matrix_t *A, matrix_t *B) {
  if (!mtrx_eq_dim(A, B))
    return NULL;

  matrix_t *C = mtrx_empty(A->rows, A->columns);

  for (size_t r = 0; r < A->rows; ++r) {
    for (size_t c = 0; c < A->columns; ++c)
      C->values[r][c] = A->values[r][c] / B->values[r][c];
  }

  return C;
}


matrix_t *mtrx_pw_pow(matrix_t *A, matrix_t *B) {
  if (!mtrx_eq_dim(A, B))
    return NULL;

  matrix_t *C = mtrx_empty(A->rows, A->columns);

  for (size_t r = 0; r < A->rows; ++r) {
    for (size_t c = 0; c < A->columns; ++c)
      C->values[r][c] = pow(A->values[r][c], B->values[r][c]);
  }

  return C;
}


/*-------------------------- Matrix Manipulation ----------------------------*/

void mtrx_row_swap(matrix_t *matrix, size_t r1, size_t r2) {
  scalar_t temp;

  for (size_t i = 0; i < matrix->columns; ++i) {
    temp = matrix->values[r1][i];
    matrix->values[r1][i] = matrix->values[r2][i];
    matrix->values[r2][i] = temp;
  }
}


void mtrx_col_swap(matrix_t *matrix, size_t c1, size_t c2) {

  scalar_t temp;

  for (size_t i = 0; i < matrix->rows; ++i) {
    temp = matrix->values[i][c1];
    matrix->values[i][c1] = matrix->values[i][c2];
    matrix->values[i][c2] = temp;
  }
}


void mtrx_scale_row(matrix_t *matrix, size_t row, scalar_t scalar) {
  for (size_t i = 0; i < matrix->columns; ++i)
    matrix->values[row][i] *= scalar;
}


void mtrx_scale_col(matrix_t *matrix, size_t col, scalar_t scalar) {
  for (size_t i = 0; i < matrix->rows; ++i)
    matrix->values[i][col] *= scalar;
}


/*----------------- Row and Column Accessors and Mutators -------------------*/

vector_t *mtrx_get_row(matrix_t *matrix, size_t row) {
  vector_t *vector = vctr_empty(matrix->columns);

  for (size_t i = 0; i < vector->length; ++i)
    vector->values[i] = matrix->values[row][i];

  return vector;
}


vector_t *mtrx_get_col(matrix_t *matrix, size_t col) {
  vector_t *vector = vctr_empty(matrix->rows);

  for (size_t i = 0; i < vector->length; ++i)
    vector->values[i] = matrix->values[i][col];

  return vector;
}


void mtrx_set_row(matrix_t *matrix, vector_t *vector, size_t row) {
  if (matrix->columns != vector->length)
    return;

  // Replace values in matrix row with those of the vector.
  for (size_t i = 0; i < vector->length; ++i)
    matrix->values[i][row] = vector->values[i];
}


void mtrx_set_col(matrix_t *matrix, vector_t *vector, size_t col) {
  if (matrix->rows != vector->length)
    return;

  // Replace values in matrix column with those of the vector.
  for (size_t i = 0; i < vector->length; ++i)
    matrix->values[i][col] = vector->values[i];
}


/*----------------------------- Sub-matrices --------------------------------*/

matrix_t *mtrx_sub_matrix(matrix_t *A, indexer_t *rows, indexer_t *columns) {
  matrix_t *sub = mtrx_zeros(rows->length, columns->length);

  for (size_t i = 0; i < rows->length; ++i) {
    for (size_t j = 0; j < columns->length; ++j)
      sub->values[i][j] = A->values[rows->values[i]][columns->values[j]];
  }
  return sub;
}


matrix_t *mtrx_sub_block(matrix_t *matrix, size_t r1, size_t r2, size_t c1,
    size_t c2) {

  size_t num_rows = r2 - r1;
  size_t num_cols = c2 - c1;
  matrix_t *child_matrix = mtrx_zeros(num_rows, num_cols);

  for (size_t i = 0; i < num_rows; ++i) {
    for (size_t j = 0; j < num_cols; ++j)
      child_matrix->values[i][j] = matrix->values[r1 + i][c1 + j];
  }
  return child_matrix;
}


/*------------------------------- Operations --------------------------------*/

matrix_t *mtrx_transpose(matrix_t *matrix) {
  matrix_t *transpose = mtrx_zeros(matrix->columns, matrix->rows);

  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j)
      transpose->values[j][i] = matrix->values[i][j];
  }
  return transpose;
}


scalar_t mtrx_det(matrix_t *matrix) {

  // Check that the matrix is square.
  if (!mtrx_is_sqr(matrix))
    return 0;

  lu_factors_t *lu_factors = lu_decomposition(matrix);

  scalar_t det = 1;
  for (size_t i = 0; i < lu_factors->L->rows; ++i)
    det *= lu_factors->L->values[i][i];

  if (lu_factors->swaps % 2)
    det *= -1;

  lu_factors_destroy(lu_factors);

  return det;
}


matrix_t *mtrx_inv(matrix_t *matrix) {
  if (!mtrx_is_sqr(matrix))
    return NULL;

  matrix_t *inverse = mtrx_zeros(matrix->rows, matrix->columns);
  vector_t *id_col = vctr_zeros(matrix->rows);
  id_col->values[0] = 1;

  // Use LU decomposition to decompose A into an upper and lower triangular matrix.
  lu_factors_t *lu_factors = lu_decomposition(matrix);

  // Solve the system.
  vector_t *C = mtrx_mult_vctr(lu_factors->P, id_col);
  vector_t *D = forward_sub(lu_factors->L, C);
  vector_t *X = back_sub(lu_factors->U, D);
  mtrx_set_col(inverse, X, 0);

  for (size_t i = 1; i < id_col->length; ++i) {
    id_col->values[i - 1] = 0;
    id_col->values[i] = 1;

    C = mtrx_mult_vctr(lu_factors->P, id_col);
    D = forward_sub(lu_factors->L, C);
    X = back_sub(lu_factors->U, D);
    mtrx_set_col(inverse, X, i);

    vctr_destroy(C);
    vctr_destroy(D);
    vctr_destroy(X);
  }

  vctr_destroy(id_col);
  lu_factors_destroy(lu_factors);

  return inverse;
}

/*----------------------------- System Solving ------------------------------*/

vector_t *mtrx_solve(matrix_t *A, vector_t *B) {

  // Use LU decomposition to decompose A into an upper and lower triangular
  // matrix.
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


bool mtrx_is_diag_dom(matrix_t *matrix) {

  // Iterate through each row of the matrix and check if that row fits the
  // criteria for diagonal dominance.
  for (size_t i = 0; i < matrix->rows; ++i) {
    scalar_t sum = 0;

    // Sum the values of the elements in the row that are not on the
    // diagonal.
    for (size_t j = 0; j < matrix->columns; ++j) {
      if (j != i)
        sum += fabs(matrix->values[i][j]);
    }

    // Absolute value of diagonal entry was not greater than or equal to
    // the sum of absolute values of the other entries. Matrix is not
    // diagonally dominant.
    if (fabs(matrix->values[i][i]) < sum)
      return false;
  }

  return true;
}


matrix_t *mtrx_make_diag_dom(matrix_t *matrix) {
  if (!mtrx_is_sqr(matrix))
    return false;

  int32_t *mark = (int32_t *)malloc(matrix->rows * sizeof(int32_t));
  for (size_t i = 0; i < matrix->rows; ++i)
    mark[i] = -1;

  for (size_t i = 0; i < matrix->rows; ++i) {
    size_t largest_col = 0;
    scalar_t sum = 0;


    for (size_t j = 1; j < matrix->columns; ++j) {
      if (fabs(matrix->values[i][j]) > fabs(matrix->values[i][largest_col]))
        largest_col = j;
      sum += fabs(matrix->values[i][j]);
    }

    // Remove largest absolute value from the sum.
    sum -= fabs(matrix->values[i][largest_col]);

    // Check if this row's largest absolute value was not greater than or equal
    // to the sum of the absolute values of the other elements. If so, the
    // matrix cannot be rearranged into a diagonally dominant form.
    if (fabs(matrix->values[i][largest_col]) < sum) {
      free(mark);
      return NULL;
    }

    // Another row already has its largest value in this column, so the
    // matrix cannot be made diagonally dominant.
    if (mark[largest_col] != -1) {
      free(mark);
      return NULL;
    }
    mark[largest_col] = i;
  }

  // Create a the rearranged diagonally dominant matrix.
  matrix_t *diag_dom_matrix = mtrx_empty(matrix->rows, matrix->columns);
  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j)
      diag_dom_matrix->values[i][j] = matrix->values[mark[i]][j];
  }

  free(mark);

  return diag_dom_matrix;
}


void mtrx_solve_gs(matrix_t *A, vector_t *B, vector_t *X, double tolerance) {
  uint32_t max_iter = 1000000;
  uint32_t iter = 0;

  bool done = false;

  while (!done && iter < max_iter) {
    done = true;
    for (size_t i = 0; i < X->length; ++i) {
      scalar_t prev = X->values[i];
      X->values[i] = B->values[i];
      for (size_t j = 0; j < i; ++j)
        X->values[i] -= A->values[i][j] * X->values[j];
      for (size_t j = i + 1; j < X->length; ++j)
        X->values[i] -= A->values[i][j] * X->values[j];
      X->values[i] /= A->values[i][i];

      // Check if this element satisfies the tolerance requirement.
      if (fabs(X->values[i] - prev) > tolerance)
        done = false;
    }
    ++iter;
  }
}

