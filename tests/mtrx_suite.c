#include "../mtrx.h"
#include "clar.h"



/*--------------------------- Initialization --------------------------------*/

void test_mtrx_suite__mtrx_empty(void) {
  matrix_t *matrix = mtrx_empty(5, 10);
  cl_assert_(matrix->rows == 5, "Rows not assigned properly!");
  cl_assert_(matrix->columns == 10, "Columns not assigned properly.");
  mtrx_destroy(matrix);
}

void test_mtrx_suite__mtrx_zeros(void) {
  matrix_t *matrix = mtrx_zeros(5, 5);
  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j)
      cl_assert_(matrix->values[i][j] == 0, "Values should all be zero.");
  }
  mtrx_destroy(matrix);
}


void test_mtrx_suite__mtrx_ones(void) {
  matrix_t *matrix = mtrx_ones(5, 5);
  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j)
      cl_assert_(matrix->values[i][j] == 1, "Values should all be one.");
  }
  mtrx_destroy(matrix);
}


void test_mtrx_suite__mtrx_id(void) {
  matrix_t *matrix = mtrx_id(10);
  for (size_t i = 0; i < matrix->rows; ++i)
    cl_assert_(matrix->values[i][i] == 1, "Values on diagonal should be one.");
  mtrx_destroy(matrix);
}


void test_mtrx_suite__mtrx_diag(void) {
  vector_t *vector = vctr_empty(10);
  for (size_t i = 0; i < vector->length; ++i)
    vector->values[i] = i;
  matrix_t *matrix = mtrx_diag(vector);

  cl_assert_(matrix->rows == vector->length, "Matrix number of rows should \
      equal vectors length.");
  cl_assert_(matrix->columns == vector->length, "Matrix number of columns \
      should equal vectors length.");

  for (size_t i = 0; i < vector->length; ++i)
    cl_assert_(matrix->values[i][i] == i, "Values on diagonal should match \
        vector.");

  vctr_destroy(vector);
  mtrx_destroy(matrix);

}


void test_mtrx_suite__mtrx_copy(void) {
  matrix_t *matrix = mtrx_rnd(10, 10, 10);
  matrix_t *copy = mtrx_copy(matrix);

  cl_assert_(copy->rows == matrix->rows, "Matrices should have equal number of \
      rows.");
  cl_assert_(copy->columns == matrix->columns, "Matrices should have equal \
      number of columns.");

  cl_assert_(mtrx_eq(matrix, copy), "Values should all be equal.");

  mtrx_destroy(matrix);
  mtrx_destroy(copy);
}


/*------------------------------ Comparisons --------------------------------*/

void test_mtrx_suite__mtrx_is_sqr(void) {
  matrix_t *matrix = mtrx_empty(10, 10);
  cl_assert_(mtrx_is_sqr(matrix), "Matrix is square.");
  mtrx_destroy(matrix);
}


void test_mtrx_suite__mtrx_is_not_sqr(void) {
  matrix_t *matrix = mtrx_empty(10, 11);
  cl_assert_(!mtrx_is_sqr(matrix), "Matrix is not square.");
  mtrx_destroy(matrix);
}


void test_mtrx_suite__mtrx_eq_dim(void) {
  matrix_t *A = mtrx_empty(6, 10);
  matrix_t *B = mtrx_empty(6, 10);
  cl_assert_(mtrx_eq_dim(A, B), "Matrices should have equal dimensions.");
  mtrx_destroy(A);
  mtrx_destroy(B);
}


void test_mtrx_suite__mtrx_not_eq_dim(void) {
  matrix_t *A = mtrx_empty(6, 10);
  matrix_t *B = mtrx_empty(10, 6);
  cl_assert_(!mtrx_eq_dim(A, B), "Matrices should not have equal dimensions.");
  mtrx_destroy(A);
  mtrx_destroy(B);
}


/*------------------------------ Max and Min --------------------------------*/

void test_mtrx_suite__mrtx_max(void) {
  matrix_t *matrix = mtrx_zeros(10, 10);
  matrix->values[5][5] = 3;
  cl_assert_(mtrx_max(matrix) == 3, "Max value incorrect.");
  mtrx_destroy(matrix);
}


void test_mtrx_suite__mtrx_min(void) {
  matrix_t *matrix = mtrx_ones(10, 10);
  matrix->values[9][9] = -5;
  cl_assert_(mtrx_min(matrix) == -5, "Min value incorrect.");
  mtrx_destroy(matrix);
}
