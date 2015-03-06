#include "../mtrx.h"
#include "clar.h"
#include <stdio.h>



void assert_all_elements_equal(matrix_t *matrix, scalar_t value) {
  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j)
      cl_assert_(matrix->values[i][j] == value, "Value error!");
  }
}


/*--------------------------- Initialization --------------------------------*/

void test_mtrx_suite__mtrx_empty(void) {
  matrix_t *matrix = mtrx_empty(5, 10);
  cl_assert_(matrix->rows == 5, "Rows not assigned properly!");
  cl_assert_(matrix->columns == 10, "Columns not assigned properly.");
  mtrx_destroy(matrix);
}

void test_mtrx_suite__mtrx_zeros(void) {
  matrix_t *matrix = mtrx_zeros(5, 5);
  assert_all_elements_equal(matrix, 0);
  mtrx_destroy(matrix);
}


void test_mtrx_suite__mtrx_ones(void) {
  matrix_t *matrix = mtrx_ones(5, 5);
  assert_all_elements_equal(matrix, 1);
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


/*--------------------------- Matrix Arithmetic -----------------------------*/

void test_mtrx_suite__mtrx_add_ones(void) {
  matrix_t *A = mtrx_ones(10, 10);
  matrix_t *B = mtrx_ones(10, 10);
  matrix_t *C = mtrx_add(A, B);
  assert_all_elements_equal(C, 2);
  mtrx_destroy(A);
  mtrx_destroy(B);
  mtrx_destroy(C);
}

void test_mtrx_suite__mtrx_add_rand(void) {
  matrix_t *A = mtrx_rnd(10, 10, 20);
  matrix_t *B = mtrx_rnd(10, 10, 20);
  matrix_t *C = mtrx_add(A, B);
  for (size_t i = 0; i < A->rows; ++i) {
    for (size_t j = 0; j < A->columns; ++j)
      cl_assert_(C->values[i][j] == A->values[i][j] + B->values[i][j],
          "Addition error!");
  }
  mtrx_destroy(A);
  mtrx_destroy(B);
  mtrx_destroy(C);
}

void test_mtrx_suite__mtrx_add_manual(void) {
  matrix_t *A = mtrx_empty(3, 2);
  A->values[0][0] = -3;
  A->values[0][1] = -2;
  A->values[1][0] = -1;
  A->values[1][1] = 0;
  A->values[2][0] = 1;
  A->values[2][1] = 2;

  matrix_t *B = mtrx_empty(3, 2);
  B->values[0][0] = -5;
  B->values[0][1] = 6;
  B->values[1][0] = -1;
  B->values[1][1] = 3;
  B->values[2][0] = 2;
  B->values[2][1] = -1;

  matrix_t *C = mtrx_add(A, B);
  cl_assert_(C->values[0][0] == -8, "Addition error!");
  cl_assert_(C->values[0][1] == 4, "Addition error!");
  cl_assert_(C->values[1][0] == -2, "Addition error!");
  cl_assert_(C->values[1][1] == 3, "Addition error!");
  cl_assert_(C->values[2][0] == 3, "Addition error!");
  cl_assert_(C->values[2][1] == 1, "Addition error!");

  mtrx_destroy(A);
  mtrx_destroy(B);
  mtrx_destroy(C);
}

void test_mtrx_suite__mtrx_add_null(void) {
  matrix_t *A = mtrx_empty(4, 2);
  matrix_t *B = mtrx_empty(2, 4);
  matrix_t *C = mtrx_add(A, B);
  cl_assert_(C == NULL, "Cannot add matrices with unequal dimensions.");
  mtrx_destroy(A);
  mtrx_destroy(B);
}

void test_mtrx_suite__mtrx_subtract_ones(void) {
  matrix_t *A = mtrx_ones(10, 10);
  matrix_t *B = mtrx_ones(10, 10);
  matrix_t *C = mtrx_subtract(A, B);
  assert_all_elements_equal(C, 0);
  mtrx_destroy(A);
  mtrx_destroy(B);
  mtrx_destroy(C);
}

void test_mtrx_suite__mtrx_subtract_rand(void) {
  matrix_t *A = mtrx_rnd(10, 10, 20);
  matrix_t *B = mtrx_rnd(10, 10, 20);
  matrix_t *C = mtrx_subtract(A, B);
  for (size_t i = 0; i < A->rows; ++i) {
    for (size_t j = 0; j < A->columns; ++j)
      cl_assert_(C->values[i][j] == A->values[i][j] - B->values[i][j],
          "Subtraction error!");
  }
  mtrx_destroy(A);
  mtrx_destroy(B);
  mtrx_destroy(C);
}

void test_mtrx_suite__mtrx_subtract_manual(void) {
  matrix_t *A = mtrx_empty(3, 2);
  A->values[0][0] = -3;
  A->values[0][1] = -2;
  A->values[1][0] = -1;
  A->values[1][1] = 0;
  A->values[2][0] = 1;
  A->values[2][1] = 2;

  matrix_t *B = mtrx_empty(3, 2);
  B->values[0][0] = -5;
  B->values[0][1] = 6;
  B->values[1][0] = -1;
  B->values[1][1] = 3;
  B->values[2][0] = 2;
  B->values[2][1] = -1;

  matrix_t *C = mtrx_subtract(A, B);
  cl_assert_(C->values[0][0] == 2, "Subtract error!");
  cl_assert_(C->values[0][1] == -8, "Subtract error!");
  cl_assert_(C->values[1][0] == 0, "Subtract error!");
  cl_assert_(C->values[1][1] == -3, "Subtract error!");
  cl_assert_(C->values[2][0] == -1, "Subtract error!");
  cl_assert_(C->values[2][1] == 3, "Subtract error!");

  mtrx_destroy(A);
  mtrx_destroy(B);
  mtrx_destroy(C);
}

void test_mtrx_suite__mtrx_subtract_null(void) {
  matrix_t *A = mtrx_empty(4, 2);
  matrix_t *B = mtrx_empty(2, 4);
  matrix_t *C = mtrx_subtract(A, B);
  cl_assert_(C == NULL, "Cannot subtract matrices with unequal dimensions.");
  mtrx_destroy(A);
  mtrx_destroy(B);
}

void test_mtrx_suite__mtrx_scale_ones(void) {
  matrix_t *A = mtrx_ones(10, 10);
  matrix_t *B = mtrx_scale(A, 100);
  assert_all_elements_equal(B, 100);
  mtrx_destroy(A);
  mtrx_destroy(B);
}

void test_mtrx_suite__mtrx_scale_zeros(void) {
  matrix_t *A = mtrx_zeros(10, 10);
  matrix_t *B = mtrx_scale(A, 100);
  assert_all_elements_equal(B, 0);
  mtrx_destroy(A);
  mtrx_destroy(B);
}

void test_mtrx_suite__mtrx_mult_simple(void) {
  matrix_t *A = mtrx_empty(3, 3);
  A->values[0][0] = 1;
  A->values[0][1] = 2;
  A->values[0][2] = 3;
  A->values[1][0] = 4;
  A->values[1][1] = 5;
  A->values[1][2] = 6;
  A->values[2][0] = 7;
  A->values[2][1] = 8;
  A->values[2][2] = 9;

  matrix_t *B = mtrx_mult(A, A);
  cl_assert_(B->values[0][0] == 30, "Multiplication error!");
  cl_assert_(B->values[0][1] == 36, "Multiplication error!");
  cl_assert_(B->values[0][2] == 42, "Multiplication error!");
  cl_assert_(B->values[1][0] == 66, "Multiplication error!");
  cl_assert_(B->values[1][1] == 81, "Multiplication error!");
  cl_assert_(B->values[1][2] == 96, "Multiplication error!");
  cl_assert_(B->values[2][0] == 102, "Multiplication error!");
  cl_assert_(B->values[2][1] == 126, "Multiplication error!");
  cl_assert_(B->values[2][2] == 150, "Multiplication error!");

  mtrx_destroy(A);
  mtrx_destroy(B);
}

void test_mtrx_suite__mtrx_mult_complex(void) {
  matrix_t *A = mtrx_rnd(500, 500, 100);
  matrix_t *B = mtrx_rnd(500, 500, 100);
  matrix_t *C = mtrx_mult(A, B);

	for (size_t i = 0; i < B->columns; ++i) {
		for (size_t j = 0; j < A->rows; ++j) {
      scalar_t value = 0;
			for (size_t k = 0; k < A->columns; ++k)
				value += A->values[j][k] * B->values[k][i];
      cl_assert_(C->values[j][i] == value, "Multiplication error!");
		}
	}
  mtrx_destroy(A);
  mtrx_destroy(B);
  mtrx_destroy(C);
}
