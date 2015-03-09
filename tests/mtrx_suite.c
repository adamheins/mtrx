#include "../mtrx.h"
#include "clar.h"
#include <stdio.h>
#include <math.h>



/*------------------------------- Helpers -----------------------------------*/

void assert_all_elements_equal(matrix_t *matrix, scalar_t value) {
  for (size_t i = 0; i < matrix->rows; ++i) {
    for (size_t j = 0; j < matrix->columns; ++j)
      cl_assert_(matrix->values[i][j] == value, "Value error!");
  }
}

bool mtrx_eq_within_tol(matrix_t *A, matrix_t *B, scalar_t tol) {
  if (!mtrx_eq_dim(A, B))
    return false;
  for (size_t i = 0; i < A->rows; ++i) {
    for (size_t j = 0; j < A->columns; ++j) {
      scalar_t diff = fabs(A->values[i][j] - B->values[i][j]);
      if (diff > tol)
        return false;
    }
  }
  return true;
}


/*--------------------------- Initialization --------------------------------*/

void test_mtrx_suite__mtrx_init(void) {
  matrix_t *master = mtrx_empty(2, 2);
  master->values[0][0] = 1.5;
  master->values[0][1] = 2;
  master->values[1][0] = -3.2;
  master->values[1][1] = 4.365;

  // Normal case.
  matrix_t *matrix = mtrx_init("1.5 2; -3.2 4.365");
  cl_assert_(mtrx_eq_within_tol(master, matrix, 0.00001), "Matrices inequal!");
  mtrx_destroy(matrix);

  // Weird spacing case.
  matrix = mtrx_init("  1.5 2  ; -3.2    4.365     ");
  cl_assert_(mtrx_eq_within_tol(master, matrix, 0.00001), "Matrices inequal!");
  mtrx_destroy(matrix);

  mtrx_destroy(master);
}

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


void text_mtrx_suite__mtrx_mult_vctr(void) {
  matrix_t *matrix = mtrx_rnd(50, 50, 100);
  vector_t *vector = vctr_rnd(50, 100);
  vector_t *result = mtrx_mult_vctr(matrix, vector);

	for (size_t i = 0; i < matrix->rows; ++i) {
    scalar_t value = 0;
		for (size_t j = 0; j < matrix->columns; ++j)
			value += matrix->values[i][j] * vector->values[j];
    cl_assert_(result->values[i] == value, "Multiplication error!");
	}
  mtrx_destroy(matrix);
  vctr_destroy(vector);
  vctr_destroy(result);
}


/*------------------------- Pointwise Arithmetic ----------------------------*/

void test_mtrx_suite__mtrx_pw_mult(void) {
  matrix_t *matrix = mtrx_empty(2, 2);
  matrix->values[0][0] = 1;
  matrix->values[0][1] = 2;
  matrix->values[1][0] = 3;
  matrix->values[1][1] = 4;

  matrix_t *result = mtrx_pw_mult(matrix, matrix);
  cl_assert_(result->values[0][0] == 1, "PW multiplication error!");
  cl_assert_(result->values[0][1] == 4, "PW multiplication error!");
  cl_assert_(result->values[1][0] == 9, "PW multiplication error!");
  cl_assert_(result->values[1][1] == 16, "PW multiplication error!");

  mtrx_destroy(matrix);
  mtrx_destroy(result);
}


void test_mtrx_suite__mtrx_pw_div(void) {
  matrix_t *A = mtrx_empty(2, 2);
  A->values[0][0] = 10;
  A->values[0][1] = 20;
  A->values[1][0] = 30;
  A->values[1][1] = 40;

  matrix_t *B = mtrx_empty(2, 2);
  B->values[0][0] = 5;
  B->values[0][1] = 2;
  B->values[1][0] = 3;
  B->values[1][1] = 5;

  matrix_t *result = mtrx_pw_div(A, B);
  cl_assert_(result->values[0][0] == 2, "PW division error!");
  cl_assert_(result->values[0][1] == 10, "PW division error!");
  cl_assert_(result->values[1][0] == 10, "PW division error!");
  cl_assert_(result->values[1][1] == 8, "PW division error!");

  mtrx_destroy(A);
  mtrx_destroy(B);
  mtrx_destroy(result);
}


void test_mtrx_suite__mtrx_pow(void) {
  matrix_t *matrix = mtrx_empty(2, 2);
  matrix->values[0][0] = 1;
  matrix->values[0][1] = 2;
  matrix->values[1][0] = 3;
  matrix->values[1][1] = 4;

  matrix_t *result = mtrx_pw_pow(matrix, matrix);
  cl_assert_(result->values[0][0] == 1, "PW pow error!");
  cl_assert_(result->values[0][1] == 4, "PW pow error!");
  cl_assert_(result->values[1][0] == 27, "PW pow error!");
  cl_assert_(result->values[1][1] == 256, "PW pow error!");

  mtrx_destroy(matrix);
  mtrx_destroy(result);
}


/*-------------------------- Matrix Manipulation ----------------------------*/

void test_mtrx_suite__mtrx_row_swap(void) {
  matrix_t *original = mtrx_rnd(10, 10, 100);
  matrix_t *copy = mtrx_copy(original);

  mtrx_row_swap(copy, 4, 5);

  for (size_t i = 0; i < original->columns; ++i) {
    cl_assert_(copy->values[4][i] == original->values[5][i],
        "Row swap error!");
    cl_assert_(copy->values[5][i] == original->values[4][i],
        "Row swap error!");
  }
  mtrx_destroy(original);
  mtrx_destroy(copy);
}


void test_mtrx_suite__mtrx_col_swap(void) {
  matrix_t *original = mtrx_rnd(10, 10, 100);
  matrix_t *copy = mtrx_copy(original);

  mtrx_col_swap(copy, 4, 5);

  for (size_t i = 0; i < original->rows; ++i) {
    cl_assert_(copy->values[i][4] == original->values[i][5],
        "Column swap error!");
    cl_assert_(copy->values[i][5] == original->values[i][4],
        "Column swap error!");
  }
  mtrx_destroy(original);
  mtrx_destroy(copy);
}


void test_mtrx_suite__mtrx_scale_row(void) {
  matrix_t *original = mtrx_rnd(10, 10, 100);
  matrix_t *copy = mtrx_copy(original);

  mtrx_scale_row(copy, 4, 2.0);

  for (size_t i = 0; i < original->columns; ++i)
    cl_assert_(copy->values[4][i] == original->values[4][i] * 2.0,
        "Row scaling error!");

  mtrx_destroy(original);
  mtrx_destroy(copy);
}


void test_mtrx_suite__mtrx_scale_col(void) {
  matrix_t *original = mtrx_rnd(10, 10, 100);
  matrix_t *copy = mtrx_copy(original);

  mtrx_scale_col(copy, 4, 2.0);

  for (size_t i = 0; i < original->rows; ++i)
    cl_assert_(copy->values[i][4] == original->values[i][4] * 2.0,
        "Column scaling error!");

  mtrx_destroy(original);
  mtrx_destroy(copy);
}


/*----------------- Row and Column Accessors and Mutators -------------------*/

void text_mtrx_suite__mtrx_get_row(void) {
  matrix_t *matrix = mtrx_rnd(10, 10, 100);
  vector_t *row = mtrx_get_row(matrix, 2);

  for (size_t i = 0; i < row->length; ++i)
    cl_assert_(row->values[i] == matrix->values[2][i], "get_row error!");

  mtrx_destroy(matrix);
  vctr_destroy(row);
}


void text_mtrx_suite__mtrx_get_col(void) {
  matrix_t *matrix = mtrx_rnd(10, 10, 100);
  vector_t *col = mtrx_get_col(matrix, 2);

  for (size_t i = 0; i < col->length; ++i)
    cl_assert_(col->values[i] == matrix->values[i][2], "get_col error!");

  mtrx_destroy(matrix);
  vctr_destroy(col);
}


void text_mtrx_suite__mtrx_set_row(void) {
  matrix_t *matrix = mtrx_rnd(10, 10, 100);
  vector_t *row = vctr_rnd(10, 100);

  mtrx_set_row(matrix, row, 2);

  for (size_t i = 0; i < row->length; ++i)
    cl_assert_(matrix->values[2][i] == row->values[i], "set_row error!");

  mtrx_destroy(matrix);
  vctr_destroy(row);
}


void text_mtrx_suite__mtrx_set_col(void) {
  matrix_t *matrix = mtrx_rnd(10, 10, 100);
  vector_t *col = vctr_rnd(10, 100);

  mtrx_set_col(matrix, col, 2);

  for (size_t i = 0; i < col->length; ++i)
    cl_assert_(matrix->values[i][2] == col->values[i], "set_col error!");

  mtrx_destroy(matrix);
  vctr_destroy(col);
}


/*----------------------------- Sub-matrices --------------------------------*/

void test_mtrx_suite__mtrx_sub_matrix(void) {
  matrix_t *matrix = mtrx_rnd(4, 4, 100);

  indexer_t *indexer = indexer_init(2);
  indexer->values[0] = 1;
  indexer->values[1] = 3;

  matrix_t *sub_matrix = mtrx_sub_matrix(matrix, indexer, indexer);

  cl_assert_(sub_matrix->values[0][0] == matrix->values[1][1],
      "Sub-matrix error!");
  cl_assert_(sub_matrix->values[0][1] == matrix->values[1][3],
      "Sub-matrix error!");
  cl_assert_(sub_matrix->values[1][0] == matrix->values[3][1],
      "Sub-matrix error!");
  cl_assert_(sub_matrix->values[1][1] == matrix->values[3][3],
      "Sub-matrix error!");

  mtrx_destroy(matrix);
  mtrx_destroy(sub_matrix);
}


void test_mtrx_suite__mtrx_sub_block(void) {
  matrix_t *matrix = mtrx_rnd(10, 10, 100);
  matrix_t *sub_block = mtrx_sub_block(matrix, 1, 4, 1, 4);

  for (size_t i = 0; i < sub_block->rows; ++i) {
    for (size_t j = 0; j < sub_block->columns; ++j)
      cl_assert_(sub_block->values[i][j] == matrix->values[i + 1][j + 1],
          "sub-block error!");
  }

  mtrx_destroy(matrix);
  mtrx_destroy(sub_block);
}


/*------------------------------- Operations --------------------------------*/

void test_mtrx_suite__mtrx_transpose(void) {
  matrix_t *matrix = mtrx_rnd(10, 10, 100);
  matrix_t *transpose = mtrx_transpose(matrix);

	for (size_t i = 0; i < matrix->rows; ++i) {
		for (size_t j = 0; j < matrix->columns; ++j)
			cl_assert_(transpose->values[j][i] == matrix->values[i][j],
        "Transpose error!");
	}

  mtrx_destroy(matrix);
  mtrx_destroy(transpose);
}


void test_mtrx_suite__mtrx_det(void) {
  matrix_t *matrix = mtrx_empty(3, 3);
  matrix->values[0][0] = 6;
  matrix->values[0][1] = 1;
  matrix->values[0][2] = 1;

  matrix->values[1][0] = 4;
  matrix->values[1][1] = -2;
  matrix->values[1][2] = 5;

  matrix->values[2][0] = 2;
  matrix->values[2][1] = 8;
  matrix->values[2][2] = 7;

  scalar_t det = mtrx_det(matrix);

  // Do to the inaccuracies of floating point values, this result may not be
  // exact. Thus the determinant is rounded to an integer here.
  cl_assert_(round(det) == -306, "Determinant error!");

  mtrx_destroy(matrix);
}


void test_mtrx_suite__mtrx_inv(void) {
  matrix_t *matrix = mtrx_empty(3, 3);
  matrix->values[0][0] = 1;
  matrix->values[0][1] = 3;
  matrix->values[0][2] = 1;

  matrix->values[1][0] = 1;
  matrix->values[1][1] = 1;
  matrix->values[1][2] = 2;

  matrix->values[2][0] = 2;
  matrix->values[2][1] = 3;
  matrix->values[2][2] = 4;

  matrix_t *inv = mtrx_inv(matrix);
  cl_assert_(inv->values[0][0] == 2, "Inverse error!");
  cl_assert_(inv->values[0][1] == 9, "Inverse error!");
  cl_assert_(inv->values[0][2] == -5, "Inverse error!");

  cl_assert_(inv->values[1][0] == 0, "Inverse error!");
  cl_assert_(inv->values[1][1] == -2, "Inverse error!");
  cl_assert_(inv->values[1][2] == 1, "Inverse error!");

  cl_assert_(inv->values[2][0] == -1, "Inverse error!");
  cl_assert_(inv->values[2][1] == -3, "Inverse error!");
  cl_assert_(inv->values[2][2] == 2, "Inverse error!");

  mtrx_destroy(matrix);
  mtrx_destroy(inv);
}


/*----------------------------- System Solving ------------------------------*/

void test_mtrx_suite__mtrx_solve(void) {
  matrix_t *matrix = mtrx_rnd(10, 10, 100);
  vector_t *vector = vctr_rnd(10, 100);

  vector_t *sol = mtrx_solve(matrix, vector);
  vector_t *compare = mtrx_mult_vctr(matrix, sol);

  for (size_t i = 0; i < vector->length; ++i)
    cl_assert_(round(compare->values[i]) == round(vector->values[i]),
        "System solved incorrectly!");

  mtrx_destroy(matrix);
  vctr_destroy(vector);
  vctr_destroy(sol);
  vctr_destroy(compare);
}


void test_mtrx_suite__mtrx_is_diag_dom(void) {
  matrix_t *matrix = mtrx_init("3 -2 1; 1 -3 2; -1 2 4");
  cl_assert_(mtrx_is_diag_dom(matrix), "Matrix should be diag dom.");
  mtrx_destroy(matrix);

  matrix = mtrx_init("-2 2 1; 1 3 2; 1 -2 0");
  cl_assert_(!mtrx_is_diag_dom(matrix), "Matrix should not be diag dom.");
  mtrx_destroy(matrix);
}


void test_mtrx_suite__mtrx_make_diag_dom(void) {
  matrix_t *matrix = mtrx_init("1 -3 2; 3 -2 1; -1 2 4");
  cl_assert_(!mtrx_is_diag_dom(matrix), "Matrix should not be diag dom.");
  matrix_t *diag_dom_matrix = mtrx_make_diag_dom(matrix);
  cl_assert_(mtrx_is_diag_dom(diag_dom_matrix), "Matrix should be diag dom.");
  mtrx_destroy(matrix);
  mtrx_destroy(diag_dom_matrix);
}


void test_mtrx_suite__mtrx_solve_gs(void) {
  matrix_t *A = mtrx_init("4 -1 -1; -2 6 1; -1 1 7");
  vector_t *B = vctr_init("3 6 -9");
  vector_t *X = vctr_zeros(3);

  mtrx_solve_gs(A, B, X, 0.000001);

  vector_t *compare = mtrx_mult_vctr(A, X);

  for (size_t i = 0; i < compare->length; ++i) {
    cl_assert_(round(compare->values[i]) == round(B->values[i]),
        "GS method error.");
  }

  mtrx_destroy(A);
  vctr_destroy(B);
  vctr_destroy(X);
  vctr_destroy(compare);
}
