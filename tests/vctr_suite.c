#include "../mtrx.h"
#include "clar.h"
#include <stdio.h>
#include <math.h>

bool vctr_eq_within_tol(vector_t *A, vector_t *B, scalar_t tol) {
  if (!vctr_eq_len(A, B))
    return NULL;
  for (size_t i = 0; i < A->length; ++i) {
    scalar_t diff = fabs(A->values[i] - B->values[i]);
    if (diff > tol)
      return false;
  }
  return true;
}


bool vctr_all_elements_equal(vector_t *vector, scalar_t value) {
  for (size_t i = 0 ; i < vector->length; ++i) {
    if (vector->values[i] != value)
      return false;
  }
  return true;
}


/*--------------------------- Initialization --------------------------------*/

void test_vctr_suite__vctr_init(void) {
  vector_t *master = vctr_empty(3);
  master->values[0] = 1.5;
  master->values[1] = -2;
  master->values[2] = 103.27;

  // Normal case.
  vector_t *vector = vctr_init("1.5 -2 103.27");
  cl_assert_(vctr_eq_within_tol(master, vector, 0.00001), "Vectors not equal.");
  vctr_destroy(vector);

  // Weird spacing case.
  vector = vctr_init("   1.5  -2 103.27  ");
  cl_assert_(vctr_eq_within_tol(master, vector, 0.00001), "Vectors not equal.");
  vctr_destroy(vector);

  vctr_destroy(master);
}


void test_vctr_suite__vctr_empty(void) {
  vector_t *vector = vctr_empty(5);
  cl_assert_(vector->length == 5, "Length error.");
  vctr_destroy(vector);
}


void test_vctr_suite__vctr_zeros(void) {
  vector_t *vector = vctr_zeros(5);
  cl_assert_(vctr_all_elements_equal(vector, 0),
      "All elements should be zero.");
  vctr_destroy(vector);
}


void test_vctr_suite__vctr_ones(void) {
  vector_t *vector = vctr_ones(5);
  cl_assert_(vctr_all_elements_equal(vector, 1),
      "All elements should be one.");
  vctr_destroy(vector);
}


void test_vctr_suite__vctr_copy(void) {
  vector_t *vector = vctr_rnd(10, 10);
  vector_t *copy = vctr_copy(vector);

  cl_assert_(copy->length == vector->length, "Lengths should be equal.");
  cl_assert_(vctr_eq(vector, copy), "Vectors should be equal.");

  vctr_destroy(vector);
  vctr_destroy(copy);
}


/*------------------------------ Comparisons --------------------------------*/

void test_vctr_suite__vctr_eq_len(void) {
  vector_t *A = vctr_empty(5);
  vector_t *B = vctr_ones(5);

  cl_assert_(A->length == B->length, "Vectors should have equal lengths.");

  vctr_destroy(A);
  vctr_destroy(B);
}


void test_vctr_suite__vctr_eq(void) {
  vector_t *vector = vctr_rnd(10, 10);
  vector_t *copy = vctr_copy(vector);

  cl_assert_(vctr_eq_within_tol(vector, copy, 0), "Vectors should be equal.");

  vctr_destroy(vector);
  vctr_destroy(copy);
}


/*------------------------------ Max and Min --------------------------------*/

void test_vctr_suite__vctr_max(void) {
  vector_t *vector = vctr_ones(10);
  vector->values[3] = 20;
  cl_assert_(vctr_max(vector) == 20, "Max value incorrect.");
  vctr_destroy(vector);
}


void test_vctr_suite__vctr_min(void) {
  vector_t *vector = vctr_ones(10);
  vector->values[3] = -5;
  cl_assert_(vctr_min(vector) == -5, "Min value incorrect.");
  vctr_destroy(vector);
}


/*------------------------------- Operations --------------------------------*/

void test_vctr_suite__vctr_dot_prod(void) {
  vector_t *A = vctr_init("1 2 3");
  vector_t *B = vctr_init("4 5 6");

  scalar_t result = vctr_dot_prod(A, B);
  cl_assert_(result == 32, "Dot product incorrect.");

  vctr_destroy(A);
  vctr_destroy(B);
}


void test_vctr_suite__vctr_cross_prod(void) {
  vector_t *A = vctr_init("3 -3 1");
  vector_t *B = vctr_init("4 9 2");

  vector_t *result = vctr_cross_prod(A, B);
  cl_assert_(result->values[0] == -15, "Cross product incorrect.");
  cl_assert_(result->values[1] == -2, "Cross product incorrect.");
  cl_assert_(result->values[2] == 39, "Cross product incorrect.");

  vctr_destroy(A);
  vctr_destroy(B);
  vctr_destroy(result);
}


void test_vctr_suite__vctr_mag(void) {
  vector_t *vector = vctr_init("1 2 3 5 5");
  scalar_t mag = vctr_mag(vector);
  cl_assert_(mag == 8,  "Magnitude incorrect.");
  vctr_destroy(vector);
}
