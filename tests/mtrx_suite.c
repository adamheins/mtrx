#include "../mtrx.h"
#include "clar.h"


void test_mtrx__init_id_matrix(void) {
  matrix_t *matrix = mtrx_id(10);
  for (size_t i = 0; i < 10; ++i)
    cl_assert_(matrix->values[i][i] == 1, "Values on diagonal should be one.");
  mtrx_destroy(matrix);
}

